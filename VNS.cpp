#include "HeuristicDecoder.h"

using namespace std;

const int K_MAX = 50;
HeuristicDecoder *heu;
double maximumTime;
clock_t startTime;
int channels20MHz[25] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24};
int arrRoot[] = {24, 41, 42, 43, 44};
double chanThroughput[45];
bool inSolution[45];

bool isStoppingCriteriaReached() {
  return (((double) (clock() - startTime)) / CLOCKS_PER_SEC) >= maximumTime;
}

void addDrop(Solution &current, Link link) {
  vector<int> usedChannels(current.getScheduledChannels());

  Solution aux(current);
  bool flag = true;

  for (int ch : usedChannels) {
    Solution copy(current);
    Link newLink(ch);
    copy.insert(newLink);

    if ((copy > aux) || flag) {
      flag = false;
      aux = copy;
    }
  }

  current = aux;
}

void betaAddDrop(Solution &current, int beta = 1) {
  vector<int> zeroLinks = current.getZeroLinks();
  beta = min(beta, int(zeroLinks.size()));
  for (int k = 0; k < beta; k++) {
    int rndIndex = rng.randInt(zeroLinks.size() - 1);
    int idConnection = zeroLinks[rndIndex];

    Link toInsert(idConnection);
    addDrop(current, toInsert);

    swap(zeroLinks[rndIndex], zeroLinks.back());
    zeroLinks.pop_back();
  }

  current.setZeroLinks(zeroLinks);
}

void betaReinsert(Solution &current, int beta = 1) {
  deque<Link> scheduledLinks(current.getScheduledLinks());
  beta = min(beta, int(scheduledLinks.size()));

  for (int i = 0; i < beta; i++) {
    int rndIndex = rng.randInt(scheduledLinks.size() - 1);

    Link aux = scheduledLinks[rndIndex];
    assert(current.removeLink(aux));

    addDrop(current, Link(aux.ch));
  }

  current.setScheduledLinks(scheduledLinks);

  //A versao de Mauricio tem um addDrop() ao final de todo Reinsert. Verificar.
}

vector<Channel> buildMultipleSolution(const Solution &curr, int rootChannel) {
  vector<Link> slinks;
  vector<Channel> ret(45);

  for (const Link &x : curr.getScheduledLinks()) {
    if (overlap[x.ch][rootChannel]) {
      slinks.emplace_back(x);
    }
  }

  sort(slinks.begin(), slinks.end());

  int idx = 0;
  while (idx < slinks.size()) {
    int a = idx;
    int b = upper_bound(slinks.begin() + a, slinks.end(), slinks[a]) - slinks.begin();

    if (b < slinks.size()) {
      int i = b;
      if (father[slinks[a].ch] == father[slinks[i].ch]) {
        int j = upper_bound(slinks.begin() + i, slinks.end(), slinks[i]) - slinks.begin();

        for (int k = a; k < j; k++) {
          Link aux(slinks[k]);
          aux.setChannel(father[slinks[k].ch]);

          slinks.emplace_back(aux);
        }

        idx = j;
      } else {
        idx = i;
      }
    } else {
      idx = b;
    }
  }

  for (int i = 0; i < slinks.size(); i++) {
    ret[slinks[i].ch].links.emplace_back(slinks[i].id);
    ret[slinks[i].ch].interference += slinks[i].interference;
  }

  return ret;
}

double solve(vector<Channel> &refCh, int chId) {
  double ret = refCh[chId].throughput;
  inSolution[chId] = true;

  if (refCh[chId].childs[0] != -1 && refCh[chId].childs[1] != -1) {
    double candidate = solve(refCh, refCh[chId].childs[0]) + solve(refCh, refCh[chId].childs[1]);

    if (candidate > refCh[chId].throughput) {
      ret = candidate;
      inSolution[chId] = false;
    }
  }

  return ret;
}

void recoverSolutionFromDP(vector<Channel> &sChannel, int chId, bool clean) {
  if (clean)
    inSolution[chId] = false;

  if (inSolution[chId])
    clean = true;

  if (sChannel[chId].childs[0] != -1 && sChannel[chId].childs[1] != -1) {
    recoverSolutionFromDP(sChannel, sChannel[chId].childs[0], clean);
    recoverSolutionFromDP(sChannel, sChannel[chId].childs[1], clean);
  }
}

double computeLinkThroughput(int idConnection, double localInterference, double bw) {
  double throughput = 0.0;

  double sinr;
  if (localInterference == 0.0) {
    sinr = 1e9;
  } else {
    sinr = (powerSender / pow(distanceMatrix[idConnection][idConnection], alfa)) / (localInterference + noise);
  }

  int mxDataRate = (bw == 20) ? 9 : 8;
  int i = 0;

  while (i < mxDataRate && sinr > SINR[i++][bw]);

  throughput = dataRates[i][bw];
  return throughput;
}

void computeChannelsThroughput(vector<Channel> &usedChannels) {
  for (Channel &ch : usedChannels) {
    ch.throughput = 0.0;

    for (int i = 0; i < ch.links.size(); i++) {
      int idConnection = ch.links[i];
      double localInterference = 0.0;

      for (int j = 0; j < ch.links.size(); j++) {
        int idOtherConnection = ch.links[j];
        localInterference += interferenceMatrix[idOtherConnection][idConnection];
      }

      ch.throughput += computeLinkThroughput(idConnection, localInterference, whichBw(ch.id));
    }
  }
}

double calcDP(vector<Channel> &chan) {
  double OF = 0.0;
  for (const int x : arrRoot) {
    OF += solve(chan, x);
  }
  return OF;
}

Solution reconstructSolution(const Solution &S, vector<Channel> &ret, double ansDP) {
  for (int x : arrRoot)
    recoverSolutionFromDP(ret, x, false);

  Solution final_;

  for (const Channel &chan : ret) {
    if (inSolution[chan.id]) {

      for (int j = 0; j < chan.links.size(); j++) {
        //final_.insert() TODO: inserir link chan.link[j] em final_
      }
    }
  }

  assert(final_.getNumberOfScheduledLinks() <= heu->getQuantConnections());
  final_.computeObjective(false);

  assert(ansDP == final_.getObjective());
}

vector<Channel> initDP(const Solution &S) {
  vector<Channel> usedChannels;
  for (int i = 0; i < 5; i++) {
    vector<Channel> aux = buildMultipleSolution(S, arrRoot[i]);
    usedChannels.insert(usedChannels.end(), aux.begin(), aux.end());
  }

  computeChannelsThroughput(usedChannels);
  return usedChannels;
}

void deleteFromChannel(const int k, Channel &channel, const Solution &curr) {
  for (int i = 0; i < channel.links.size(); i++) {
    if (channel.links[i] == k) {
      swap(channel.links[i], channel.links.back());
      break;
    }
  }

  channel.links.pop_back();
  computeChannelsThroughput({channel});
}

void insertInChannel(const int k, Channel &channel, const Solution &curr) {
  channel.links.push_back(k);
  computeChannelsThroughput({channel});
}

void setDP(vector<Channel> &rep) {
  memset(chanThroughput, 0, sizeof chanThroughput);
  memset(inSolution, 0, sizeof inSolution);

  for (int i = 0; i < rep.size(); i++) {
    chanThroughput[rep[i].id] = rep[i].throughput;
  }
}

//currSol must be a \delta solution, and multiple is the solution returned by buildMultipleSolution()
void VNS_Reinsert(const Solution &curr) {
  vector<Channel> multiple;
  for (int i = 0; i < 5; i++) {
    vector<Channel> aux = buildMultipleSolution(curr, arrRoot[i]);
    multiple.insert(multiple.end(), aux.begin(), aux.end());
  }
  computeChannelsThroughput(multiple);

  bool improved = false;
  do {
    improved = false;
    for (int k = 0; k < heu->getQuantConnections(); k++) {
      vector<Channel> multipleClean(multiple);

      for (int i = 0; i < multipleClean.size(); i++) {
        for (int j = 0; j < multipleClean[i].links.size(); j++) {
          if (multipleClean[i].links[j] == k) {
            //TODO: remover link desse canal
            deleteFromChannel(k, multipleClean[i], curr);
            break;
          }
        }
      }

      vector<Channel> multipleContaining(multipleClean);

      for (Channel &ch_ : multipleContaining) {
        //TODO: inserir link 'i' em todos os canais ch_
        insertInChannel(k, ch_, curr);
      }

      double bestFO = -1;
      short bestChannel = -1;

      vector<int> scheduledChannels = curr.getScheduledChannels();

      for (const int aux : scheduledChannels) {
        setDP(multipleClean);

        int currChan = aux;

        while (currChan != -1) {
          chanThroughput[currChan] = multipleContaining[currChan].throughput;
          currChan = father[currChan];
        }

        double FO = calcDP(multiple);
        if (FO > bestFO) {
          bestFO = FO;
          bestChannel = aux;
        }
      }

      if (bestFO > curr.getObjective(false)) { //TODO: curr?
        improved = true;

        for (int i = 0; i < multiple.size(); i++) {
          for (int j = 0; j < multiple[i].links.size(); j++) {
            if (multiple[i].links[j] == k) {
              deleteFromChannel(k, multiple[i], curr);
            }
          }
        }

        while (bestChannel != -1) {
          insertInChannel(k, multiple[bestChannel], curr);
          bestChannel = father[bestChannel];
        }
      }
    }
  } while (improved);

  reconstructSolution(multiple);
}

bool allChannels20MHz(const Solution &check) {
  for (const int &x : check.getScheduledChannels()) {
    if (whichBw(x) > 20)
      return false;
  }

  return true;
}

Solution convert(const Solution &aux) {
  Solution ret(aux);

  while (!allChannels20MHz(ret)) {
    for (Link &x : ret.getScheduledLinks()) {
      if (whichBw(x.getChannel()) > 20) {
        split(ret, ret, x.getChannel(), true);
        break;
      }
    }
  }
  return ret;
}

Solution localSearch(Solution &current) {
  current = convert(current);
  vector<Channel> chDP = initDP(current);
  double dpOF = calcDP(chDP);

  Solution localOptima = reconstructSolution(current, chDP, dpOF);
  return localOptima;
}

Solution pertubation(Solution S, int k, const int NUMBER_OF_LINKS) {
  double threeshold = rng.rand();
  if (threeshold >= .5 && (S.getNumberOfScheduledLinks() < heu->getQuantConnections())) {
    betaAddDrop(S, lround(NUMBER_OF_LINKS * k * 0.01));
  } else {
    betaReinsert(S, lround(NUMBER_OF_LINKS * k * 0.01));
  }
  return S;
}

void init(const string &openingFile = "", double timeLimit = 10) {
#ifdef DEBUG_CLION
  puts("WITH DEBUG");
  freopen("/Users/jjaneto/Downloads/codes_new/BRKGA_FF_Best/Instancias/D250x250/U_8/U_8_1.txt", "r", stdin);
#else
  if (!openingFile.empty()) {
    fprintf(stderr, "trying to open input file %s\n", openingFile.c_str());
    freopen(openingFile.c_str(), "r", stdin);

    if (stdin == NULL) {
      fprintf(stderr, "error opening input file\n");
      exit(-1);
    }
  }
#endif

  heu = new HeuristicDecoder();
  maximumTime = timeLimit;

  fprintf(stdout, "will execute for %lf seconds\n", maximumTime);
}

int main(int argc, char *argv[]) {
  if (argc != 3) {
    fprintf(stderr, "wrong arguments\n");
    exit(-1);
  }

  FILE *solutionFile = fopen(argv[2], "a");

  if (solutionFile == NULL) {
    fprintf(stderr, "error opening solutionFile file\n");
    exit(-1);
  }

  init(argv[1], 10);
  const int NUMBER_OF_LINKS = heu->getQuantConnections();
  startTime = clock();

  //--------------------------------------------------------

  int K_MUL = max(1, int(NUMBER_OF_LINKS / 100));

  Solution curr = heu->generateSolution();
  Solution S_global = localSearch(curr);
  Solution localMax = curr;
  while (!isStoppingCriteriaReached()) {
    int k = 1;
    while (k <= K_MAX && !isStoppingCriteriaReached()) {
      curr = localMax;

      if (rng.randInt(1)) {
        //AddDrop
        betaAddDrop(curr, k * K_MUL);
      } else {
        //Reinsert
        betaReinsert(curr, k * K_MUL);
      }

      Solution localAux = localSearch(curr);

      if (localAux > localMax) {
        localMax = curr;
        k = 1;
      } else {
        k++;
      }

      if (localMax > S_global) {
        S_global = localMax;
      }
    }
  }

  printf("%lf\n", S_global.getObjective());

  fprintf(solutionFile, "OBJECTIVE %lf\n", S_global.getObjective());
  for (const Link &x : S_global.getScheduledLinks()) {
    fprintf(solutionFile, "%d %d\n", x.id, x.ch);
  }

  fclose(solutionFile);
  delete heu;
  return 0;
}