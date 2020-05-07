#include "HeuristicDecoder.h"

using namespace std;

const int MAX_SPECTRUM = 3;
const int MAX_CHANNELS = 45;
double maximumTime;
clock_t startTime;
int parent[MAX_SPECTRUM][MAX_CHANNELS];
int child[MAX_SPECTRUM][MAX_CHANNELS][2];
int chanThroughput[MAX_SPECTRUM][MAX_CHANNELS];
bool inSolution[MAX_SPECTRUM][MAX_CHANNELS];

bool stop() {
  return (((double) (clock() - startTime)) / CLOCKS_PER_SEC) >= maximumTime;
}

bool allChannels20MHz(const Solution &sol) {
  for (const Spectrum &spectrum : sol.spectrums) {
    for (const Channel &channel : spectrum.channels) {
      if (channel.bandwidth != 20) {
        return false;
      }
    }
  }
  return true;
}

Solution convert20MHz(Solution sol) {
  for (Spectrum &spectrum : sol.spectrums) {
    sort(spectrum.channels.begin(), spectrum.channels.end());
  }

  while (!allChannels20MHz(sol)) {
    for (Spectrum &spectrum : sol.spectrums) {
      for (int c = 0; c < spectrum.channels.size(); c++) {
        if (spectrum.channels[c].bandwidth == 20)
          continue;

        int newBw = spectrum.channels[c].bandwidth / 2;
        Channel child1(newBw), child2(newBw);

        vector<Connection> connections = spectrum.channels[c].connections;

        for (Connection &conn : connections) {
          if (rng.randInt(1)) {
            child1.connections.emplace_back(conn);
          } else {
            child2.connections.emplace_back(conn);
          }
        }

        spectrum.channels[c] = child1;
        spectrum.channels.insert(spectrum.channels.begin() + c, child2); //TODO: it was .end(). Did it right?
      }
    }
  }

  return sol;
}

void setDP(const Solution &sol) { //TODO
  memset(chanThroughput, 0, sizeof chanThroughput);
  memset(inSolution, false, sizeof inSolution);

  for (int s = 0; s < sol.spectrums.size(); s++) {
    for (int c = 0; c < sol.spectrums[s].channels.size(); c++) {
      //TODO: verify if at this point, the throughput is already updated
      chanThroughput[s][c] = sol.spectrums[s].channels[c].throughput;
    }
  }
}

double solve(int i, int j) {
  inSolution[i][j] = true;
  double ret = chanThroughput[i][j];
  if (child[i][j][0] != -1 && child[i][j][1] != -1) {
    double a1 = solve(i, child[i][j][0]);
    double a2 = solve(i, child[i][j][1]);

    if (a1 + a2 > ret) {
      ret = a1 + a2;
      inSolution[i][j] = false;
    }
  }

  return ret;
}

double calcDP(Solution &sol) { //TODO
  double OF = 0.0;
  for (int s = 0; s < sol.spectrums.size(); s++) {
    for (int c = 0; c < sol.spectrums[s].channels.size(); c++) {
      if (parent[s][c] == -1) {
        OF += solve(s, c);
      }
    }
  }

  return OF;
}

bool reinsert(Solution &sol, Connection conn, ii from, ii to, bool force = false) {
  if (to == from)
    return false;

  Solution copy(sol);
  bool improved = false;

  Channel oldChan = deleteFromChannel(copy.spectrums[from.first].channels[from.second], conn.id);
  Channel newChan = insertInChannel(copy.spectrums[to.first].channels[to.second], conn.id);

  swap(copy.spectrums[from.first].channels[from.second], oldChan);
  swap(copy.spectrums[to.first].channels[to.second], newChan);

  if (copy > sol)
    improved = true;

  if ((copy > sol) || force) {
    sol = copy;
  }

  return improved;
}

void K_AddDrop(Solution &sol, int K) {
  //K = min(K, conexoes_nao_alocadas); //TODO
  ii channelFrom = {3, 0};
  for (int i = 0; i < K; i++) {
    int idx = rng.randInt(sol.spectrums[3].channels[0].connections.size() - 1);
    Connection conn = sol.spectrums[3].channels[0].connections[idx];

    int a = rng.randInt(sol.spectrums.size() - 1);
    int b = rng.randInt(sol.spectrums[a].channels.size() - 1);
    ii chanelTo = {a, b};

    reinsert(sol, conn, channelFrom, chanelTo, true);
  }
}

void K_RemoveAndInserts(Solution &sol, int K) {
  int k = 0;
  while (k < K) {
    int a = rng.randInt(sol.spectrums.size() - 1);
    int b = rng.randInt(sol.spectrums[a].channels.size() - 1);

    if (sol.spectrums[a].channels[b].connections.empty())
      continue;

    k++;

    int z = rng.randInt(sol.spectrums[a].channels[b].connections.size() - 1);
    Connection conn = sol.spectrums[a].channels[b].connections[z];

    ii from = {a, b};
    reinsert(sol, conn, from, {3, 0}, true);
  }
  K_AddDrop(sol, K);
}

Solution multipleRepresentation(Solution ret) {
  memset(parent, -1, sizeof parent);
  memset(child, -1, sizeof child);

  for (int s = 0; s < ret.spectrums.size(); s++) {
    Spectrum &spectrum = ret.spectrums[s];
    int curr = 0;
    while (curr + 1 < spectrum.channels.size()) {
      if (spectrum.channels[curr].bandwidth == spectrum.channels[curr + 1].bandwidth) {
        Channel merged = spectrum.channels[curr];
        merged.bandwidth *= 2;

        for (const Connection &conn : spectrum.channels[curr].connections) {
          merged.connections.emplace_back(conn);
        }

        int p = spectrum.channels.size();
        child[s][p][0] = curr;
        child[s][p][1] = curr + 1;
        parent[s][curr] = p;
        parent[s][curr + 1] = p;

        spectrum.channels.emplace_back(merged);
        curr += 2;
      } else {
        curr += 1;
      }
    }
  }

  computeThroughput(ret);
  return ret;
}

void recoverSolution(int i, int j, bool clean) {
  if (clean)
    inSolution[i][j] = false;

  if (inSolution[i][j])
    clean = true;

  if (child[i][j][0] != -1 && child[i][j][1] != -1) {
    recoverSolution(i, child[i][j][0], clean);
    recoverSolution(j, child[i][j][1], clean);
  }
}

Solution explicitSolution(const Solution &curr) {
  Solution ret;

  for (int s = 0; s < curr.spectrums.size(); s++) {
    for (int c = 0; c < curr.spectrums[s].channels.size(); c++) {
      if (parent[s][c] == -1) {
        recoverSolution(s, c, false);
      }
    }
  }

  for (int s = 0; s < MAX_SPECTRUM; s++) {
    for (int c = 0; c < MAX_CHANNELS; c++) {
      if (inSolution[s][c]) {
        ret.spectrums[s].channels.emplace_back(curr.spectrums[s].channels[c]);
      }
    }
  }

  return ret;
}

Solution newVNS_Reinsert(Solution &multiple, Solution &curr) {
  bool improved = false;
  do {
    improved = false;

    for (int i = 0; i < nConnections; i++) {
      Solution multipleClean(multiple);

      for (Spectrum &spectrum : multipleClean.spectrums) {
        for (Channel &channel : spectrum.channels) {
          for (const Connection &conn : channel.connections) {
            if (conn.id == i) {
              channel = deleteFromChannel(channel, conn.id);
              break;
            }
          }
        }
      }

      Solution multipleContaining(multipleClean);
      for (Spectrum &spectrum : multipleClean.spectrums) {
        for (Channel &channel : spectrum.channels) {
          channel = insertInChannel(channel, i);
        }
      }

      double bestOF = -1;
      ii bestChannel = {-1, -1};

      for (int s = 0; s < curr.spectrums.size(); s++) {
        for (int c = 0; c < curr.spectrums[s].channels.size(); c++) {
          setDP(multipleClean);
          int currChan = c;

          while (currChan != -1) {
            chanThroughput[s][currChan] = multipleContaining.spectrums[s].channels[currChan].throughput;
            currChan = parent[s][currChan];
          }

          double OF = calcDP(multiple);
          if (OF > bestOF) {
            bestOF = OF;
            bestChannel = {s, c};
          }
        }
      }

      if (bestOF > multiple.totalThroughput) {
        improved = true;

        for (int s = 0; s < curr.spectrums.size(); s++) {
          for (int c = 0; c < curr.spectrums[s].channels.size(); c++) {
            Channel &channel = curr.spectrums[s].channels[c];
            for (const Connection &conn : channel.connections) {
              if (conn.id == i) {
                channel = deleteFromChannel(channel, conn.id);
                break;
              }
            }
          }
        }

        int newSpec = bestChannel.first;
        int newChan = bestChannel.second;

        while (newChan != -1) {
          multiple.spectrums[newSpec].channels[newChan] = insertInChannel(multiple.spectrums[newSpec].channels[newChan],
                                                                          i);
          newChan = parent[newSpec][newChan];
        }
        computeThroughput(multiple);
      }

    }
  } while (improved);
  return explicitSolution(multiple);
}

Solution VNS(Solution initsol) {
  Solution delta = convert20MHz(initsol);
  Solution rep = multipleRepresentation(delta);

  setDP(rep);
  calcDP(rep);

  Solution explicitSol = explicitSolution(rep);
  Solution star = explicitSol;
  Solution localMax = delta;

  int K_MUL = max(1, nConnections / 100);
  int K_MAX = 50;
  while (!stop()) {
    int k = 1;
    while (k <= K_MAX && !stop()) {
      delta = localMax;
      if (rng.randInt(1)) { //AddDrop
        K_AddDrop(delta, k * K_MUL);
      } else { //Reinsert
        K_RemoveAndInserts(delta, k * K_MUL);
      }
      Solution multiple = multipleRepresentation(delta);
      setDP(multiple);
      calcDP(multiple);

      Solution aux = newVNS_Reinsert(multiple, delta);
      delta = convert20MHz(aux);

      if (aux > localMax) {
        k = 1;
        localMax = aux;
      } else {
        k++;
      }

      if (localMax > star) {
        star = aux;
      }
    }
  }

  return star;  //todo
}

void init(const string &openingFile = "", double timeLimit = 10) {
#ifdef DEBUG_CLION //TODO: remind to remove the MACRO before real tests
  puts("============== WITH DEBUG ==============");
  freopen("/Users/jjaneto/Downloads/codes_new/BRKGA_FF_Best/Instancias/D250x250/U_8/U_8_1.txt", "r", stdin);
#else
  if (!openingFile.empty()) {
    fprintf(stderr, "trying to open input file %s\n", openingFile.c_str());
    freopen(openingFile.c_str(), "r", stdin);
  }
#endif

  if (stdin == NULL) {
    fprintf(stderr, "error opening input file (stdin)\n");
    exit(13);
  }

  loadData();
  maximumTime = timeLimit;

  fprintf(stdout, "will execute for %lf seconds\n", maximumTime);
}

int main(int argc, char *argv[]) {
  if (argc != 4) {
    fprintf(stderr, "wrong arguments. Must be: stdin, solutionFile, objectiveFile\n");
    exit(13);
  }

  FILE *solutionFile = fopen(argv[2], "a");
  if (solutionFile == NULL) {
    fprintf(stderr, "error opening solutionFile file\n");
    exit(13);
  }

  FILE *objectivesFile = fopen(argv[3], "a");
  if (objectivesFile == NULL) {
    fprintf(stderr, "error opening objectivesFile file\n");
    exit(13);
  }

  init(argv[1], 10);

  Solution aux = createSolution();
  fprintf(objectivesFile, "%lf\n", aux.totalThroughput);
  printf("%lf\n", aux.totalThroughput);
  aux.printSolution(solutionFile);

  fclose(solutionFile);
  fclose(objectivesFile);
  return 0;
}