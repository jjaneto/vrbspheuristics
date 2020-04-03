#include "HeuristicDecoder.h"

const int K_MAX = 50;

HeuristicDecoder *heu;
vector<Link> nonScheduledLinks;
double maximumTime;
clock_t startTime;

bool isStoppingCriteriaReached() {
  return (((double) (clock() - startTime)) / CLOCKS_PER_SEC) >= maximumTime;
}

int betaAddDrop(const Solution current, Link link) {
  vector<int> channels = current.getScheduledChannels();
  int bestChannel = -1;

  Solution dummy;

  for (int ch : channels) {
    Solution aux(current);
    Link copyLink(link);
    copyLink.setChannel(ch);

    aux.insert(copyLink);

    if (aux > dummy) { //FIXME: Should I return the best solution from this loop or only if its better than current?
      dummy = aux;
      bestChannel = ch;
    }
  }

  return bestChannel;
}

void betaAddDrop(Solution &current, int beta = 1) {
  for (int i = 0; i < beta; i++) {
    assert(nonScheduledLinks.empty());
    int rndIndex = rng.randInt(nonScheduledLinks.size() - 1);
//    printf("betaAddDrop no link %d\n", nonScheduledLinks[rndIndex].id);
    int bestChannel = betaAddDrop(current, nonScheduledLinks[rndIndex]);
    assert(bestChannel != -1);

    Link toInsert(nonScheduledLinks[rndIndex]);
    toInsert.setChannel(bestChannel);
    current.insert(toInsert);

    swap(nonScheduledLinks[rndIndex], nonScheduledLinks.back());
    nonScheduledLinks.pop_back();
  }
}

void betaReinsert(Solution &current, int beta = 1) {
  unordered_set<int> usedLinks;
  for (int i = 0; i < beta; i++) {
    int rndIndex = rng.randInt(current.getNumberOfScheduledLinks() - 1);

    while (usedLinks.count(rndIndex))
      rndIndex = rng.randInt(current.getNumberOfScheduledLinks() - 1);

    Link rmvLink = current.removeLinkByIndex(rndIndex);
    int bestChannel = betaAddDrop(current, rmvLink);
    assert(bestChannel != -1);

    rmvLink.setChannel(bestChannel);
    current.insert(rmvLink);
  }
}

Solution solve(Solution S, int channel) { //TODO: need to implement origChannel variable.
  if (bwIdx(channel) > 0) {
    Solution S1, S2;

    for (const Link &link_ : S.getScheduledLinks()) {
      char moveTo = PATH_TO[link_.ch][link_.origChannel][bwIdx(link_.ch)];
      assert(moveTo == 'L' || moveTo == 'R');

      Link copy(link_);
      if (moveTo == 'L') {
        copy.setChannel(mapChtoCh[copy.ch].first);
        S1.insert(copy);
      } else if (moveTo == 'R') {
        copy.setChannel(mapChtoCh[copy.ch].second);
        S2.insert(copy);
      }
    }

    Solution solToLeft = solve(S1, mapChtoCh[channel].first);
    Solution solToRight = solve(S2, mapChtoCh[channel].second);

    if (solToLeft.getObjective() + solToRight.getObjective()) {
      deque<Link> newLinks(solToLeft.getScheduledLinks());
      for (const Link &link_ : solToRight.getScheduledLinks()) {
        newLinks.emplace_back(link_);
      }

      S.setScheduledLinks(newLinks);
      S.computeObjective(false);
    }
  }

  return S;
}

Solution buildRootSolution(const Solution &curr, int root) {
    Solution ret;

    for (const Link &x : curr.getScheduledLinks()) {
      if (overlap[x.getChannel()][root]) {
        Link newLink(x);
        newLink.setChannel(root);
        ret.insert(newLink);
      }
    }

    return ret;
}

Solution solve(Solution S) {
  Solution arr[5];

  arr[0] = solve(buildRootSolution(S, 45), 45);
  arr[1] = solve(buildRootSolution(S, 44), 44);
  arr[2] = solve(buildRootSolution(S, 43), 43);
  arr[3] = solve(buildRootSolution(S, 42), 42);
  arr[4] = solve(buildRootSolution(S, 25), 25);

  Solution final_;
  for (const Solution &s_ : arr) {
    final_.addLinks(s_.getScheduledLinks());
  }

  return final_;
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

Solution localSearch(Solution current) {
  bool improve = true;
  Solution S_star = convert(current);
  Solution S(S_star), S_2(S_star);

  while (improve) {
    for (Link active : S.getScheduledLinks()) {
      //TODO: Do I need SCopy here? Kinda lost
      //TODO: see split() function. Variable origChannel with wrong values.

      for (int channel20 : channels20MHz) {
        if (channel20 == active.getChannel())
          continue;

        Solution S_1(S);
//      active.printLink();
        assert(S_1.removeLink(active)); //Testando se a funcao vai funcionar corretamente

        Link aux(active);
        aux.setChannel(channel20);
        S_1.insert(aux);
        solve(S_1);

        S_2 = (S_1 > S_2) ? S_1 : S_2;
      }

      if (S_2 > S_star) {
        S_star = S_2;
        S = S_2;
      } else {
        improve = false;
      }
    }
  }

  return S;
}

Solution pertubation(Solution S, int k, const int NUMBER_OF_LINKS) {
  double threeshold = rng.rand();
  if (threeshold >= .5) {
    betaAddDrop(S, lround(NUMBER_OF_LINKS * k * 0.01));
  } else {
    betaReinsert(S, lround(NUMBER_OF_LINKS * k * 0.01));
  }
  return S;
}

void init() {
#ifdef DEBUG_CLION
  puts("WITH DEBUG");
  freopen("/Users/jjaneto/Downloads/codes_new/BRKGA_FF_Best/Instancias/D250x250/U_8/U_8_1.txt", "r", stdin);
#endif

  heu = new HeuristicDecoder();
  for (int i = 0; i < heu->getQuantConnections(); i++) {
    nonScheduledLinks.emplace_back(Link(i));
  }
  maximumTime = 10;
}

int main(int argc, char *argv[]) {
  init();
  const int NUMBER_OF_LINKS = heu->getQuantConnections();
  startTime = clock();
  //--------------------------------------------------------
  Solution S_dummy = heu->generateSolution();
  Solution S = localSearch(S_dummy);
  while (!isStoppingCriteriaReached()) {
    int k = 1;
    while (k < K_MAX) {
      Solution S_1 = pertubation(S, k, NUMBER_OF_LINKS);
      Solution S_2 = localSearch(S_1);

      if (S_2.getObjective() > S_1.getObjective()) {
        S = S_2;
        k = 1;
      } else {
        k++;
      }
    }
  }

  delete heu;
  return 0;
}