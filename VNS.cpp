#include <cstdio>
#include <vector>
#include <map>
#include <set>
#include <deque>
#include <algorithm>
#include <cstring>
#include <cmath>
#include "HeuristicDecoder.h"
#include "MTRand.h"

const int K_MAX = 50;

MTRand rng;
HeuristicDecoder heu;
vector<Link> nonScheduledLinks;
double maximumTime;
clock_t startTime;

bool isStoppingCriteriaReached() {
  return (((double) (clock() - startTime)) / CLOCKS_PER_SEC) >= maximumTime;
}

Solution betaAddDrop(const Solution current, Link link) {
  vector<int> channels = current.getScheduledChannels();

  Solution newly;
  for (int ch : channels) {
    Solution aux(current);
    Link copyLink(link);
    copyLink.setChannel(ch);

    aux.insert(copyLink);

    if (aux > newly) {
      newly = aux;
    }
  }

  //TODO: and...? What goes next?
  return newly;
}

void betaAddDrop(const Solution current, int beta = 1) {
  for (int i = 0; i < beta; i++) {
    int rndIndex = rng.randInt(nonScheduledLinks.size() - 1);
    betaAddDrop(current, nonScheduledLinks[rndIndex]); //TODO: devo retornar algo?

    swap(nonScheduledLinks[rndIndex], nonScheduledLinks.back());
    nonScheduledLinks.pop_back();
  }
}

void betaReinsert(const Solution current, int beta = 1) {
  assert(beta >= current.getNumberOfScheduledLinks());

  Solution currentCopy(current), currentDummy(current);
  for (int i = 0; i < beta; i++) {
    int rndIndex = rng.randInt(currentDummy.getNumberOfScheduledLinks() - 1); //TODO: verificar se o rand eh inclusive

    Link rmvLink = currentDummy.removeLinkByIndex(rndIndex);
//    betaAddDrop(currentDummy, rmvLink); //FIXME: I'm passing rmvLink by value; therefore, nothing changes after here.
//
//    Link newLink(rmvLink); //TODO: tem que setar o novo canal.
//    currentCopy.exchangeLinks(rmvLink.getId(), newLink);

    currentCopy = betaAddDrop(currentDummy, rmvLink);
  }
}

Solution solve(Solution S, int channel) {
  //TODO: best <- canal.throughput (?????)
  double best = S.getChannelThroughput(channel);

  if (whichBw(channel) > 20) {
    //TODO
//    double bestC1 = solve(S, mapChtoCh[channel].first); //FIXME: it's not double, but Solution
//    double bestC2 = solve(S, mapChtoCh[channel].second);

//    if (bestC1 + bestC2 > best) {
//      best = bestC1 + bestC2;
//      TODO: canal.inSolution = false
//    }
  }

  return Solution();
}

bool allChannels20MHz(const Solution &check) {
  for (const int &x : check.getScheduledChannels()) {
    if (whichBw(x) > 20)
      return false;
  }

  return true;
}

Solution convert(const Solution &aux) { //TODO: colocar todo mundo em canal de 20MHz
  Solution ret;

  while (allChannels20MHz(aux)) {
    for (Link &x : aux.getScheduledLinks()) {
      if (whichBw(x.getChannel()) > 20) {
        split(ret, ret, x.getChannel());
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
    for (Link active : current.getScheduledLinks()) {
      Solution SCopy = S;
      Solution S_1 = SCopy;
      S_1.removeLink(active);
      for (int channel20 : channels20MHz) {
        Link aux(active);
        aux.setChannel(channel20);
        S_1.insert(aux);
        solve(S_1, channel20);

        S_2 = (S_1 > S_2) ? S_1 : S_2;

        //S_1.removeLink(active)
      }

      if (S_2 > S_star) {
        S_star = S_2;
        S = S_2;
      } else {
        improve = false;
      }
    }
  }

  return S; //FIXME: o que eh para retornar?
}

Solution pertubation(Solution S, int k, const int NUMBER_OF_LINKS) {
  betaAddDrop(S, NUMBER_OF_LINKS % k);
  betaReinsert(S, NUMBER_OF_LINKS % k);
  Solution ret;
  return ret; //TODO: O que eh para retornar?
}

void init() {
  for (int i = 0; i < heu.getQuantConnections(); i++) {
    nonScheduledLinks.emplace_back(Link(i));
  }

  maximumTime = 10;

}

int main(int argc, char *argv[]) {
  init();
  const int NUMBER_OF_LINKS = heu.getQuantConnections();
  startTime = clock();
  //--------------------------------------------------------
  Solution S_dummy = heu.generateSolution();
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
        k++; //TODO: isso mesmo?
      }
    }
  }
  return 0;
}