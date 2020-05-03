//
// Created by José Joaquim on 03/03/20.
//

#ifndef BRKGA_FF_BEST_HEURISTICDECODER_H
#define BRKGA_FF_BEST_HEURISTICDECODER_H

#include <list>
#include <vector>
#include <algorithm>
#include <fstream>
#include <cmath>
#include <iostream>
#include <cstdio>
#include <deque>
#include <vector>
#include <stack>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <map>
#include <cstring>
#include <cassert>
#include <numeric>

#include <thread>
#include <chrono>

#include "MTRand.h"

extern int channels20MHz[25];
//extern int channels40MHz[12] = {25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36};
//extern int channels80MHz[6] = {37, 38, 39, 40, 41, 42};
//extern int channels160MHz[2] = {43, 44};

extern int overlap[45][45];
extern std::string PATH_TO[46][46];
extern std::unordered_map<int, std::pair<int, int>> mapChtoCh;
extern MTRand rng;

extern int nConnections;
extern double dataRates[10][4];
extern double distanceMatrix[2048][2048], interferenceMatrix[2048][2048];
extern double senders[2048][2], receivers[2048][2];
extern std::vector<std::vector<double>> SINR;
extern double powerSender, alfa, noise, ttm;
extern std::map<int, std::vector<int>> chToLinks;

struct Channel {
    double throughput;
    double interference;
    short parent;
    short childs[2];
    short id;
    short bw;
    std::vector<short> links;
};

struct Link {
    int _idR, _idS;
    int id;
    int ch;
    int bw;
    int origChannel;
    double distanceSenderReceiver;
    double interference;
    double SINR;
    int MCS;

    Link();

    Link(int id);

    Link(const Link &x);

    int getId() const;

    int getChannel() const;

    void operator=(const Link &x);

    friend bool operator==(const Link &oq, const Link &o2);

    friend bool operator<(const Link &o1, const Link &o2);

    friend bool operator>(const Link &o1, const Link &o2);

    void setChannel(int ch);

    void printLink() const;
};

class Solution {
public:
    double objective;
    std::deque<Link> scheduled_links;

    std::vector<int> zeroLinks;

    bool objectiveFlag;

    Solution();

    Solution(const Solution &o1);

    Solution(double objective, const std::deque<Link> &scheduledLinks);

    Solution(double objective);

    void setZeroLinks(const std::vector<int> &zeroLinks);

    std::vector<int> getZeroLinks();

    void computeObjective(bool show = false);

    void computeInterference();

    void clearChannel(int channel);

    void insert(const Link &l);

    void setObjectiveFlag(bool value);

    bool removeLink(Link link);

    int getNumberOfScheduledLinks() const;

    Link removeLinkByIndex(int index);

    void exchangeLinks(int idOldLink, Link newLink);

    std::vector<int> getScheduledChannels() const;

    double getChannelThroughput(int channel) const;

    void setChannelOfLink(int id, int channel);

    void setScheduledLinks(const std::deque<Link> &newLinks);

    void addLinks(const std::deque<Link> &links);

    double getObjective(bool force = false) const;

    std::deque<Link> getScheduledLinks() const;

    std::deque<Link> getLinksInChannel(int ch) const;

    friend bool operator<(const Solution &o1, const Solution &o2);

    friend bool operator>(const Solution &o1, const Solution &o2);

    friend bool operator>=(const Solution &o1, const Solution &o2);

    friend bool operator<=(const Solution &o1, const Solution &o2);

    friend bool operator==(const Solution &o1, const Solution &o2);

    void operator=(const Solution &o) {
      objective = o.objective;
      scheduled_links = o.scheduled_links;
      objectiveFlag = o.objectiveFlag;
    }
};

class HeuristicDecoder {
private:
    clock_t TempoFO_StarInic;
    std::fstream dataFileRead;
    std::string fileName;

public:

    HeuristicDecoder();

    double decode(std::vector<double> &chromosome) const;

    int getQuantConnections();

    void setInitialTime();

    Solution generateSolution();
};

int whichBw(int ch);

int bwIdx(int bw);

void split(Solution &dest, Solution &src, int ch, bool ok = false);

void dfs(int u, int pai = -1, std::string path = "");

#endif //BRKGA_FF_BEST_HEURISTICDECODER_H
