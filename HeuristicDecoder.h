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
#include <set>
#include <map>

#include <thread>
#include <chrono>

using namespace std;

struct Link {
    int _idR, _idS;
    int id;
    int ch;
    int bw;
    double distanceSenderReceiver;
    double interference;
    double SINR;
    int MCS;

    Link();

    Link(int id);

    Link(const Link &x);

    void operator=(const Link &x);

    friend bool operator==(const Link &oq, const Link &o2);

    void setChannel(int ch);

    void printLink() const;
};

class Solution {
public:
    double objective;
    deque<Link> scheduled_links;

    Solution();

    Solution(const Solution &o1);

    Solution(double objective, const deque<Link> &scheduledLinks);

    Solution(double objective);

    void computeObjective(bool show = false);

    void computeInterference();

    void clearChannel(int channel);

    void insert(const Link &l);

    double getObjective() const;

    deque<Link> getScheduledLinks() const;

    deque<Link> getLinksInChannel(int ch) const;

    friend bool operator<(const Solution &o1, const Solution &o2);

    friend bool operator>(const Solution &o1, const Solution &o2);

    friend bool operator==(const Solution &o1, const Solution &o2);

    void operator=(const Solution &o) {
      objective = o.objective;
      scheduled_links = o.scheduled_links;
    }
};

class HeuristicDecoder {
private:
    clock_t TempoFO_StarInic;
    fstream dataFileRead;
    string fileName;

public:

    HeuristicDecoder();

    double decode(std::vector<double> &chromosome) const;

    int getQuantConnections();

    void setInitialTime();
};

#endif //BRKGA_FF_BEST_HEURISTICDECODER_H
