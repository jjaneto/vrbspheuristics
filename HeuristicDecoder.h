//
// Created by José Joaquim on 03/03/20.
//

#ifndef BRKGA_FF_BEST_HEURISTICDECODER_H
#define BRKGA_FF_BEST_HEURISTICDECODER_H

#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <cstdio>
#include <vector>
#include <cstring>
#include <cassert>
#include <numeric>
#include <random>
#include <algorithm>

#include <thread>
#include <chrono>

#include "MTRand.h"

const int MAX_CONN = 2048;

extern int nConnections, nSpectrums;
extern double dataRates[10][4];
extern double distanceMatrix[MAX_CONN][MAX_CONN], interferenceMatrix[MAX_CONN][MAX_CONN];
extern double senders[MAX_CONN][2], receivers[MAX_CONN][2];
extern std::vector<std::vector<double>> SINR;
extern std::pair<int, int> zeroChannel;
extern double powerSender, alfa, noise, ttm;
extern MTRand rng;
extern const int MAX_SPECTRUM, MAX_CHANNELS;

using ii = std::pair<int, int>;

struct Connection {
    int id;
    double throughput;
    double interference;
    double SINR;
    double distanceSR;

    Connection(int id, double throughput, double interference, double distanceSR);

    Connection(int id);

    bool operator<(const Connection &other) const {
        return distanceSR < other.distanceSR;
    }

    bool operator>(const Connection &other) const {
        return !operator<(other);
    }
};

struct Channel {
    double throughput; //TODO: So far, I do not have any use for this variable.
    double interference;
    int bandwidth;
    std::vector<Connection> connections;

    bool operator<(const Channel &other) const {
        return bandwidth < other.bandwidth;
    }

    bool operator>(const Channel &other) const {
        return !operator<(other);
    }

    Channel(double throughput, double interference, int bandwidth, const std::vector<Connection> &connections);

    Channel(int bandwidth);
};

struct Spectrum {
    int maxFrequency;
    int usedFrequency;
    std::vector<Channel> channels;

    Spectrum(int maxFrequency, int usedFrequency, const std::vector<Channel> &channels);
};

class Solution {
public:
    std::vector<Spectrum> spectrums;
    double totalThroughput;

    bool throughputFlag; //FIXME: conflict with the constructors

    Solution(const std::vector<Spectrum> &spectrums, double totalThroughput, bool throughputFlag);

    Solution();

    void printSolution(FILE *file = nullptr);

    //AVAILABLE ONLY IN THE META-HEURISTICS FILES
    double decode(std::vector<double> variables) const;

    friend bool operator>(const Solution &o1, const Solution &o2);

    friend bool operator<(const Solution &o1, const Solution &o2);
};

bool checkOne(const Solution &s);

bool checkTwo(const Solution &s);

int bwIdx(int bw);

void loadData();

void rawInsert(Solution &sol, int conn, ii where);

void rawRemove(Solution &sol, int conn, ii where);

Channel insertInChannel(Channel newChannel, int conn);

Channel deleteFromChannel(const Channel &channel, int conn);

double computeThroughput(Solution &curr, bool force = false);

bool double_equals(double a, double b, double epsilon = 0.000000001);

void computeChannelsThroughput(std::vector<Channel> &channels);

double computeConnectionThroughput(Connection &conn, int bandWidth, bool force = false);

Solution createSolution();

#endif //BRKGA_FF_BEST_HEURISTICDECODER_H
