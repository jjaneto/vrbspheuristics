#include "HeuristicDecoder.h"
#include <set>

using namespace std;

const int MAX_SPECTRUM = 4;
const int MAX_CHANNELS = 45;
double maximumTime;
clock_t startTime;
int parent[MAX_SPECTRUM][MAX_CHANNELS];
int child[MAX_SPECTRUM][MAX_CHANNELS][2];
double chanThroughput[MAX_SPECTRUM][MAX_CHANNELS];
bool inSolution[MAX_SPECTRUM][MAX_CHANNELS];

bool stop() {
    return (((double) (clock() - startTime)) / CLOCKS_PER_SEC) >= maximumTime;
}

bool allChannels20MHz(const Solution &sol) {
    for (const Spectrum &spectrum : sol.spectrums) {
        for (const Channel &channel : spectrum.channels) {
            if (channel.bandwidth > 20) {
                return false;
            }
        }
    }
    return true;
}

Solution convert20MHz(Solution sol) {
    for (Spectrum &spectrum : sol.spectrums) {
        sort(spectrum.channels.rbegin(), spectrum.channels.rend());
    }

    while (!allChannels20MHz(sol)) {
        for (Spectrum &spectrum : sol.spectrums) {
            for (int c = 0; c < spectrum.channels.size(); c++) {
                if (spectrum.channels[c].bandwidth <= 20)
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

void setDP(Solution &sol, bool ok = false) {
    memset(chanThroughput, 0, sizeof chanThroughput);
    memset(inSolution, false, sizeof inSolution);
//    for (int i = 0; i < MAX_SPECTRUM; i++) {
//        for (int j = 0; j < MAX_CHANNELS; j++) {
//            chanThroughput[i][j] = 0.0;
//            inSolution[i][j] = false;
//        }
//    }

    for (int s = 0; s < sol.spectrums.size(); s++) {
//        computeChannelsThroughput(sol.spectrums[s].channels);
        for (int c = 0; c < sol.spectrums[s].channels.size(); c++) {
            chanThroughput[s][c] = sol.spectrums[s].channels[c].throughput;
            if (ok) {
                printf("throughput of {%d, %d} eh %.3lf\n", s, c, chanThroughput[s][c]);
            }
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

double calcDP(Solution &sol, bool ok = false) {
    double OF = 0.0;
    for (int s = 0; s < sol.spectrums.size(); s++) {
        for (int c = 0; c < sol.spectrums[s].channels.size(); c++) {
            if (parent[s][c] == -1) {
                double ret = solve(s, c);
                if (ok) {
                    printf("will add %.3lf, starting from {%d, %d}\n", ret, s, c);
                }
                OF += ret;
            }
        }
    }

    return OF;
}

bool reinsert(Solution &sol, double &objective, Connection conn, ii from, ii to, bool force = false) {
    if (to == from)
        return false;

//  Solution copy(sol);

    Channel &oldChan = sol.spectrums[from.first].channels[from.second];
    Channel &newChan = sol.spectrums[to.first].channels[to.second];

    Channel copyOldChan = deleteFromChannel(oldChan, conn.id);
    Channel copyNewChan = insertInChannel(newChan, conn.id);

//  Channel oldChan = deleteFromChannel(copy.spectrums[from.first].channels[from.second], conn.id);
//  Channel newChan = insertInChannel(copy.spectrums[to.first].channels[to.second], conn.id);
//
//  swap(copy.spectrums[from.first].channels[from.second], oldChan);
//  swap(copy.spectrums[to.first].channels[to.second], newChan);

    double newObjective =
            objective - oldChan.throughput - newChan.throughput + copyNewChan.throughput + copyOldChan.throughput;

    bool improved = false;
    if (newObjective > objective)
        improved = true;

    if ((newObjective > objective) || force) {
        oldChan = copyOldChan;
        newChan = copyNewChan;
//    computeThroughput(sol); //TODO: do I have to do this?
        //
        objective = newObjective;
    }

    return improved;
}

void K_AddDrop(Solution &sol, double &_FO_delta, int K) {
    ii chFrom = {sol.spectrums.size() - 1, 0};
    K = min(K, int(sol.spectrums[chFrom.first].channels[chFrom.second].connections.size()));
    for (int i = 0; i < K; i++) {
        int idx = rng.randInt(sol.spectrums[chFrom.first].channels[chFrom.second].connections.size() - 1);
        Connection conn = sol.spectrums[chFrom.first].channels[chFrom.second].connections[idx];

        int a = rng.randInt(sol.spectrums.size() - 1);
        int b = rng.randInt(sol.spectrums[a].channels.size() - 1);
        ii channelTo = {a, b};

        reinsert(sol, _FO_delta, conn, chFrom, channelTo, true);
    }
}

void K_RemoveAndInserts(Solution &sol, double &_FO_delta, int K) {
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
        reinsert(sol, _FO_delta, conn, from, {sol.spectrums.size() - 1, 0}, true);
    }
    K_AddDrop(sol, _FO_delta, K);
}

Solution multipleRepresentation(Solution ret) {
    memset(parent, -1, sizeof parent);
    memset(child, -1, sizeof child);

    for (int s = 0; s < ret.spectrums.size(); s++) {
        Spectrum &spectrum = ret.spectrums[s];
        int curr = 0;
        while (curr < spectrum.channels.size()) {
            if (curr + 1 < spectrum.channels.size() &&
                spectrum.channels[curr].bandwidth == spectrum.channels[curr + 1].bandwidth) {
                Channel merged = spectrum.channels[curr];
                merged.bandwidth *= 2;

                for (const Connection &conn : spectrum.channels[curr + 1].connections) {
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

    if (inSolution[i][j]) {
        clean = true;
    }

    if (child[i][j][0] != -1 && child[i][j][1] != -1) {
        recoverSolution(i, child[i][j][0], clean);
        recoverSolution(i, child[i][j][1], clean);
    }
}

Solution explicitSolution(const Solution &curr) {
    Solution ret;
    ret.spectrums.resize(curr.spectrums.size(), Spectrum(0.0, 0.0, vector<Channel>()));

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
    computeThroughput(ret);
    return ret;
}

bool fixChannels(Solution &sol, double &_FO_delta) {
    bool improved = false;
    ii zeroChannel = {sol.spectrums.size() - 1, 0};
    do {
        improved = false;
        for (int s = 0; s < sol.spectrums.size(); s++) {
            for (int c = 0; c < sol.spectrums[s].channels.size(); c++) {
                if (make_pair(s, c) == zeroChannel)
                    continue;

                int idx = 0;
                while (idx < sol.spectrums[s].channels[c].connections.size()) {
                    Connection &conn = sol.spectrums[s].channels[c].connections[idx];
                    if (double_equals(computeConnectionThroughput(conn, sol.spectrums[s].channels[c].bandwidth), 0.0)) {
                        reinsert(sol, _FO_delta, conn, {s, c}, zeroChannel, true);
                        improved = true;
                    } else {
                        idx++;
                    }
                }
            }
        }
    } while (improved);
    computeThroughput(sol);
    return improved;
}

Solution newVNS_Reinsert(Solution &multiple, double &_FO_delta, Solution &curr) {
    bool improved = false;
    int cnt = 0;
    double initLoopTime = clock(), loopCleanContaining = 0.0, loopBestCh = 0.0, updtSol = 0.0, loopEachConn = 0.0, aux = 0.0;
    do {
        double aux0 = clock();
        improved = false;
        for (int i = 0; i < nConnections; i++) {
            double loopEachConn0 = clock();
            Solution multipleClean(multiple);
            double loopCleanContaining0 = clock();
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
            for (Spectrum &spectrum : multipleContaining.spectrums) {
                for (Channel &channel : spectrum.channels) {
                    channel = insertInChannel(channel, i);
                }
            }
            loopCleanContaining = max(loopCleanContaining,
                                      (((double) (clock() - loopCleanContaining0) / CLOCKS_PER_SEC)));
            double bestOF = -1;
            ii bestChannel = {-1, -1};

            double loopBestCh0 = clock();
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
            loopBestCh = max(loopBestCh, (((double) (clock() - loopBestCh0) / CLOCKS_PER_SEC)));

            bool enter = false;
            if (bestOF > _FO_delta) {
                enter = true;
//                printf("%lf %lf diff %lf\n", bestOF, _FO_delta, bestOF - _FO_delta);
                double updtSol0 = clock();
                improved = true;
                _FO_delta = bestOF;

                for (int s = 0; s < multiple.spectrums.size(); s++) {
                    for (int c = 0; c < multiple.spectrums[s].channels.size(); c++) {
                        Channel &channel = multiple.spectrums[s].channels[c];
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
                    multiple.spectrums[newSpec].channels[newChan] = insertInChannel(
                            multiple.spectrums[newSpec].channels[newChan],
                            i);
                    newChan = parent[newSpec][newChan];
                }
                updtSol = max(updtSol, (((double) (clock() - updtSol0) / CLOCKS_PER_SEC)));
            }
//            if ((((double) (clock() - loopEachConn0) / CLOCKS_PER_SEC)) > loopEachConn) {
//                printf("new high, %lf %lf did upgrade %d\n", (((double) (clock() - loopEachConn0) / CLOCKS_PER_SEC)),
//                       loopEachConn, enter);
//            }
            loopEachConn = max(loopEachConn, (((double) (clock() - loopEachConn0) / CLOCKS_PER_SEC)));
//            printf("this loop takes %lf seconds\n", (((double) (clock() - loopEachConn0) / CLOCKS_PER_SEC)));
        }
        printf("demorou %lf\n", (((double) (clock() - aux0) / CLOCKS_PER_SEC)));
        aux = max(aux, (((double) (clock() - aux0) / CLOCKS_PER_SEC)));
    } while (improved);

//    printf("(1) loop clean contaning lasts %lf\n", loopCleanContaining);
//    printf("(2) search for the best channel lasts %lf\n", loopBestCh);
//    printf("(3) loop for conn lasts %lf\n", loopEachConn);
//    printf("(4) loop for each connection lasts %lf\n", aux);
//    printf("(X) update solution lasts %lf\n", updtSol);
//    printf("===== total do-while lasts %lf =====\n", (((double) (clock() - initLoopTime)) / CLOCKS_PER_SEC));
//    this_thread::sleep_for(chrono::seconds(5));
    return explicitSolution(multiple);
}

bool checkOne(const Solution &s) {
    set<int> x;

    for (int i = 0; i < s.spectrums.size(); i++) {
        for (int j = 0; j < s.spectrums[i].channels.size(); j++) {
            for (const Connection &connection : s.spectrums[i].channels[j].connections) {
                if (x.find(connection.id) != x.end()) {
                    fprintf(stderr, "(one) Duplicated %d in {%d, %d}\n", connection.id, i, j);
                    return false;
                }

                x.insert(connection.id);
            }
        }
    }

    return true;
}

bool checkTwo(const Solution &s) {
    for (int i = 0; i < s.spectrums.size(); i++) {
        for (int j = 0; j < s.spectrums[i].channels.size(); j++) {
            set<int> x;
            for (const Connection &connection : s.spectrums[i].channels[j].connections) {
                if (x.find(connection.id) != x.end()) {
                    fprintf(stderr, "(two) Duplicated %d in {%d, %d}\n", connection.id, i, j);
                    return false;
                }

                x.insert(connection.id);
            }
        }
    }

    return true;
}

Solution VNS(FILE **solutionFile, Solution initSol) {
    Solution delta = convert20MHz(initSol);
    Solution rep = multipleRepresentation(delta);

    setDP(rep);
    double retOF = calcDP(rep);

    Solution explicitSol = explicitSolution(rep);
    initSol.printSolution();
//  assert(explicitSol.totalThroughput >= initSol.totalThroughput);
    assert(double_equals(retOF, explicitSol.totalThroughput));
    assert(checkOne(explicitSol));

    Solution star = explicitSol;
    Solution localMax = delta;

    double _FO_localMax = retOF;
    double _FO_delta = retOF;
    double _FO_star = initSol.totalThroughput;

    int K_MUL = max(1, nConnections / 100);
    int K_MAX = 50;
    startTime = clock();
    while (!stop()) {
        int k = 1;
        while (k <= K_MAX && !stop()) {
            delta = localMax;
            _FO_delta = _FO_localMax;

            if (rng.randInt(1)) { //AddDrop
                K_AddDrop(delta, _FO_delta, k * K_MUL);
                fixChannels(delta, _FO_delta);
            } else { //Reinsert
                K_RemoveAndInserts(delta, _FO_delta, k * K_MUL);
                fixChannels(delta, _FO_delta);
            }

            Solution multiple = multipleRepresentation(delta);
            setDP(multiple);
            assert(checkTwo(multiple));
            _FO_delta = calcDP(multiple);
//      printf("hmm %.3lf\n", _FO_delta);

//      puts("entrei");
            explicitSol = newVNS_Reinsert(multiple, _FO_delta, delta);
//      puts("e sai");
            fixChannels(explicitSol, _FO_delta);
            delta = convert20MHz(explicitSol);

//      printf("(1) %.3lf %.3lf\n", explicitSol.totalThroughput, localMax.totalThroughput);
            if (_FO_delta > _FO_localMax) {
                k = 1;
                _FO_localMax = _FO_delta;
                localMax = delta;
            } else {
                k++;
            }
//      printf("(2) %.3lf %.3lf\n", explicitSol.totalThroughput, star.totalThroughput);
            if (_FO_localMax > _FO_star) {
//      if (_FO_localMax > star.totalThroughput) {
#ifndef DEBUG_CLION
                fprintf(*solutionFile, "->> %.3lf %.3lf %lf\n", explicitSol.totalThroughput, star.totalThroughput,
                        (((double) (clock() - startTime)) / CLOCKS_PER_SEC));
#endif
                _FO_star = _FO_localMax;
                star = explicitSol;
            }
        }
    }

    assert(checkOne(star));
    return star;
}

void init(int argc, char **argv, FILE **solutionFile = nullptr, FILE **objectivesFile = nullptr) {
#ifdef DEBUG_CLION //TODO: remind to remove the MACRO before real tests
    puts("============== WITH DEBUG ==============");
    freopen("/Users/joaquimnt_/git/vrbspheuristics/Instancias/D250x250/U_2048/U_2048_1.txt", "r", stdin);

    maximumTime = 10;
#else
    if (argc != 5) {
      fprintf(stderr, "wrong arguments. Provided %d, Must be: stdin, solutionFile, objectiveFile, timeLimit\n", argc);
      exit(13);
    }

    const string openingFile(argv[1]);
    if (!openingFile.empty()) {
      fprintf(stderr, "trying to open input file %s\n", openingFile.c_str());
      freopen(openingFile.c_str(), "r", stdin);
    }

    *solutionFile = fopen(argv[2], "a");
    if (solutionFile == NULL) {
      fprintf(stderr, "error opening solutionFile file\n");
      exit(13);
    }

    *objectivesFile = fopen(argv[3], "a");
    if (objectivesFile == NULL) {
      fprintf(stderr, "error opening objectivesFile file\n");
      exit(13);
    }

    maximumTime = stoi(argv[4]);
#endif

    if (stdin == nullptr) {
        fprintf(stderr, "error opening input file (stdin)\n");
        exit(13);
    }

    loadData();
    fprintf(stdout, "will execute for %lf seconds\n", maximumTime);
}

int main(int argc, char *argv[]) {
    FILE *solutionFile = nullptr, *objectivesFile = nullptr;
    init(argc, argv, &solutionFile, &objectivesFile);

    Solution aux = createSolution();
    Solution ans = VNS(&solutionFile, aux);

    printf("%.3lf ==> %.3lf\n", aux.totalThroughput, ans.totalThroughput);

#ifdef DEBUG_CLION
//  aux.printSolution();
    ans.printSolution();
#else
    if (solutionFile != nullptr) {
      ans.printSolution(solutionFile);
    } else {
      fprintf(stderr, "solutionFile is null!\n");
      exit(13);
    }

    if (objectivesFile != nullptr) {
      fprintf(objectivesFile, "%lf\n", ans.totalThroughput);
    } else {
      fprintf(stderr, "objectivesFiles is null!\n");
      exit(13);
    }

    fclose(solutionFile);
    fclose(objectivesFile);
#endif
    return 0;
}