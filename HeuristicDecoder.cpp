//
// Created by José Joaquim on 03/03/20.
//

#include "HeuristicDecoder.h"

using namespace std;

const int X_c = 0;
const int Y_c = 1;
const double EPS = 1e-9;

int nConnections, nSpectrums;
double dataRates[10][4];
double distanceMatrix[MAX_CONN][MAX_CONN], interferenceMatrix[MAX_CONN][MAX_CONN];
double senders[MAX_CONN][2], receivers[MAX_CONN][2];
std::vector<std::vector<double>> SINR;
double powerSender, alfa, noise, ttm;
std::vector<Spectrum> initConfiguration;
std::random_device rd;
auto whatever = std::default_random_engine{rd()};
MTRand rng;
std::pair<int, int> zeroChannel;

inline double distance(double X_si, double Y_si, double X_ri, double Y_ri) {
  return hypot((X_si - X_ri), (Y_si - Y_ri));
}

void distanceAndInterference() {
  for (int i = 0; i < nConnections; i++) {
    double X_si = receivers[i][X_c];
    double Y_si = receivers[i][Y_c];


    for (int j = 0; j < nConnections; j++) {

      double X_rj = senders[j][X_c];
      double Y_rj = senders[j][Y_c];

      double dist = distance(X_si, Y_si, X_rj, Y_rj);

      distanceMatrix[i][j] = dist;

      if (i == j) {
        interferenceMatrix[i][j] = 0.0;
      } else {
        double value = (dist != 0.0) ? powerSender / pow(dist, alfa) : 1e9;
        interferenceMatrix[i][j] = value;
      }
    }
  }
}

double convertDBMToMW(double _value) {
  double result = 0.0, b;

  b = _value / 10.0;// dBm dividido por 10
  result = pow(10.0, b);//Converte de DBm para mW

  return result;
}

void convertTableToMW(const vector<vector<double> > &_SINR, vector<vector<double> > &_SINR_Mw) {
  double result, b;
  for (int i = 0; i < _SINR_Mw.size(); i++) {
    for (int j = 0; j < _SINR_Mw[i].size(); j++) {

      if (_SINR[i][j] != 0) {
        b = _SINR[i][j] / 10.0;// dBm divided by 10
        result = pow(10.0, b);//Convert DBM to mW

        _SINR_Mw[i][j] = result;
      } else {
        _SINR_Mw[i][j] = 0;
      }
    }
  }
}

void initSpectrums() {
  for (Spectrum &sp : initConfiguration) {
    while (sp.maxFrequency - sp.usedFrequency > 0) {
      int bw = 160;
      while (bw > (sp.maxFrequency - sp.usedFrequency) && bw > 20) {
        bw /= 2;
      }

      if (bw <= (sp.maxFrequency - sp.usedFrequency)) {
        sp.usedFrequency += bw;
        sp.channels.emplace_back(0.0, 0.0, bw, vector<Connection>());
      }
    }

    assert(sp.maxFrequency - sp.usedFrequency >= 0);
  }

  zeroChannel = {initConfiguration.size(), 0};
}

void loadData() {
  SINR.assign(10, vector<double>(4, 0));
  int aux1;
  scanf("%d", &aux1);
  scanf("%d %lf %lf %lf %lf %d", &nConnections, &ttm, &alfa, &noise, &powerSender, &nSpectrums);

  for (int i = 0; i < nSpectrums; i++) {
    int s;
    scanf("%d", &s);
    initConfiguration.emplace_back(s, 0, vector<Channel>());
  }

  if (noise != 0) {
    noise = convertDBMToMW(noise);
  }

  for (int i = 0; i < 10; i++) {
    for (int j = 0; j < 4; j++) {
      scanf("%lf", &dataRates[i][j]);
    }
  }

  for (int i = 0; i < 10; i++) {
    for (int j = 0; j < 4; j++) {
      scanf("%lf", &SINR[i][j]);
    }
  }

  convertTableToMW(SINR, SINR);

  for (int i = 0; i < nConnections; i++) {
    double x, y;
    scanf("%lf %lf", &x, &y);
    receivers[i][X_c] = x;
    receivers[i][Y_c] = y;
  }

  for (int i = 0; i < nConnections; i++) {
    double x, y;
    scanf("%lf %lf", &x, &y);
    senders[i][X_c] = x;
    senders[i][Y_c] = y;
  }

  memset(interferenceMatrix, 0, sizeof interferenceMatrix);
  memset(distanceMatrix, 0, sizeof distanceMatrix);

  initSpectrums();
  distanceAndInterference();

//  printf("===== SINR ======\n");
//  for (int i = 0; i < 10; i++) {
//    for (int j = 0; j < 4; j++) {
//      printf("%7.5lf ", SINR[i][j]);
//    }
//    puts("");
//  }
//
//  printf("===== DR ======\n");
//  for (int i = 0; i < 10; i++) {
//    for (int j = 0; j < 4; j++) {
//      printf("%7.5lf ", dataRates[i][j]);
//    }
//    puts("");
//  }

//  for (int i = 0; i < 10; i++) {
//    for (int j = 0; j < 10; j++) {
//      printf("%.10lf ", interferenceMatrix[i][j]);
//    }
//    puts("");
//  }
}

int bwIdx(int bw) {
  if (bw == 40) {
    return 1;
  } else if (bw == 80) {
    return 2;
  } else if (bw == 160) {
    return 3;
  }
  return 0;
}

bool double_equals(double a, double b, double epsilon) {
  return std::abs(a - b) < epsilon;
}

bool operator>(const Solution &o1, const Solution &o2) {
  return operator<(o2, o1);
}

bool operator<(const Solution &o1, const Solution &o2) {
  assert(o1.throughputFlag && o2.throughputFlag);
  return o1.totalThroughput < o2.totalThroughput;
}

double computeConnectionThroughput(Connection &conn, int bandwidth, bool force) {
  int mcs = -1;
  int maxDataRate = bandwidth == 20 ? 8 : 9;

  if (double_equals(conn.interference, 0.0)) {
    mcs = maxDataRate;
    conn.throughput = dataRates[mcs][bwIdx(bandwidth)];
  } else {
    double conn_SINR = (powerSender / pow(conn.distanceSR, alfa)) / (conn.interference + noise);
    conn.SINR = conn_SINR;

    while (mcs + 1 <= maxDataRate && conn_SINR > SINR[mcs + 1][bwIdx(bandwidth)])
      mcs++;

    if (force) {
      printf("   ->> SINR EH %.10lf, MCS %d, BW %d\n", conn_SINR, mcs, bandwidth);
    }


    conn.throughput = dataRates[mcs][bwIdx(bandwidth)];
  }

  return conn.throughput;
}

void testConnection(Connection &conn, int bandwidth, bool force = true) {
  printf("conn %d:\n", conn.id);
  printf("  interference: %.10lf\n", conn.interference);
  printf("  bw: %d\n", bandwidth);
  printf("  distanceSR: %.10lf\n", conn.distanceSR);
  computeConnectionThroughput(conn, bandwidth, true);
}

Channel insertInChannel(Channel newChannel, int idConn) {
  Connection conn(idConn, 0.0, 0.0, distanceMatrix[idConn][idConn]);

  for (Connection &connection : newChannel.connections) {
//    connection.interference += interferenceMatrix[conn.id][connection.id];
//    conn.interference += interferenceMatrix[connection.id][conn.id];
    connection.interference += interferenceMatrix[connection.id][conn.id];
    conn.interference += interferenceMatrix[conn.id][connection.id];
  }

  newChannel.connections.emplace_back(conn);
  newChannel.throughput = 0.0;
  for (Connection &connection : newChannel.connections) {
    computeConnectionThroughput(connection, newChannel.bandwidth);
    newChannel.throughput += connection.throughput;
  }

  return newChannel;
}

Channel deleteFromChannel(const Channel &channel, int idConn) {
  Channel newChannel(channel.bandwidth);

  for (const Connection &conn : channel.connections) {
    if (conn.id != idConn) {
      newChannel.connections.emplace_back(conn);
    }
  }

  newChannel.throughput = 0.0;
  for (Connection &conn : newChannel.connections) {
    conn.interference -= interferenceMatrix[conn.id][idConn];
    computeConnectionThroughput(conn, newChannel.bandwidth);
    newChannel.throughput += conn.throughput;
  }

  return newChannel;
}

double computeThroughput(Solution &curr, bool force) {
  double OF = 0.0;

  for (int s = 0; s < curr.spectrums.size(); s++) {
    for (int c = 0; c < curr.spectrums[s].channels.size(); c++) {
      if (make_pair(s, c) == zeroChannel)
        break;

      double &chThroughput = curr.spectrums[s].channels[c].throughput;
      chThroughput = 0.0;
      for (Connection &conn : curr.spectrums[s].channels[c].connections) {
        conn.interference = 0.0;
        conn.throughput = 0.0;
        for (Connection &otherConn : curr.spectrums[s].channels[c].connections) {
          if (conn.id == 0 && force) {
//            printf("vou somar %.10lf\n", interferenceMatrix[conn.id][otherConn.id]);
          }
          conn.interference += interferenceMatrix[conn.id][otherConn.id];
        }
        chThroughput += computeConnectionThroughput(conn, curr.spectrums[s].channels[c].bandwidth, force);
      }
      OF += chThroughput;
    }
  }

  curr.totalThroughput = OF;
  curr.throughputFlag = true;
  return OF;
}

void insertInSpectrum(Solution &sol, vector<Channel> &channels, int specId) {
  sol.throughputFlag = false;
  for (const Channel &ch : channels) {
    sol.spectrums[specId].channels.emplace_back(ch);
  }

  computeThroughput(sol);
}

void rawInsert(Solution &sol, int conn, ii where) {
  sol.spectrums[where.first].channels[where.second].connections.emplace_back(Connection(conn));
}

//double testConnThroughput(int idConn, Channel newChannel) {
//  double ret = 0.0;
//  Connection conn(idConn, 0.0, 0.0, distanceMatrix[idConn][idConn]);
//
//  for (Connection &connection : newChannel.connections) {
//    connection.interference += interferenceMatrix[conn.id][connection.id];
//    conn.interference += interferenceMatrix[connection.id][conn.id];
//  }
//
//  newChannel.connections.emplace_back(conn);
//  newChannel.throughput = 0.0;
//  for (Connection &connection : newChannel.connections) {
//    computeConnectionThroughput(connection, newChannel.bandwidth);
//    newChannel.throughput += connection.throughput;
//  }
//
//  return ret;
//}

void computeChannelThroughput(vector<Channel> &channels) {
  for (Channel &channel : channels) {
    double of = 0.0;
    for (Connection &conn : channel.connections) {
      of += computeConnectionThroughput(conn, channel.bandwidth);
    }
    channel.throughput = of;
  }
}

int GLOBAL_CONT;

double some(const Connection &c1, const Connection &c2) {
  return c1.throughput + c2.throughput;
}

vector<Channel> split(Channel toSplit, int spectrum) {
  int newBw = toSplit.bandwidth / 2;
  vector<Channel> ret;
  ret.emplace_back(Channel(0.0, 0.0, newBw, vector<Connection>()));
  ret.emplace_back(Channel(0.0, 0.0, newBw, vector<Connection>()));

  vector<Connection> connToInsert = toSplit.connections;
  sort(connToInsert.rbegin(), connToInsert.rend());

  double bestThroughputSoFar = 0.0;
  for (int i = 0; i < connToInsert.size(); i++) {
    double bestThroughputIteration = 0.0, testando = 0.0;
    Channel bestChannel(0.0, 0.0, newBw, vector<Connection>());
    int bestIdx = -1;

    for (int c = 0; c < ret.size(); c++) {
      Channel inserted = insertInChannel(ret[c], connToInsert[i].id);
      double resultThroughput = inserted.throughput;

      double g_throughput = (bestThroughputSoFar - ret[c].throughput) + resultThroughput;
      testando = max(testando, g_throughput);
      bool one = g_throughput > bestThroughputIteration;
      bool two = double_equals(g_throughput, bestThroughputIteration) &&
                 (inserted.connections.size() < bestChannel.connections.size());

      if (one || two) {
        bestThroughputIteration = g_throughput;
        bestChannel = inserted;
        bestIdx = c;
      }
    }

    assert(bestIdx >= 0);
    ret[bestIdx] = bestChannel;
    bestThroughputSoFar = bestThroughputIteration;
  }
  computeChannelThroughput(ret);

//  if (GLOBAL_CONT == 283) {
//    puts(" SPLIT ");
//    double l_throughput = 0.0;
//    for (int c = 0; c < ret.size(); c++) {
//      double value = 0.0;
//      for (Connection &conn : ret[c].connections) {
//        value += conn.throughput;
//      }
//      printf("      no subcanal %d: %.3lf\n", c, value);
//      for (Connection &conn : ret[c].connections) {
//        printf("{%d, %.10lf, %.10lf, %lf} ", conn.id, conn.interference, conn.SINR, conn.throughput);
//      }
//      puts("");
//    }
//    printf("    THROUGHPUT RESULTOU EM %.3lf %.3lf\n", l_throughput, bestThroughputSoFar);
////    exit(44);
//  }

  assert(double_equals(ret[0].throughput + ret[1].throughput, bestThroughputSoFar));
  return ret;
}

void printChannel(const Channel &chan, ii where = {-1, -1}) {
  printf("Atual estado do canal {%d, %d} (%d MHz):\n", where.first, where.second, chan.bandwidth);
  printf("   Throughput: %.4lf\n", chan.throughput);
  printf("   Conexoes:");
  double sum = 0.0;
  for (const Connection &conn : chan.connections) {
    sum += conn.throughput;
    printf(" {%d, %.4lf, %.4lf} ", conn.id, conn.distanceSR, conn.throughput);
  }
  puts("");

  assert(double_equals(sum, chan.throughput));
}

//void test_debug(const Solution &sol, int conn) {
////  printf("Testando insercao de %d\n", conn);
//  double bestFound = 0.0;
//  for (int s = 0; s < sol.spectrums.size(); s++) {
//    for (int c = 0; c < sol.spectrums[s].channels.size(); c++) {
////      printChannel(sol.spectrums[s].channels[c], {s, c});
//      Channel aux = insertInChannel(sol.spectrums[s].channels[c], conn);
//      bestFound = max(bestFound, sol.totalThroughput - sol.spectrums[s].channels[c].throughput + aux.throughput);
////      printf("===== depois de inserir =====\n");
////      printChannel(aux, {s, c});
//      if (sol.spectrums[s].channels[c].bandwidth >= 40) {
////        printf("===== A divisao do canal =====\n");
//        vector<Channel> channelSplit = split(sol.spectrums[s].channels[c]);
////        printChannel(channelSplit[0]);
////        printChannel(channelSplit[1]);
//
//        Channel splitAux0 = insertInChannel(channelSplit[0], conn);
//        Channel splitAux1 = insertInChannel(channelSplit[1], conn);
////        puts(" inserindo em [0]:");
////        printChannel(splitAux0);
////        puts(" inserindo em [1]");
////        printChannel(splitAux1);
//        bestFound = max(bestFound,
//                        sol.totalThroughput - sol.spectrums[s].channels[c].throughput + channelSplit[0].throughput +
//                        splitAux1.throughput);
//        bestFound = max(bestFound,
//                        sol.totalThroughput - sol.spectrums[s].channels[c].throughput + channelSplit[1].throughput +
//                        splitAux0.throughput);
//      }
//
////      puts("---------- fim do canal");
//    }
//  }
//
////  printf("    ->>> MAIOR THROUGHPUT FOI %lf <<<-\n", bestFound);
//  if (bestFound > sol.totalThroughput) {
//    printf("OPAAA %lf %lf\n", bestFound, sol.totalThroughput);
//  }
//  assert(bestFound <= sol.totalThroughput);
//}

Solution createSolution() {
  Solution ret(initConfiguration, 0.0, true);
  vector<int> links;
  for (int i = 0; i < nConnections; i++)
    links.emplace_back(i);

  shuffle(links.begin(), links.end(), whatever);

  Solution retCopy(ret);
  vector<int> zeroConn;
  double bestThroughputSoFar = 0.0;
  for (int conn : links) {
//  for (int conn = 0; conn < 800; conn++) {
    GLOBAL_CONT = conn;
//    if (GLOBAL_CONT == 283)
//      printf("na conexao %d, throughput atual %.3lf\n", GLOBAL_CONT, bestThroughputSoFar);
    double bestThroughputIteration = bestThroughputSoFar;
    vector<Channel> bestSplitChannels;
    Channel bestChannel(0.0, 0.0, 0, vector<Connection>());
    Channel oldChannel(0.0, 0.0, 0, vector<Connection>());
    ii where = {-1, -1};
    bool isSplit = false, inserted = false;
    int bestBandwidthSplit = -1, bandwidthSplit = -1;


    for (int s = 0; s < retCopy.spectrums.size(); s++) {
      for (int c = 0; c < retCopy.spectrums[s].channels.size(); c++) {
//        if (GLOBAL_CONT == 283) {
//          printf("tentando {%d, %d}: %.3lf\n", s, c, retCopy.spectrums[s].channels[c].throughput);
//        }
        const Channel &currentChannel = retCopy.spectrums[s].channels[c];
        bandwidthSplit = -1;

        Channel channelInserted = insertInChannel(currentChannel, conn);

//        if (GLOBAL_CONT == 283) {
//          printf(" inserir canal completo %.3lf\n", channelInserted.throughput);
//        }
        vector<Channel> channelSplit;

        if (channelInserted.bandwidth >= 40) {
//          channelSplit = split(channelInserted);
          channelSplit = split(currentChannel, s);

          Channel split0 = insertInChannel(channelSplit[0], conn);
          Channel split1 = insertInChannel(channelSplit[1], conn);

//          if (GLOBAL_CONT == 283) {
//            printf("inserindo em split0 (%d conns): %.3lf --> total %.3lf\n", split0.connections.size(),
//                   split0.throughput, split0.throughput + channelSplit[1].throughput);
//            printf("inserindo em split1 (%d conns): %.3lf --> total %.3lf\n", split1.connections.size(),
//                   split1.throughput, split1.throughput + channelSplit[0].throughput);
//          }
//          bool cond01 = split0.throughput > split1.throughput;
//          bool cond02 = double_equals(split0.throughput, split1.throughput) &&
//                        (split0.connections.size() < split1.connections.size());
//
//          bool cond11 = split1.throughput > split0.throughput;
//          bool cond12 = double_equals(split0.throughput, split1.throughput) &&
//                        (split1.connections.size() < split0.connections.size());
//
//          if (cond01 || cond02) {
//            channelSplit[0] = split0;
//          } else if (cond11 || cond12) {
//            channelSplit[1] = split1;
//          } else {
//            channelSplit[0] = split0;
//          }

          int op = 0;
          double conf1 = split0.throughput + channelSplit[1].throughput;
          double conf2 = split1.throughput + channelSplit[0].throughput;
          if (conf1 > conf2) {
            channelSplit[0] = split0;
          } else if (conf2 > conf1) {
            op = 1;
            channelSplit[1] = split1;
          } else if (split0.connections.size() < split1.connections.size()) {
            op = 2;
            channelSplit[0] = split0;
          } else {
            op = 3;
            channelSplit[1] = split1;
          }

//          printf(" WINNING OP %d\n", op);
        }

        double ammount = channelInserted.throughput;
        bool betterSplit = false;
        if (!channelSplit.empty()) {
          ammount = max(channelSplit[0].throughput + channelSplit[1].throughput, ammount);
          betterSplit = (channelSplit[0].throughput + channelSplit[1].throughput) > channelInserted.throughput;
        }

        bandwidthSplit = betterSplit ? currentChannel.bandwidth / 2 : currentChannel.bandwidth;
        double totalThroughputIteration = bestThroughputSoFar - currentChannel.throughput + ammount;

//        if (GLOBAL_CONT == 283) {
//          printf("%.3lf %.3lf %.3lf\n", totalThroughputIteration, bestThroughputIteration, ammount);
//        }
        bool cond1 = totalThroughputIteration > bestThroughputIteration;
        bool cond2 =
                double_equals(totalThroughputIteration, bestThroughputIteration) && bandwidthSplit > bestBandwidthSplit;
        bool cond3 = double_equals(totalThroughputIteration, bestThroughputIteration) &&
                     (bandwidthSplit == bestBandwidthSplit) && betterSplit == false && isSplit == true;

//        if (GLOBAL_CONT == 283) {
//          printf("%d %d %d\n", cond1, cond2, cond3);
//        }

        if (cond1 || cond2 || cond3) {
          bestThroughputIteration = totalThroughputIteration;
          inserted = true;
          where = {s, c};

          bestBandwidthSplit = bandwidthSplit;
          if (betterSplit) {
            isSplit = true;
            bestSplitChannels = channelSplit;
          } else {
            isSplit = false;
            bestChannel = channelInserted;
          }
//          if (!channelSplit.empty()) {
//            if (betterSplit) {
//              isSplit = true;
//              bestSplitChannels = channelSplit;
//            } else {
//              isSplit = false;
//              bestChannel = channelInserted;
//            }
//          } else {
//            bestChannel = channelInserted;
//          }
        }
      }
    }

    double aux = bestThroughputSoFar;
    bestThroughputSoFar = bestThroughputIteration;

    if (inserted) {
//      printf("Coloquei conn %d em {%d, %d}, throughput %.3lf -> %.3lf\n", conn, where.first, where.second, aux,
//             bestThroughputSoFar);
      if (isSplit) {
        retCopy.spectrums[where.first].channels.erase(retCopy.spectrums[where.first].channels.begin() + where.second);
        retCopy.spectrums[where.first].channels.insert(retCopy.spectrums[where.first].channels.begin() + where.second,
                                                       bestSplitChannels.begin(), bestSplitChannels.end());
//        printf("      -> com split. Throughput (1) %.3lf, (2) %.3lf ==> %.3lf\n", bestSplitChannels[0].throughput,
//               bestSplitChannels[1].throughput, bestSplitChannels[0].throughput + bestSplitChannels[1].throughput);
      } else {
//        printf("      -> sem split.\n");
        retCopy.spectrums[where.first].channels[where.second] = bestChannel;
      }
      computeThroughput(retCopy);

      if (!double_equals(retCopy.totalThroughput, bestThroughputSoFar)) {
        printf("vixe %.3lf %.3lf\n", retCopy.totalThroughput, bestThroughputSoFar);
      }
      assert(double_equals(retCopy.totalThroughput, bestThroughputSoFar));
    } else {
//      printf("====== Nao inseri conexao %d ======\n", conn);
      zeroConn.push_back(conn);
    }

//    if (GLOBAL_CONT == 283)
//      exit(44);
  }

  computeThroughput(retCopy);
  assert(double_equals(bestThroughputSoFar, retCopy.totalThroughput));

  Channel aux(0.0, 0.0, 0, vector<Connection>());
  retCopy.spectrums.emplace_back(Spectrum(0, 0, {aux}));
  for (int i = 0; i < zeroConn.size(); i++) {
    rawInsert(retCopy, zeroConn[i], {zeroChannel.first, zeroChannel.second});
  }
//  printf("nao coloquei %lu conexoes\n", zeroConn.size());

  ret = retCopy;
  return ret;
}

Connection::Connection(int id, double throughput, double interference, double distanceSR) : id(id),
                                                                                            throughput(throughput),
                                                                                            interference(interference),
                                                                                            distanceSR(distanceSR) {}

Channel::Channel(double throughput, double interference, int bandwidth, const vector<Connection> &connections)
        : throughput(throughput), interference(interference), bandwidth(bandwidth), connections(connections) {}

Channel::Channel(int bandwidth) : bandwidth(bandwidth) {
  interference = 0.0;
  throughput = 0.0;
  connections = vector<Connection>();
}

Spectrum::Spectrum(int maxFrequency, int usedFrequency, const vector<Channel> &channels) : maxFrequency(
        maxFrequency), usedFrequency(usedFrequency), channels(channels) {}

Solution::Solution(const vector<Spectrum> &spectrums, double total, bool flag) : spectrums(spectrums),
                                                                                 totalThroughput(total),
                                                                                 throughputFlag(flag) {}

Solution::Solution() {
  throughputFlag = true;
  totalThroughput = 0.0;
}

Connection::Connection(int id) : id(id) {
  interference = 0.0;
  throughput = 0.0;
  distanceSR = distanceMatrix[id][id];
}

void Solution::printSolution(FILE *file) {
  if (file == nullptr) {
    file = stdout;
  }

  int cont = 0;
  fprintf(file, "%lf\n", totalThroughput);
  for (int i = 0; i < spectrums.size(); i++) {
    if (i == 3)
      break;

    fprintf(file, "In spec %d:\n", i);
    for (int j = 0; j < spectrums[i].channels.size(); j++) {
      fprintf(file, "  In channel %d (%d MHz): ", j, spectrums[i].channels[j].bandwidth);
      for (Connection &conn : spectrums[i].channels[j].connections) {
        cont++;
        fprintf(file, "{%d, %.10lf, %.10lf, %lf} ", conn.id, conn.interference, conn.SINR, conn.throughput);
      }
      fprintf(file, "\n");
    }
    fprintf(file, "\n");
  }

  fprintf(file, "TOTAL OF %d CONNECTIONS\n", cont);
}