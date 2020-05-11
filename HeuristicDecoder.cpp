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
    double conn_SINR = (powerSender / pow(distanceMatrix[conn.id][conn.id], alfa)) / (conn.interference + noise);
    conn.SINR = conn_SINR;

    while (mcs + 1 <= maxDataRate && conn_SINR > SINR[mcs + 1][bwIdx(bandwidth)])
      mcs++;


    conn.throughput = dataRates[mcs][bwIdx(bandwidth)];
  }

  return conn.throughput;
}

Channel insertInChannel(const Channel &channel, int idConn) {
  Channel newChannel = channel;
  Connection conn(idConn, 0.0, 0.0);

  for (Connection &connection : newChannel.connections) {
    connection.interference += interferenceMatrix[conn.id][connection.id];
    conn.interference += interferenceMatrix[connection.id][conn.id];
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
    conn.interference -= interferenceMatrix[idConn][conn.id];
    computeConnectionThroughput(conn, newChannel.bandwidth);
    newChannel.throughput += conn.throughput;
  }

  return newChannel;
}


double computeThroughput(Solution &curr, bool force) {
  double OF = 0.0;

  for (int s = 0; s < curr.spectrums.size(); s++) {
    for (int c = 0; c < curr.spectrums[s].channels.size(); c++) {
      double &chThroughput = curr.spectrums[s].channels[c].throughput;
      chThroughput = 0.0;
      for (Connection &conn : curr.spectrums[s].channels[c].connections) {
        conn.interference = 0.0;
        conn.throughput = 0.0;
        for (Connection &otherConn : curr.spectrums[s].channels[c].connections) {
          conn.interference += interferenceMatrix[otherConn.id][conn.id];
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

Solution split(Solution newSol, ii where) {
  Channel toSplit = newSol.spectrums[where.first].channels[where.second];
  Channel child1(toSplit.bandwidth / 2), child2(toSplit.bandwidth / 2);

  swap(newSol.spectrums[where.first].channels[where.second], newSol.spectrums[where.first].channels.back());
  newSol.spectrums[where.first].channels.pop_back();

  vector<Connection> scheduled_conn(toSplit.connections);

  if (scheduled_conn.size() == 1) {
    if (rng.randInt(1)) {
      child1.connections.push_back(scheduled_conn[0]);
    } else {
      child2.connections.push_back(scheduled_conn[0]);
    }
  } else if (scheduled_conn.size() == 2) {
    if (rng.randInt(1)) {
      child1.connections.push_back(scheduled_conn[0]);
      child2.connections.push_back(scheduled_conn[1]);
    } else {
      child1.connections.push_back(scheduled_conn[1]);
      child2.connections.push_back(scheduled_conn[0]);
    }
  } else {//Three or more connections
    int a = -1, b = -1;
    double a_in = 0.0, b_in = 0.0;
    for (int i = 0; i < scheduled_conn.size(); i++) {
      if (scheduled_conn[i].interference >= a_in) {
        a = i;
        a_in = scheduled_conn[i].interference;
      }
    }

    for (int i = 0; i < scheduled_conn.size(); i++) {
      if (i == a)
        continue;

      if (scheduled_conn[i].interference >= b_in) {
        b = i;
        b_in = scheduled_conn[i].interference;
      }
    }

    assert(a >= 0 && b >= 0 && a != b);

    if (rng.randInt(1)) {
      child1.connections.emplace_back(scheduled_conn[a]);
      child2.connections.emplace_back(scheduled_conn[b]);
    } else {
      child1.connections.emplace_back(scheduled_conn[b]);
      child2.connections.emplace_back(scheduled_conn[a]);
    }

    for (int i = 2; i < scheduled_conn.size(); i++) {
      if (rng.randInt(1)) {
        child1.connections.emplace_back(scheduled_conn[i]);
      } else {
        child2.connections.emplace_back(scheduled_conn[i]);
      }
    }
  }

  vector<Channel> aux = {child1, child2};
  insertInSpectrum(newSol, aux, where.first);
  return newSol;
}

Solution createSolution() {
  Solution ret(initConfiguration, 0.0, true);

  vector<int> links;
  for (int i = 0; i < nConnections; i++)
    links.emplace_back(i);

  shuffle(links.begin(), links.end(), whatever);

  vector<int> notInserted;
  Solution retCopy(ret);
  for (int conn : links) {
    Solution copySolution(retCopy);

    for (int s = 0; s < copySolution.spectrums.size(); s++) {
      for (int c = 0; c < copySolution.spectrums[s].channels.size(); c++) {
        Solution copySolution2(copySolution);
        Channel newChannel = insertInChannel(copySolution2.spectrums[s].channels[c], conn);

        swap(newChannel, copySolution2.spectrums[s].channels[c]);
        computeThroughput(copySolution2);

        Solution candidate1 = copySolution2;
        Solution candidate2;

        if (copySolution2.spectrums[s].channels[c].bandwidth >= 40)
          candidate2 = split(candidate1, {s, c});

        if (candidate2 > candidate1) { //Better to split
          candidate1 = candidate2;
        }

        if (candidate1 > retCopy) {
          retCopy = candidate1;
        }
      }
    }
  }

  ret = retCopy;
  Channel zeroChannel(0.0, 0.0, 0, vector<Connection>());
  for (const int x : notInserted) {
    if (!x) {
      zeroChannel.connections.emplace_back(Connection(x));
    }
  }

  ret.spectrums.emplace_back(Spectrum(0, 0, {zeroChannel}));
  return ret;
}

Connection::Connection(int id, double throughput, double interference) : id(id),
                                                                         throughput(throughput),
                                                                         interference(interference) {}

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
}

void Solution::printSolution(FILE *file) {
  if (file == nullptr) {
    file = stdout;
  }

  fprintf(file, "%lf\n", totalThroughput);
  for (int i = 0; i < spectrums.size(); i++) {
    fprintf(file, "In spec %d:\n", i);
    for (int j = 0; j < spectrums[i].channels.size(); j++) {
      fprintf(file, "  In channel %d (%d MHz): ", j, spectrums[i].channels[j].bandwidth);
      for (Connection &conn : spectrums[i].channels[j].connections) {
        fprintf(file, "{%d, %.10lf, %.10lf, %lf} ", conn.id, conn.interference, conn.SINR, conn.throughput);
      }
      fprintf(file, "\n");
    }
    fprintf(file, "\n");
  }
}