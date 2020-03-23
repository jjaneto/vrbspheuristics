//
// Created by José Joaquim on 03/03/20.
//

#include "HeuristicDecoder.h"

const int X_c = 0;
const int Y_c = 1;
const double EPS = 1e-9;

int nConnections;
int L;
double dataRates[10][4];
double distanceMatrix[2048][2048], interferenceMatrix[2048][2048];
double senders[2048][2], receivers[2048][2];
vector<vector<double>> SINR;
double powerSender, alfa, noise, ttm;
map<int, vector<int>> chToLinks;

int overlap[45][45] = {{1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0},
                       {0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0},
                       {0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0},
                       {0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0},
                       {0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0},
                       {0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0},
                       {0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0},
                       {0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0},
                       {0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1},
                       {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1},
                       {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1},
                       {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1},
                       {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1},
                       {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1},
                       {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1},
                       {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1},
                       {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0},
                       {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0},
                       {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0},
                       {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0},
                       {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0},
                       {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0},
                       {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0},
                       {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0},
                       {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                       {1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0},
                       {0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0},
                       {0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0},
                       {0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0},
                       {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1},
                       {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1},
                       {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1},
                       {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1},
                       {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0},
                       {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0},
                       {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0},
                       {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0},
                       {1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0},
                       {0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0},
                       {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1},
                       {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1},
                       {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0},
                       {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0},
                       {1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0},
                       {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1}};


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

      double value = (dist != 0.0) ? powerSender / pow(dist, alfa) : 1e9;
      interferenceMatrix[i][j] = value;
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

inline void mapSplitChannels() {
  mapChtoCh[44] = {38, 39};
  mapChtoCh[38] = {26, 27};
  mapChtoCh[39] = {28, 29};
  mapChtoCh[26] = {1, 2};
  mapChtoCh[27] = {3, 4};
  mapChtoCh[28] = {5, 6};
  mapChtoCh[29] = {7, 8};
  mapChtoCh[45] = {40, 41};
  mapChtoCh[40] = {30, 31};
  mapChtoCh[41] = {32, 33};
  mapChtoCh[30] = {9, 10};
  mapChtoCh[31] = {11, 12};
  mapChtoCh[32] = {13, 14};
  mapChtoCh[33] = {15, 16};
  mapChtoCh[42] = {34, 35};
  mapChtoCh[34] = {17, 18};
  mapChtoCh[35] = {19, 20};
  mapChtoCh[43] = {36, 37};
  mapChtoCh[36] = {21, 22};
  mapChtoCh[37] = {23, 24};
}

void readFile() {
  SINR.assign(10, vector<double>(4, 0));
  double aux1;
  scanf("%lf", &aux1);
  scanf("%d %lf %lf %lf %lf %lf %lf %lf %lf", &nConnections, &ttm, &alfa, &noise, &powerSender, &aux1, &aux1,
        &aux1, &aux1);

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

  distanceAndInterference();
}

int whichBw(int ch) {
  if (ch >= 26 && ch <= 37)
    return 40;
  else if (ch >= 38 && ch <= 43)
    return 80;
  else if (ch == 44 || ch == 45)
    return 160;

  return 20;
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

bool double_equals(double a, double b, double epsilon = 0.000000001) {
  return std::abs(a - b) < epsilon;
}

Link::Link() {
  _idR = _idS = id = -1;
  ch = bw = -1;
  interference = SINR = 0.0;
  MCS = -1;
  distanceSenderReceiver = -1.0;
}

Link::Link(int id) : _idR(id), _idS(id), id(id) {
  ch = bw = -1;
  interference = SINR = 0.0;
  MCS = -1;
  distanceSenderReceiver = distanceMatrix[_idS][_idR];
}

Link::Link(const Link &x) {
  id = x.id;
  ch = x.ch;
  bw = x.bw;
  _idS = x._idS;
  _idR = x._idR;
  interference = x.interference;
  MCS = x.MCS;
  SINR = x.SINR;
  distanceSenderReceiver = x.distanceSenderReceiver;
}

void Link::operator=(const Link &x) {
  id = x.id;
  ch = x.ch;
  bw = x.bw;
  interference = x.interference;
  MCS = x.MCS;
  _idR = x._idR;
  _idS = x._idS;
  distanceSenderReceiver = x.distanceSenderReceiver;
  SINR = x.SINR;
}

int Link::getId() const {
  return this->id;
}

int Link::getChannel() const {
  return this->ch;
}

bool operator==(const Link &o1, const Link &o2) {
  if (o1._idR != o2._idR) {
    return false;
  }

  if (o1.id != o2.id) {
    return false;
  }

  if (o1.ch != o2.ch) {
    return false;
  }

  if (o1.bw != o2.bw) {
    return false;
  }

  if (o1.interference != o2.interference) {
    return false;
  }

  if (o1.SINR != o2.SINR) {
    return false;
  }

  if (o1.MCS != o2.MCS) {
    return false;
  }

  if (o1.distanceSenderReceiver != o2.distanceSenderReceiver) {
    return false;
  }

  return true;
}

void Link::setChannel(int ch) {
  this->ch = ch;
  this->bw = whichBw(ch);
}

void Link::printLink() const {
  printf("============= COMPARE =============\n");
  printf("link %d idS %d idR %d distance %.5lf ch %d, bw %d, interference %.10lf, SINR %.10lf, MCS %d\n", id, _idR,
         _idS, distanceSenderReceiver, ch, bw, interference, SINR, MCS);
  printf("=============== END ===============\n");
}

Solution::Solution(const Solution &o1) {
  objective = o1.objective;
//  scheduled_links = o1.scheduled_links;
  for (int i = 0; i < o1.scheduled_links.size(); i++) {
    scheduled_links.push_back(o1.scheduled_links[i]);
  }
}

Solution::Solution() {
  objective = 0.0;
}

//Solution &Solution::operator=(const Solution &o1) {
//  this->objective = o1.objective;
//  this->scheduled_links = o1.scheduled_links;
//  return *this;
//}

bool operator<(const Solution &o1, const Solution &o2) {
  return o1.objective < o2.objective;
}

bool operator>(const Solution &o1, const Solution &o2) {
  return o1.objective > o2.objective;
}

bool operator==(const Solution &o1, const Solution &o2) {

  if (o1.objective != o2.objective)
    return false;

  if (o1.scheduled_links.size() != o2.scheduled_links.size())
    return false;

  bool cond = true;

  const deque<Link> arr1 = o1.getScheduledLinks();
  const deque<Link> arr2 = o2.getScheduledLinks();

  for (int i = 0; i < int(arr1.size()) && cond; i++) {
    bool go = false;
    for (int j = 0; j < int(arr2.size()); j++) {
      if (arr1[i] == arr2[j]) {
        go = true;
      }
    }

    cond &= go;
  }

  return cond;
}

void Solution::computeInterference() {
  for (Link &u : scheduled_links) {
    u.interference = 0.0;
    u.SINR = 0.0;
    for (Link &v : scheduled_links) {

      if (u == v) {
        continue;
      }

      if (overlap[u.ch - 1][v.ch - 1]) {
        u.interference += interferenceMatrix[v._idS][u._idR];
      }
    }

    if (double_equals(u.interference, 0.0)) {
      u.SINR = 1e9;
    } else {
      u.SINR = (powerSender / pow(distanceMatrix[u._idS][u._idR], alfa)) / (u.interference + noise);
    }
  }
}

void Solution::computeObjective(bool show) {
  computeInterference();
  objective = 0.0;
  for (Link &x : scheduled_links) {
    int mxDataRate = (x.bw == 20) ? 9 : 10;
    bool go = false;

    for (int _mcs = 0; _mcs < mxDataRate; _mcs++) {
      if (show) {
        printf("will compare %.3lf and %.3lf (_mcs %d, bw %d, idx bw %d)\n", SINR[_mcs][bwIdx(x.bw)], x.SINR, _mcs,
               x.bw, bwIdx(x.bw));
      }
      if (SINR[_mcs][bwIdx(x.bw)] > x.SINR) {
        if (show) {
          printf("     uhu!\n");
        }
        x.MCS = _mcs - 1;
        go = true;

        if (x.MCS == -1) {
          x.MCS = 0;
        }

        break;
      } else {
        if (show) {
          printf("     ouch!\n");
        }
      }
    }

    if (!go) {
      if (show) {
        printf("enter here link %d mxDataRate %d\n", x.id, mxDataRate);
      }
      x.MCS = mxDataRate - 1;
    } else {
      if (show) {
        printf("NOT enter here link %d mxDataRate %d\n", x.id, mxDataRate);
      }
    }
  }

  for (Link &x : scheduled_links) {
    if (show) {
      printf("===> link %d interference %.10lf SINR %.10lf MCS %d bw %d GIVES %.3lf\n", x.id, x.interference,
             x.SINR, x.MCS, x.bw, dataRates[x.MCS][bwIdx(x.bw)]);
    }
    objective += dataRates[x.MCS][bwIdx(x.bw)];
  }
}

void Solution::insert(const Link &l) {
  scheduled_links.emplace_back(l);
  //TODO: Do I really need this?
  computeInterference();

  objective = 0.0;
  computeObjective();
}

void Solution::clearChannel(int ch) {
  set<int> MARK;
  for (Link &l : scheduled_links) {
    if (l.ch == ch) {
      MARK.insert(l.id);
    }
  }

  auto it = scheduled_links.begin();
  while (it != scheduled_links.end()) {
    if (MARK.count(it->id)) {
      it = scheduled_links.erase(it);
    } else {
      it++;
    }
  }
}

double Solution::getObjective() const {
  return objective;
}

deque<Link> Solution::getScheduledLinks() const {
  return scheduled_links;
}

deque<Link> Solution::getLinksInChannel(int ch) const {
  deque<Link> ret;
  for (const Link &l : scheduled_links) {
    if (l.ch == ch) {
      ret.emplace_back(l);
    }
  }
  return ret;
}

HeuristicDecoder::HeuristicDecoder() {
  readFile();

  chToLinks[25] = vector<int>();
  chToLinks[42] = vector<int>();
  chToLinks[43] = vector<int>();
  chToLinks[44] = vector<int>();
  chToLinks[45] = vector<int>();
  mapSplitChannels();
}

void split(Solution &dest, Solution &src, int ch) {
  int ch1 = mapChtoCh[ch].first, ch2 = mapChtoCh[ch].second;
  src.computeObjective();
  deque<Link> links = src.getLinksInChannel(ch);

  if (links.size() < 2) {
    return;
  }

  Link largest1;
  for (const Link &l : links) {
    if (l.interference > largest1.interference) {
      largest1 = l;
    }
  }

  Link largest2;
  for (const Link &l : links) {
    if (!(l == largest1) && (l.interference > largest2.interference)) {
      largest2 = l;
    }
  }

  auto it = links.begin();
  while (it != links.end()) {
    if ((*it == largest1) || (*it == largest2)) {
      it = links.erase(it);
    } else {
      it++;
    }
  }

  Solution current(src);
  current.clearChannel(ch);

//  fprintf(stderr, "apos apagar tem %lu links (array links %lu)\n", current.getScheduledLinks().size(), links.size());
  largest1.setChannel(ch1);
  largest2.setChannel(ch2);
  current.insert(largest1);
  current.insert(largest2);

  deque<int> randomPos;
  for (int i = 0; i < int(links.size()); i++) {
    randomPos.emplace_back(i);
  }

  while (!randomPos.empty()) {
    int idx = rand() % randomPos.size();

    Solution copy1(current), copy2(current);
    Link lcp1(links[idx]), lcp2(links[idx]);
    lcp1.setChannel(ch1), lcp2.setChannel(ch2);

    copy1.insert(lcp1);
    copy2.insert(lcp2);

    current = (copy1 > copy2) ? copy1 : copy2;

    swap(randomPos[idx], randomPos.back());
    randomPos.pop_back();
  }

  dest = current;
}

inline void decideBest(Solution &f, const Solution &u, const Solution &v) {
  if (u > f || v > f) {
    if (u > v) {
//      fprintf(stderr, "best solution is S1 with objective %.2lf\n", u.getObjective());
    } else if (v > u) {
//      fprintf(stderr, "DIVIDIR FO MELHOR S2 with objective %.2lf\n", v.getObjective());
    }
    f = (u > v) ? u : v;
  }
}

Solution HeuristicDecoder::generateSolution() {
  deque<int> links;
  for (int i = 0; i < nConnections; i++)
    links.emplace_back(i);

  Solution S;
  int cont = 0;
  while (!links.empty()) {
    int idx = rand() % (links.size());
    int link = links[idx];
    //fprintf(stderr, "link %d\n", link);
    //-------------
    Solution Scopy(S);
    for (auto &el : chToLinks) {
      Solution S1(S), S2;
      int ch = el.first;
      Link aux(link);
      aux.setChannel(ch);
      S1.insert(aux);
      //
      if (whichBw(ch) > 20) {
        split(S2, S1, ch);
      }

      //
      decideBest(Scopy, S1, S2);
      //
    }
    if (Scopy > S) {
      S = Scopy;
    }

    //-------------
    swap(links[idx], links.back());
    links.pop_back();
  }

  Solution seila;
  return seila; //TODO
}

void Solution::removeLink(Link link) { //TODO

}

double HeuristicDecoder::decode(std::vector<double> &chromosome) const {
  vector<pair<double, int> > auxVector;

  for (int i = 0; i < chromosome.size(); i++) {
    auxVector.push_back({chromosome[i], i});
  }

  sort(auxVector.rbegin(), auxVector.rend()); //ordenando em ordem nao crescente.

  vector<int> linksSequence;
  for (int i = 0; i < auxVector.size(); i++) {
    linksSequence.push_back(auxVector[i].second);
  }

  Solution S;
  for (int idx = 0; idx < linksSequence.size(); idx += 1) {
    int idLink = linksSequence[idx];
//    printf("tentando inserir link %d\n", idLink);

    Solution Scopy(S);
    for (auto &el : chToLinks) {
//      printf("    -> tentando canal %d\n", el.first);
      Solution S1(S), S2;

//      printf("         ===> S1 tem %d conexoes\n", S1.getScheduledLinks().size());

      int candidateChannel = el.first;

      Link auxLink(idLink);
      auxLink.setChannel(candidateChannel);
      S1.insert(auxLink);

      if (whichBw(candidateChannel) > 20) {
        split(S2, S1, candidateChannel);
      }

      decideBest(Scopy, S1, S2);
    }

    if (Scopy > S) {
      S = Scopy;
    }
  }

  double obj = S.getObjective();
//  printf("Tem %d conexoes %lf\n", S.getScheduledLinks().size(), obj);
//  this_thread::sleep_for(chrono::milliseconds(3000));
  return -obj;
};

int HeuristicDecoder::getQuantConnections() {
  return nConnections;
}

void HeuristicDecoder::setInitialTime() {
  TempoFO_StarInic = clock();
}

int Solution::getNumberOfScheduledLinks() const {
  return scheduled_links.size();
}

Link Solution::removeLinkByIndex(int index) {
  swap(scheduled_links[index], scheduled_links.back());
  Link ret = Link(scheduled_links.back());

  scheduled_links.pop_back();
  return ret;
}

void Solution::exchangeLinks(int idOldLink, Link newLink) {
  scheduled_links[idOldLink] = newLink; //TODO: is this what I meant to do?
}

vector<int> Solution::getScheduledChannels() const {
  set<int> setRet;
  for (int i = 0; i < scheduled_links.size(); i++) {
    setRet.insert(scheduled_links[i].getChannel());
  }

  return vector<int>(setRet.begin(), setRet.end());
}

double Solution::getChannelThroughput(int channel) const {
  double ret = 0.0;
  for (const Link &x : scheduled_links) {
    if (x.getChannel() == channel) {
      ret += channel;
    }
  }
  return ret;
}
