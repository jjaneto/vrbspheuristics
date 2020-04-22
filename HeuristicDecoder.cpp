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
  mapChtoCh[43] = {37, 38};
  mapChtoCh[37] = {25, 26};
  mapChtoCh[38] = {27, 28};
  mapChtoCh[25] = {0, 1};
  mapChtoCh[26] = {2, 3};
  mapChtoCh[27] = {4, 5};
  mapChtoCh[28] = {6, 7};
  mapChtoCh[44] = {39, 40};
  mapChtoCh[39] = {29, 30};
  mapChtoCh[40] = {31, 32};
  mapChtoCh[29] = {8, 9};
  mapChtoCh[30] = {10, 11};
  mapChtoCh[31] = {12, 13};
  mapChtoCh[32] = {14, 15};
  mapChtoCh[41] = {33, 34};
  mapChtoCh[33] = {16, 17};
  mapChtoCh[34] = {18, 19};
  mapChtoCh[42] = {35, 36};
  mapChtoCh[35] = {20, 21};
  mapChtoCh[36] = {22, 23};
}

//inline void mapSplitChannels() {
//  mapChtoCh[44] = {38, 39};
//  mapChtoCh[38] = {26, 27};
//  mapChtoCh[39] = {28, 29};
//  mapChtoCh[26] = {1, 2};
//  mapChtoCh[27] = {3, 4};
//  mapChtoCh[28] = {5, 6};
//  mapChtoCh[29] = {7, 8};
//  mapChtoCh[45] = {40, 41};
//  mapChtoCh[40] = {30, 31};
//  mapChtoCh[41] = {32, 33};
//  mapChtoCh[30] = {9, 10};
//  mapChtoCh[31] = {11, 12};
//  mapChtoCh[32] = {13, 14};
//  mapChtoCh[33] = {15, 16};
//  mapChtoCh[42] = {34, 35};
//  mapChtoCh[34] = {17, 18};
//  mapChtoCh[35] = {19, 20};
//  mapChtoCh[43] = {36, 37};
//  mapChtoCh[36] = {21, 22};
//  mapChtoCh[37] = {23, 24};
//}

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
  if (ch >= 25 && ch <= 36)
    return 40;
  else if (ch >= 37 && ch <= 42)
    return 80;
  else if (ch == 43 || ch == 44)
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
  origChannel = -1;
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
  origChannel = x.origChannel;
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
  origChannel = x.origChannel;
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

  if (o1.origChannel != o2.origChannel) {
    return false;
  }

  return true;
}

void Link::setChannel(int ch) {
  if (origChannel == -1) {
    this->origChannel = ch;
  }
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
  objectiveFlag = o1.objectiveFlag;
  for (int i = 0; i < o1.scheduled_links.size(); i++) {
    scheduled_links.push_back(o1.scheduled_links[i]);
  }
}

Solution::Solution() {
  objective = 0.0;
  objectiveFlag = true;
}

//Solution &Solution::operator=(const Solution &o1) {
//  this->objective = o1.objective;
//  this->scheduled_links = o1.scheduled_links;
//  return *this;
//}

bool operator<(const Solution &o1, const Solution &o2) {
  assert(o1.objectiveFlag && o2.objectiveFlag);
  return o1.objective < o2.objective;
}

bool operator>(const Solution &o1, const Solution &o2) {
  return operator<(o2, o1);
}

bool operator>=(const Solution &o1, const Solution &o2) {
  return !operator<(o1, o2);
}

bool operator<=(const Solution &o1, const Solution &o2) {
  return !operator>(o1, o2);
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

//      assert(u.ch - 1 >= 0 && v.ch -1 >= 0);
      if (overlap[u.ch][v.ch]) {
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

  objectiveFlag = true;
}

void Solution::insert(const Link &l) {
  scheduled_links.emplace_back(l);

  objectiveFlag = false;
  //TODO: Do I really need this?
//  computeInterference();

  objective = 0.0;
//  computeObjective();
}

void Solution::clearChannel(int ch) {
//  set<int> MARK;
//  for (Link &l : scheduled_links) {
//    if (l.ch == ch) {
//      MARK.insert(l.id);
//    }
//  }

  auto it = scheduled_links.begin();
  while (it != scheduled_links.end()) {
//    if (MARK.count(it->id)) {
    if (it->ch == ch) {
      it = scheduled_links.erase(it);
    } else {
      it++;
    }
  }
}

double Solution::getObjective(bool force) const {
  if (!force)
    assert(this->objectiveFlag);

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

  chToLinks[24] = vector<int>();
  chToLinks[41] = vector<int>();
  chToLinks[42] = vector<int>();
  chToLinks[43] = vector<int>();
  chToLinks[44] = vector<int>();
  mapSplitChannels();
  dfs(43);
  dfs(44);
  dfs(41);
  dfs(42);
  //dfs(24); //TODO: need this?
}

void Solution::setChannelOfLink(int id, int channel) {
  for (Link &x : scheduled_links) {
    if (x.id == id) {
      x.setChannel(channel);
    }
  }
  objectiveFlag = false;
}

void split(Solution &dest, Solution &src, int ch, bool ok) {
  int ch1 = mapChtoCh[ch].first, ch2 = mapChtoCh[ch].second;
  src.computeObjective();
  deque<Link> links = src.getLinksInChannel(ch);

  if (links.size() < 2) {
    if (ok && !links.empty()) {
      src.setChannelOfLink(links[0].id, ch1);
      src.computeObjective();
      dest = src;
    }
    return;
  }

//  //TODO: Finish the code below. Get the pair of links with maximum mutual interference
//  deque<Link> scheduledLinks = src.getScheduledLinks();
//
//  double mxInter = -1.0;
//  int index1 = -1, index2 = -1;
//  for (int i = 0; i < int(scheduledLinks.size()); i++) {
//    for (int j = i + 1; j < int(scheduledLinks.size()); j++) {
//      if (interferenceMatrix[scheduledLinks[j].id][scheduledLinks[i].id] > mxInter) {
//        mxInter = interferenceMatrix[scheduledLinks[j].id][scheduledLinks[i].id];
//        index1 = i;
//        index2 = j;
//      }
//    }
//  }
//
//  Link largest1(scheduledLinks[index1]);
//  Link largest2(scheduledLinks[index2]);
//
//  largest1.setChannel(ch1);
//  largest2.setChannel(ch2);
//
//  Solution newly(src);
//
//  newly.clearChannel(ch);
//
//  newly.insert(largest1);
//  newly.insert(largest2);
//
//  newly.computeObjective();

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

  Solution newly(src);

  newly.clearChannel(ch);

  largest1.setChannel(ch1);
  largest2.setChannel(ch2);
  newly.insert(largest1);
  newly.insert(largest2);

  newly.computeObjective();

  while (!links.empty()) {
    int rndIndex = rng.randInt(links.size() - 1);
    Link link = links[rndIndex];

    Link dummy1(link), dummy2(link);
    dummy1.setChannel(ch1), dummy2.setChannel(ch2);

    Solution copy1(newly), copy2(newly);

    copy1.insert(dummy1);
    copy2.insert(dummy2);

    copy1.computeObjective();
    copy2.computeObjective();

    newly = (copy1 > copy2) ? copy1 : copy2;

    swap(links[rndIndex], links.back());
    links.pop_back();
  }

  dest = newly;
}

inline void decideBest(Solution &f, const Solution &u, const Solution &v) {
  assert(u.objectiveFlag && v.objectiveFlag);
//  printf("f eh %lf comparando %lf com %lf\n", f.getObjective(), u.getObjective(), v.getObjective());
  if (u > f || v > f) {
    if (v > u) {
//      puts("alguma vez?");
    }
    f = (u > v) ? u : v;
  }
}

Solution HeuristicDecoder::generateSolution() {
  deque<int> links;
  for (int i = 0; i < nConnections; i++)
    links.emplace_back(i);

  Solution S;
  set<int> rootChannels;
  rootChannels.insert(24);
  rootChannels.insert(41);
  rootChannels.insert(42);
  rootChannels.insert(43);
  rootChannels.insert(44);
  while (!links.empty()) {
    int idx = rng.randInt(links.size() - 1);
    int link = links[idx];

    //-------------
    Solution Scopy(S);
    vector<int> availableChannels(Scopy.getScheduledChannels());
    if (availableChannels.empty()) {
//      puts("asdas");
      availableChannels = {24, 42, 41, 43, 44};
    } else {
      for (const int ch : rootChannels) {
        availableChannels.push_back(ch);
      }
    }

    for (const int ch : availableChannels) {
      Solution S1(S), S2;
      Link aux(link);
      aux.setChannel(ch);
      S1.insert(aux);
      //
      if (whichBw(ch) > 20) {
        split(S2, S1, ch);
      }
      //
      Scopy.computeObjective();
      S1.computeObjective();
      S2.computeObjective();
      decideBest(Scopy, S1, S2);
      //
    }

    if (Scopy > S) {
      S = Scopy;
      if (!rootChannels.empty()) {
        for (const int ch : S.getScheduledChannels()) {
          rootChannels.erase(ch);
        }
      }
    }

    //-------------
    swap(links[idx], links.back());
    links.pop_back();
  }

  return S;
}

bool Solution::removeLink(Link link) {
  auto it = scheduled_links.begin();
  while (it != scheduled_links.end()) {
    if (*it == link) {
      scheduled_links.erase(it);
      return true;
    } else {
      it++;
    }
  }

  return false;
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
    for (auto &el : chToLinks) { //FIXME: Change chToLinks to Solution::getScheduledLinks.
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
  Link ret(scheduled_links.back());

  scheduled_links.pop_back();
  return ret;
}

void Solution::exchangeLinks(int idOldLink, Link newLink) {
  scheduled_links[idOldLink] = newLink; //TODO: is this what I meant to do?
}

void Solution::setScheduledLinks(const deque<Link> newLinks) {
  this->scheduled_links = newLinks;
  objectiveFlag = false;
}

void Solution::addLinks(const deque<Link> &links) {
  for (const Link &link_ : links) {
    insert(link_);
  }
  objectiveFlag = false;
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

void dfs(int u, int pai, string path) {
//  printf("seeing %d with father %d and path %s\n", u, pai, path.c_str());
  if (pai != -1) {
    PATH_TO[pai][u] = path;
  }

  if (mapChtoCh.find(u) != mapChtoCh.end()) {
    dfs(mapChtoCh[u].first, u, path + "L");
    dfs(mapChtoCh[u].second, u, path + "R");
  }
}