/*
 * SampleDecoder.cpp
 *
 *  Created on: Jan 14, 2011
 *      Author: rtoso
 */

#include <iostream>
#include <math.h>
#include <vector>
#include <list>
#include <algorithm>
#include <fstream>
#include "SampleDecoder.h"
#include "Structures.h"
#include "Utility.h"
#include <limits.h>
#include <float.h>
#include <iomanip>

using namespace std;

SampleDecoder::SampleDecoder(string _fileName) {

  fileName = _fileName;

  Coordenate auxCoordenate;
  Utility utilities;

  int readerInt, count;
  double x_Receptor, y_Receptor, x_Sender, y_Sender;
  double readerDouble;

  int quantSpectrums, sizeSpectrum;
  quantBandWidth = 4;
  quantChannels = 24;

  //***Aqui!!! cout<<"Nome do Arquivo: "<<fileName.c_str()<<endl;

  dataFileRead.open(fileName.c_str(), fstream::in);


  if (dataFileRead.is_open()) {
    dataFileRead.seekg(0, ios::beg);
    //Aqui!!! cout<<"Operation successfully performed Aqui!!!"<<endl;
  } else {
    cout << "Não foi possivel abrir o arquivo!" << endl;
    dataFileRead.clear();
    exit(1);
  }

  dataFileRead >> seed;
  dataFileRead >> quantConnections;
  dataFileRead >> timeSlots;
  dataFileRead >> alfa;
  dataFileRead >> noise;
  dataFileRead >> powerSender;

  dataFileRead >> quantSpectrums;
  for (int c = 0; c < quantSpectrums; c++) {
    dataFileRead >> sizeSpectrum;
  }
  penality = -780.0 * 10.0 * quantConnections;

  if (noise != 0) {
    noise = utilities.convertDBMtoMW(noise);
  }

  dataRates = utilities.createMatrixDouble(10, 4);

  for (int i = 0; i < 10; i++) {
    for (int j = 0; j < 4; j++) {
      if (dataFileRead.good()) {
        dataFileRead >> readerDouble;
        dataRates[i][j] = readerDouble * 10.0;//Antes era dataRates[i][j]= readerDouble;

        dataRates[i][j] = floor(dataRates[i][j]);
      } else {
        exit(1);
      }
    }
  }
  SINR = utilities.createMatrixDouble(10, 4);

  //***Aqui!!! cout<<"SINR"<<endl;
  for (int i = 0; i < 10; i++) {
    for (int j = 0; j < 4; j++) {
      if (dataFileRead.good()) {
        dataFileRead >> readerDouble;
        SINR[i][j] = readerDouble;
      } else {
        exit(1);
      }
    }
  }

  utilities.convertTableToMw(SINR, SINR, 10, 4);

  count = 0;
  while (dataFileRead.good() && (count < quantConnections)) {

    dataFileRead >> x_Receptor;
    dataFileRead >> y_Receptor;

    //Decidiu-se por não usar mais métodos de arrendondamento nessa parte para evitar arredondamentos diferentes dos valores testados
    //no processo de geração das instâncias. Portanto, optou-se por apenas ler os valores das coordenadas salvos
    auxCoordenate.x = x_Receptor;

    auxCoordenate.y = y_Receptor;

    coordReceptor.push_back(auxCoordenate);

    count++;
  }

  count = 0;
  while (dataFileRead.good() && (count < quantConnections)) {

    dataFileRead >> x_Sender;
    dataFileRead >> y_Sender;

    //Decidiu-se por não usar mais métodos de arrendondamento nessa parte para evitar arredondamentos diferentes dos valores testados
    //no processo de geração das instâncias. Portanto, optou-se por apenas ler os valores das coordenadas salvos
    auxCoordenate.x = x_Sender;

    auxCoordenate.y = y_Sender;


    coordSender.push_back(auxCoordenate);
    count++;
  }

  generateConnectionsAndInterfGraph();

  dataFileRead.clear();
  dataFileRead.close();
}

SampleDecoder::~SampleDecoder() {

}

void SampleDecoder::cleanMemory() {
  coordSender.clear();
  coordReceptor.clear();
  connections.clear();
  interferenceGraph.clear();
  delete dataRates;
  delete SINR;
  dataRates = NULL;
  SINR = NULL;
}

int SampleDecoder::getQuantConnections() {
  return quantConnections;
}

double SampleDecoder::getSeed() {
  return seed;
}

void SampleDecoder::setInitialTime() {
  this->TempoFO_StarInic = clock();
}

void SampleDecoder::generateConnectionsAndInterfGraph() {
  Connection auxConnection;
  Interference auxInterference;
  double coordX1, coordY1, coordX2, coordY2, xResult, yResult, result, interferenSenderReceptor;

  //***cout<<"\n ====== Generate Connections and Interference Graph  ====== \n"<<endl;

  for (int i = 0; i < quantConnections; i++) {

    //***cout<<"Connection= "<<i<<endl;

    coordX1 = coordReceptor[i].x;
    coordY1 = coordReceptor[i].y;

    coordX2 = coordSender[i].x;
    coordY2 = coordSender[i].y;

    //Calculation of the distance between Receptor[i] and Sender[i] of the Connection[i]
    xResult = (coordX1 - coordX2) * (coordX1 - coordX2);
    yResult = (coordY1 - coordY2) * (coordY1 - coordY2);

    result = sqrt(xResult + yResult);

    auxConnection.distanceSenderReceptor = result;
    auxConnection.powerSR = powerSender / pow(result,
                                              alfa);//Aqui é calculada a força do sinal transmitido entre o transmissor i e o receptor i.

    connections.push_back(auxConnection);

    //Inicia a construção do grafo de interferências
    interferenceGraph.push_back(vector<Interference>());

    for (int j = 0; j < quantConnections; j++) {

      if (i != j) {

        auxInterference.idReceptor = i;
        auxInterference.idSender = j;
        coordX2 = coordSender[j].x;
        coordY2 = coordSender[j].y;

        if (coordX1 == coordX2 && coordY1 == coordY2) {

          result = 0.0;
          //distance between Receptor[i] and Sender[j]
          auxInterference.distanceConnections = result;
          //Interference between Receptor[i] and Sender[j] ( P_jj/(distance_ij)^alfa)
          //OBS: See restriction in mathematical model
          interferenSenderReceptor = 1000000000;
          auxInterference.valueInterference = interferenSenderReceptor;

          interferenceGraph[interferenceGraph.size() - 1].push_back(auxInterference);

        } else {
          //Calculation of the distance between Receptor[i] and Sender[j] of the Connections i and j, respectively.
          xResult = (coordX1 - coordX2) * (coordX1 - coordX2);
          yResult = (coordY1 - coordY2) * (coordY1 - coordY2);

          result = sqrt(xResult + yResult);

          //distance between Receptor[i] and Sender[j]
          auxInterference.distanceConnections = result;
          //OBS: See restriction in mathematical model
          interferenSenderReceptor = powerSender / pow(result, alfa);
          auxInterference.valueInterference = interferenSenderReceptor;
          interferenceGraph[interferenceGraph.size() - 1].push_back(auxInterference);
        }

      }
    }
  }

  coordReceptor.clear();
  coordSender.clear();
}

int SampleDecoder::defineMaxDataRate(double _receivedSINR, int _idBandwidth) const {
  int idDataRate = -1;
  for (int i = 9; i >= 0; i--) {
    //Essa condição serve para evitar que algoritmo use uma taxa de transmissão que não existe na tabela referente ao padrão IEEE 802.11ac.
    if (_idBandwidth == 0 && i == 9) {
      i = 9 - 1;
    }
    if (_receivedSINR >= SINR[i][_idBandwidth]) {
      idDataRate = i;
      break;
    }
  }

  return idDataRate;
}

int SampleDecoder::setBandwidthID(int _idBandwidth) const {

  int result = 0;

  if (_idBandwidth == 20) {
    result = 0;//20MHz a largura de banda
  } else if (_idBandwidth == 40) {
    result = 1;//40MHz a largura de banda
  } else if (_idBandwidth == 80) {
    result = 2;//80MHz a largura de banda
  } else if (_idBandwidth == 160) {
    result = 3;//160MHz a largura de banda
  }

  return result;
}

double SampleDecoder::calcInterference(int _idConnection, int _time, int _spectrum, int _channel,
                                       const vector<SetSlotTimes> &_timeSlotsGroups) const {

  double interference = 0.0;
  int conflictConnection, posConnectionInterference = 0;
  //Se duas conexões transmitem ao mesmo tempo t, usam a mesma frequência f e largura de banda b
  //significa que há interferência entre as conexões _idConnection e i
  for (unsigned i = 0;
       i < _timeSlotsGroups[_time].listSpectrum[_spectrum].listChannels[_channel].listConnections.size(); i++) {

    conflictConnection = _timeSlotsGroups[_time].listSpectrum[_spectrum].listChannels[_channel].listConnections[i];

    if (_idConnection != conflictConnection) {


      if (_idConnection < conflictConnection) {
        posConnectionInterference = conflictConnection - 1;
      } else if (_idConnection > conflictConnection) {
        posConnectionInterference = conflictConnection;
      }

      interference = interference + interferenceGraph[_idConnection][posConnectionInterference].valueInterference;

    }
  }
  return interference;
}

double
SampleDecoder::objectiveFunction(const vector<SetSlotTimes> &_timeSlotsGroups, std::vector<double> &chromosome) const {
  double dataRate, interference = 0.0;
  double fo = 0.0;
  int id_DataRate, id_Bandwidth;

  //Aqui!!! cout<<"====== Objective Function with Penality ======"<<endl;

  for (unsigned i = 0; i < _timeSlotsGroups.size(); i++) {

    for (unsigned j = 0; j < _timeSlotsGroups[i].listSpectrum.size(); j++) {

      for (unsigned s = 0; s < _timeSlotsGroups[i].listSpectrum[j].listChannels.size(); s++) {

        for (unsigned l = 0; l < _timeSlotsGroups[i].listSpectrum[j].listChannels[s].listConnections.size(); l++) {

          interference = calcInterference(_timeSlotsGroups[i].listSpectrum[j].listChannels[s].listConnections[l], i, j,
                                          s, _timeSlotsGroups);

          if (interference + noise > 0) {
            dataRate = connections[_timeSlotsGroups[i].listSpectrum[j].listChannels[s].listConnections[l]].powerSR /
                       (interference + noise);
          } else if (interference + noise == 0) {
            dataRate = SINR[9][3];
          }

          id_Bandwidth = setBandwidthID(_timeSlotsGroups[i].listSpectrum[j].listChannels[s].bandwidth);
          id_DataRate = defineMaxDataRate(dataRate, id_Bandwidth);

          //Nessa estrutura IF-ELSE é trocado o valor do gene responsável por definir a largura de banda que o canal em que uma conexão pode estar
          // pelo valor da largura de banda do canal em que a conexão realmente está localizada.
          if (_timeSlotsGroups[i].listSpectrum[j].listChannels[s].bandwidth == 20 && (convertChromoBandwidth(
                  chromosome[(_timeSlotsGroups[i].listSpectrum[j].listChannels[s].listConnections[l] * 2) + 1]) !=
                                                                                      20)) {

            chromosome[(_timeSlotsGroups[i].listSpectrum[j].listChannels[s].listConnections[l] * 2) + 1] =
                    (0 + 0.25) / 2.0;

            //Aqui!!!cout<<"Bandwidth= "<<_timeSlotsGroups[i].listSpectrum[j].listChannels[s].bandwidth<<" Chromosome["<<_timeSlotsGroups[i].listSpectrum[j].listChannels[s].listConnections[l]+1<<"]= "<<chromosome[ _timeSlotsGroups[i].listSpectrum[j].listChannels[s].listConnections[l]+1 ]<<endl;

          } else if (_timeSlotsGroups[i].listSpectrum[j].listChannels[s].bandwidth == 40 && (convertChromoBandwidth(
                  chromosome[(_timeSlotsGroups[i].listSpectrum[j].listChannels[s].listConnections[l] * 2) + 1]) !=
                                                                                             40)) {

            chromosome[(_timeSlotsGroups[i].listSpectrum[j].listChannels[s].listConnections[l] * 2) + 1] =
                    (0.5 + 0.25) / 2.0;

            //Aqui!!!cout<<"Bandwidth= "<<_timeSlotsGroups[i].listSpectrum[j].listChannels[s].bandwidth<<" Chromosome["<<_timeSlotsGroups[i].listSpectrum[j].listChannels[s].listConnections[l]+1<<"]= "<<chromosome[ _timeSlotsGroups[i].listSpectrum[j].listChannels[s].listConnections[l]+1 ]<<endl;

          } else if (_timeSlotsGroups[i].listSpectrum[j].listChannels[s].bandwidth == 80 && (convertChromoBandwidth(
                  chromosome[(_timeSlotsGroups[i].listSpectrum[j].listChannels[s].listConnections[l] * 2) + 1]) !=
                                                                                             80)) {

            chromosome[(_timeSlotsGroups[i].listSpectrum[j].listChannels[s].listConnections[l] * 2) + 1] =
                    (0.75 + 0.5) / 2.0;

            //Aqui!!!cout<<"Bandwidth= "<<_timeSlotsGroups[i].listSpectrum[j].listChannels[s].bandwidth<<" Chromosome["<<_timeSlotsGroups[i].listSpectrum[j].listChannels[s].listConnections[l]+1<<"]= "<<chromosome[ _timeSlotsGroups[i].listSpectrum[j].listChannels[s].listConnections[l]+1 ]<<endl;

          } else if (_timeSlotsGroups[i].listSpectrum[j].listChannels[s].bandwidth == 160 && (convertChromoBandwidth(
                  chromosome[(_timeSlotsGroups[i].listSpectrum[j].listChannels[s].listConnections[l] * 2) + 1]) !=
                                                                                              160)) {

            chromosome[(_timeSlotsGroups[i].listSpectrum[j].listChannels[s].listConnections[l] * 2) + 1] =
                    (1.0 + 0.75) / 2.0;

          }
          fo = fo + dataRates[id_DataRate][id_Bandwidth];
        }

      }

    }

  }

  //Aqui!!! cout<<"FO= "<<fo<<endl;
  //Aqui!!! getchar();


  return (-1.0) * fo;
}

int SampleDecoder::convertChromoBandwidth(double _value) const {
  int result = 0;
  if (_value < 0.25) {
    result = 20;//20MHz a largura de banda
  } else if (_value < 0.5) {
    result = 40;//40MHz a largura de banda
  } else if (_value < 0.75) {
    result = 80;//80MHz a largura de banda
  } else if (_value < 1.0) {
    result = 160;//160MHz a largura de banda
  }

  return result;
}


void SampleDecoder::initTimeSlotsGroups(vector<SetSlotTimes> &_timeSlotsGroups) const {

  for (unsigned i = 0; i < this->timeSlots; i++) {

    _timeSlotsGroups.push_back(SetSlotTimes());

    //_timeSlotsGroups[i].totalSpectrumUsed= 0;

    for (unsigned j = 0; j < 3; j++) {
      _timeSlotsGroups[i].listSpectrum.push_back(Spectrum());
    }

    _timeSlotsGroups[i].listSpectrum[0].spectrumUsed = 0;
    _timeSlotsGroups[i].listSpectrum[1].spectrumUsed = 0;
    _timeSlotsGroups[i].listSpectrum[2].spectrumUsed = 0;

    _timeSlotsGroups[i].listSpectrum[0].maxSpectrum = 160;
    _timeSlotsGroups[i].listSpectrum[1].maxSpectrum = 240;
    _timeSlotsGroups[i].listSpectrum[2].maxSpectrum = 100;

  }
}

void SampleDecoder::clearTimeSlotsGroups(vector<SetSlotTimes> &_timeSlotsGroups) const {
  for (int i = 0; i < this->timeSlots; i++) {
    _timeSlotsGroups[i].listSpectrum.clear();
  }
  _timeSlotsGroups.clear();
}


//insertFreeChannels(timeSlotsGroups, convertChromoBandwidth(chromosome[permutation[i] + 1]),
//        positionChromosome, timeSlotsGroupsFULL, chromosome, totalThroughput);
//Método responsável por inserir as conexões de forma organizada no conjunto timeSlotsGroups
void SampleDecoder::insertFreeChannels(vector<SetSlotTimes> &_timeSlotsGroups, int _bandwidth, int _idConnection,
                                       int &totalSpectrumUsed, std::vector<double> &chromosome,
                                       double &totalThroughput) const {
  Utility util;
  Channel auxChannel, auxChannelCase2;
  int auxfreeSpectrum, auxfreeSpectrumCase2, auxTime, auxSpectrum, auxChannelListPosition;
  bool channelCreated = false;

  double dataRate;
  int id_DataRate, id_Bandwidth;
  unsigned positionNewChannel;

  auxChannelCase2.frequency = -1;
  auxfreeSpectrumCase2 = 0;
  auxfreeSpectrum = INT_MIN;

  for (unsigned i = 0; i < _timeSlotsGroups.size(); i++) {
    for (unsigned j = 0; j < _timeSlotsGroups[i].listSpectrum.size(); j++) {

      //Verifica se em tal espectro há algum canal livre se há um espaço (mesmo que não seja CONTÍGUO) do tamanho do canal solicitado pela conexão
      if (_timeSlotsGroups[i].listSpectrum[j].spectrumUsed < _timeSlotsGroups[i].listSpectrum[j].maxSpectrum) {

        auxfreeSpectrum =
                _timeSlotsGroups[i].listSpectrum[j].maxSpectrum - _timeSlotsGroups[i].listSpectrum[j].spectrumUsed;

        if (_bandwidth <= auxfreeSpectrum) {
          if (_timeSlotsGroups[i].listSpectrum[j].spectrumUsed == 0) {
            auxChannel.frequency = 0;
          } else {
            int szAux = int(_timeSlotsGroups[i].listSpectrum[j].listChannels.size() - 1);
            auxChannel.frequency =
                    _timeSlotsGroups[i].listSpectrum[j].listChannels[szAux].frequency +
                    _timeSlotsGroups[i].listSpectrum[j].listChannels[szAux].bandwidth;
          }

          dataRate = SINR[9][3];
          id_Bandwidth = setBandwidthID(_bandwidth);
          id_DataRate = defineMaxDataRate(dataRate, id_Bandwidth);
          totalThroughput = totalThroughput + dataRates[id_DataRate][id_Bandwidth];
          _timeSlotsGroups[i].listSpectrum[j].listChannels.push_back(auxChannel);
          positionNewChannel = _timeSlotsGroups[i].listSpectrum[j].listChannels.size() - 1;

          _timeSlotsGroups[i].listSpectrum[j].listChannels[positionNewChannel].bandwidth = _bandwidth;
          _timeSlotsGroups[i].listSpectrum[j].listChannels[positionNewChannel].throughput = dataRates[id_DataRate][id_Bandwidth];

          _timeSlotsGroups[i].listSpectrum[j].listChannels[positionNewChannel].listConnections.push_back(_idConnection);

          _timeSlotsGroups[i].listSpectrum[j].spectrumUsed =
                  _timeSlotsGroups[i].listSpectrum[j].spectrumUsed + _bandwidth;

          totalSpectrumUsed = totalSpectrumUsed + _bandwidth;

          //É atribuído o valor -1 novamente para evitar que possíveis canais encontrados,
          //de tamanho menor, sejam criados desnecessariamente no CASO 2
          auxChannelCase2.frequency = -1;
          j = _timeSlotsGroups[i].listSpectrum.size();
          channelCreated = true;
          break;
        } else if (auxfreeSpectrum > auxfreeSpectrumCase2) {
          if (_timeSlotsGroups[i].listSpectrum[j].spectrumUsed == 0) {
            auxChannelCase2.frequency = 0;
          } else {
            int szAux = int(_timeSlotsGroups[i].listSpectrum[j].listChannels.size() - 1);
            auxChannelCase2.frequency =
                    _timeSlotsGroups[i].listSpectrum[j].listChannels[szAux].frequency +
                    _timeSlotsGroups[i].listSpectrum[j].listChannels[szAux].bandwidth;
          }

          auxfreeSpectrumCase2 = auxfreeSpectrum;
          //IF-ELSE usado para identificar a largura de banda do canal em tal espaço livre
          if (auxfreeSpectrumCase2 >= 80) {
            auxChannelCase2.bandwidth = 80;
          } else if (auxfreeSpectrumCase2 >= 40) {
            auxChannelCase2.bandwidth = 40;
          } else {
            auxChannelCase2.bandwidth = 20;
          }

          auxTime = i;
          auxSpectrum = j;
        }
      }
    }

    if (channelCreated == true) {
      i = _timeSlotsGroups.size();
    }
  }

  if (auxChannelCase2.frequency >= 0) {
    //Update throughput vector
    dataRate = SINR[9][3];

    id_Bandwidth = setBandwidthID(auxChannelCase2.bandwidth);
    id_DataRate = defineMaxDataRate(dataRate, id_Bandwidth);
    totalThroughput = totalThroughput + dataRates[id_DataRate][id_Bandwidth];

    _timeSlotsGroups[auxTime].listSpectrum[auxSpectrum].listChannels.push_back(auxChannelCase2);
    positionNewChannel = _timeSlotsGroups[auxTime].listSpectrum[auxSpectrum].listChannels.size() - 1;
    _timeSlotsGroups[auxTime].listSpectrum[auxSpectrum].listChannels[positionNewChannel].throughput = dataRates[id_DataRate][id_Bandwidth];
    _timeSlotsGroups[auxTime].listSpectrum[auxSpectrum].listChannels[positionNewChannel].listConnections.push_back(
            _idConnection);
    _timeSlotsGroups[auxTime].listSpectrum[auxSpectrum].spectrumUsed =
            _timeSlotsGroups[auxTime].listSpectrum[auxSpectrum].spectrumUsed + auxChannelCase2.bandwidth;
    totalSpectrumUsed = totalSpectrumUsed + auxChannelCase2.bandwidth;

    //Nessa estrutura IF-ELSE é trocado o valor do gene responsável por definir a largura de banda que o canal em que uma conexão pode estar
    // pelo valor da largura de banda do canal em que a conexão realmente está localizada.
    if (auxChannelCase2.bandwidth == 20) {
      chromosome[(_idConnection * 2) + 1] = (0 + 0.25) / 2.0;
    } else if (auxChannelCase2.bandwidth == 40) {
      chromosome[(_idConnection * 2) + 1] = (0.5 + 0.25) / 2.0;
    } else if (auxChannelCase2.bandwidth == 80) {
      chromosome[(_idConnection * 2) + 1] = (0.75 + 0.5) / 2.0;
    } else if (auxChannelCase2.bandwidth == 160) {
      chromosome[(_idConnection * 2) + 1] = (1.0 + 0.75) / 2.0;
    }
  }
}


void SampleDecoder::insertBestChannel(vector<SetSlotTimes> &_timeSlotsGroups, int _bandwidth, int _idConnection,
                                      vector<double> &connectionInterference, std::vector<double> &chromosome,
                                      double &totalThroughput) const {

  FreeChannel viableChannel, auxChannel;
  double interference = 0.0, sumSINR = 0.0, valueSINR = 0.0, bestViableSINR, bestSINR;
  int posInterference = -1, id_Bandwidth, id_DataRate;
  bool viableTransmissions, channelTransmission, viableTransmission;

  Utility util;
  viableTransmission = false;

  bestViableSINR = DBL_MIN;
  bestSINR = DBL_MIN;

  double interferenceID_Connection, interferenceChannelInsertion, channelThroughput, auxThroughput, bestChannelThroughput = DBL_MIN;
  double newTotalThroughput = totalThroughput;
  unsigned sizeChannel;
  int numberConnection;

  for (unsigned i = 0; i < _timeSlotsGroups.size(); i++) {
    for (unsigned j = 0; j < _timeSlotsGroups[i].listSpectrum.size(); j++) {
      for (unsigned s = 0; s < _timeSlotsGroups[i].listSpectrum[j].listChannels.size(); s++) {
        channelTransmission = true;

        channelThroughput = 0.0;
        interferenceID_Connection = 0.0;

        for (unsigned l = 0; l < _timeSlotsGroups[i].listSpectrum[j].listChannels[s].listConnections.size(); l++) {
          interference = connectionInterference[_timeSlotsGroups[i].listSpectrum[j].listChannels[s].listConnections[l]];

          if (_timeSlotsGroups[i].listSpectrum[j].listChannels[s].listConnections[l] < _idConnection) {
            posInterference = _idConnection - 1;
          } else if (_timeSlotsGroups[i].listSpectrum[j].listChannels[s].listConnections[l] > _idConnection) {
            posInterference = _idConnection;
          }
          interference = interference +
                         interferenceGraph[_timeSlotsGroups[i].listSpectrum[j].listChannels[s].listConnections[l]][posInterference].valueInterference;
          if (interference + noise > 0) {
            valueSINR = connections[_timeSlotsGroups[i].listSpectrum[j].listChannels[s].listConnections[l]].powerSR /
                        (interference + noise);
          } else if (interference + noise == 0) {
            valueSINR = SINR[9][3];
          }


          id_Bandwidth = setBandwidthID(_timeSlotsGroups[i].listSpectrum[j].listChannels[s].bandwidth);
          id_DataRate = defineMaxDataRate(valueSINR, id_Bandwidth);

          if (id_DataRate >= 0) {

            if (_idConnection < _timeSlotsGroups[i].listSpectrum[j].listChannels[s].listConnections[l]) {
              posInterference = _timeSlotsGroups[i].listSpectrum[j].listChannels[s].listConnections[l] - 1;
            } else if (_idConnection > _timeSlotsGroups[i].listSpectrum[j].listChannels[s].listConnections[l]) {
              posInterference = _timeSlotsGroups[i].listSpectrum[j].listChannels[s].listConnections[l];
            }

            channelThroughput = channelThroughput + dataRates[id_DataRate][id_Bandwidth];
            interferenceID_Connection =
                    interferenceID_Connection + interferenceGraph[_idConnection][posInterference].valueInterference;

          } else {

            channelTransmission = false;
            l = _timeSlotsGroups[i].listSpectrum[j].listChannels[s].listConnections.size();
            break;
          }
        }

        if (channelTransmission == true) {
          if (interferenceID_Connection + noise > 0) {
            valueSINR = connections[_idConnection].powerSR / (interferenceID_Connection + noise);
          } else if (interferenceID_Connection + noise == 0) {
            valueSINR = SINR[9][3];
          }

          id_Bandwidth = setBandwidthID(_timeSlotsGroups[i].listSpectrum[j].listChannels[s].bandwidth);
          id_DataRate = defineMaxDataRate(valueSINR, id_Bandwidth);

          if (id_DataRate >= 0) {
            channelThroughput = channelThroughput + dataRates[id_DataRate][id_Bandwidth];
            auxThroughput = totalThroughput - _timeSlotsGroups[i].listSpectrum[j].listChannels[s].throughput;
            auxThroughput = auxThroughput + channelThroughput;
            if (auxThroughput > newTotalThroughput) {

              newTotalThroughput = auxThroughput;

              bestChannelThroughput = channelThroughput;
              interferenceChannelInsertion = interferenceID_Connection;
              viableChannel.slotTime = i;
              viableChannel.spectrum = j;
              viableChannel.listChannelsPosition = s;
              viableTransmission = true;
            }
          }
        }
      }
    }
  }

  unsigned chromosomeBandwith;

  if (viableTransmission == true) {
    sizeChannel = _timeSlotsGroups[viableChannel.slotTime].listSpectrum[viableChannel.spectrum]
            .listChannels[viableChannel.listChannelsPosition].listConnections.size();

    for (unsigned i = 0; i < sizeChannel; i++) {
      numberConnection = _timeSlotsGroups[viableChannel.slotTime].listSpectrum[viableChannel.spectrum]
              .listChannels[viableChannel.listChannelsPosition].listConnections[i];

      if (numberConnection < _idConnection) {
        posInterference = _idConnection - 1;
      } else if (numberConnection > _idConnection) {
        posInterference = _idConnection;
      }

      connectionInterference[numberConnection] = connectionInterference[numberConnection] +
                                                 interferenceGraph[numberConnection][posInterference].valueInterference;
    }

    connectionInterference[_idConnection] = interferenceChannelInsertion;

    totalThroughput = totalThroughput -
                      _timeSlotsGroups[viableChannel.slotTime].listSpectrum[viableChannel.spectrum]
                              .listChannels[viableChannel.listChannelsPosition].throughput;

    totalThroughput = totalThroughput + bestChannelThroughput;

    _timeSlotsGroups[viableChannel.slotTime].listSpectrum[viableChannel.spectrum]
            .listChannels[viableChannel.listChannelsPosition].throughput = bestChannelThroughput;

    _timeSlotsGroups[viableChannel.slotTime].listSpectrum[viableChannel.spectrum]
            .listChannels[viableChannel.listChannelsPosition].listConnections.push_back(_idConnection);

    chromosomeBandwith = _timeSlotsGroups[viableChannel.slotTime].listSpectrum[viableChannel.spectrum]
            .listChannels[viableChannel.listChannelsPosition].bandwidth;

    //Nessa estrutura IF-ELSE é trocado o valor do gene responsável por definir a largura de banda que o canal
    // em que uma conexão pode estar
    // pelo valor da largura de banda do canal em que a conexão realmente está localizada.
    if (chromosomeBandwith != _bandwidth) {
      if (chromosomeBandwith == 20) {
        chromosome[(_idConnection * 2) + 1] = (0 + 0.25) / 2.0;
      } else if (chromosomeBandwith == 40) {
        chromosome[(_idConnection * 2) + 1] = (0.5 + 0.25) / 2.0;
      } else if (chromosomeBandwith == 80) {
        chromosome[(_idConnection * 2) + 1] = (0.75 + 0.5) / 2.0;
      } else if (chromosomeBandwith == 160) {
        chromosome[(_idConnection * 2) + 1] = (1.0 + 0.75) / 2.0;
      }
    }
  }
}


double SampleDecoder::buildVRBSP(std::vector<double> &chromosome, vector<unsigned> &permutation) const {
  vector<SetSlotTimes> timeSlotsGroups;

  Utility util;
  double fitness = 0.0;
  int totalSpectrum = (160 + 240 + 100) * timeSlots, timeSlotsGroupsFULL, positionChromosome;

  initTimeSlotsGroups(timeSlotsGroups);
  timeSlotsGroupsFULL = 0;

  std::vector<double> connectionInterference;
  double totalThroughput = 0.0;

  connectionInterference.assign(quantConnections, 0.0);//inicializa com zeros

  unsigned i = 0;
  while ((i < permutation.size()) && (timeSlotsGroupsFULL < totalSpectrum)) {
    positionChromosome = permutation[i] / 2;//Inserção pela ordem definida pelo permutation
//    printf("(1) tentando conn %d, band %d (%lf)\n", positionChromosome, convertChromoBandwidth(chromosome[permutation[i] + 1]), chromosome[permutation[i] + 1]);
    insertFreeChannels(timeSlotsGroups, convertChromoBandwidth(chromosome[permutation[i] + 1]),
                       positionChromosome, timeSlotsGroupsFULL, chromosome, totalThroughput);
    i++;
  }

  while (i < permutation.size()) {
    positionChromosome = permutation[i] / 2;//Inserção pela ordem definida pelo permutation
//    printf("(2) tentando conn %d, band %d (%lf)\n", positionChromosome, convertChromoBandwidth(chromosome[permutation[i] + 1]), chromosome[permutation[i] + 1]);
    insertBestChannel(timeSlotsGroups, convertChromoBandwidth(chromosome[permutation[i] + 1]),
                      positionChromosome, connectionInterference, chromosome, totalThroughput);
    i++;
  }

  fitness = totalThroughput * (-1.0);
  clearTimeSlotsGroups(timeSlotsGroups);
  connectionInterference.clear();

  return fitness;
}

// Runs in \Theta(n \log n): //Exemplo consultado em brkga API
double SampleDecoder::decode(std::vector<double> &chromosome) const {
  // Here we store the fitness of the chromosome
  double myFitness = 0.0;

  std::vector<std::pair<double, unsigned> > ranking;//ranking(chromosome.size());

//  if ((((double) (clock() - TempoFO_StarInic)) / CLOCKS_PER_SEC) > 600) {
//    myFitness = quantConnections * 780.0 * 10.0;
//    return myFitness;
//  }

  // Here we compute the fitness as the inner product of the random-key vector,
  // i.e. the chromosome, and the vector of indices [0,1,2,...,n-1]
  for (unsigned i = 0; i < chromosome.size(); i = i + 2) {
    ranking.push_back(std::pair<double, unsigned>(chromosome[i], i));
  }

  // Here we sort ranking, which will then produce a permutation of [n]
  // in pair::second:
  sort(ranking.begin(), ranking.end());
  // Here we obtain the permutation of [n]:

  vector<unsigned> permutation;
  for (vector<pair<double, unsigned> >::const_iterator i = ranking.begin(); i != ranking.end(); i++) {
    permutation.push_back(i->second);
  }

  // Here we return the fitness of chromosome:
  myFitness = buildVRBSP(chromosome, permutation);
  return myFitness;
}

double SampleDecoder::solutionOjectiveFunction(const vector<SetSlotTimes> &_timeSlotsGroups) {
  double dataRate, interference = 0.0;
  double fo = 0.0;
  int id_DataRate, id_Bandwidth;

  //Aqui!!! cout<<"====== Objective Function with Penality ======"<<endl;

  for (unsigned i = 0; i < _timeSlotsGroups.size(); i++) {

    for (unsigned j = 0; j < _timeSlotsGroups[i].listSpectrum.size(); j++) {

      for (unsigned s = 0; s < _timeSlotsGroups[i].listSpectrum[j].listChannels.size(); s++) {

        for (unsigned l = 0; l < _timeSlotsGroups[i].listSpectrum[j].listChannels[s].listConnections.size(); l++) {

          interference = calcInterference(_timeSlotsGroups[i].listSpectrum[j].listChannels[s].listConnections[l], i, j,
                                          s, _timeSlotsGroups);

          if (interference + noise > 0) {
            dataRate = connections[_timeSlotsGroups[i].listSpectrum[j].listChannels[s].listConnections[l]].powerSR /
                       (interference + noise);
          } else if (interference + noise == 0) {
            dataRate = SINR[9][3];
          }

          id_Bandwidth = setBandwidthID(_timeSlotsGroups[i].listSpectrum[j].listChannels[s].bandwidth);
          id_DataRate = defineMaxDataRate(dataRate, id_Bandwidth);

          fo = fo + dataRates[id_DataRate][id_Bandwidth];
        }
      }
    }
  }

  return (-1.0) * fo;
}


void
SampleDecoder::buildSolution(const std::vector<double> &chromosome, vector<unsigned> &permutation, ofstream &nameFile) {
  vector<SetSlotTimes> timeSlotsGroups;

  Utility util;
  double fitness = 0.0;
  int totalSpectrum = (160 + 240 + 100) * timeSlots, timeSlotsGroupsFULL, positionChromosome;

  initTimeSlotsGroups(timeSlotsGroups);
  timeSlotsGroupsFULL = 0;

  std::vector<double> connectionInterference;
  double totalThroughput = 0.0;

  connectionInterference.assign(quantConnections, 0.0);//inicializa com zeros

  unsigned i = 0;
  while ((i < permutation.size()) && (timeSlotsGroupsFULL < totalSpectrum)) {
    positionChromosome = permutation[i] / 2;//Inserção pela ordem definida pelo permutation
    insertFreeChannelsSolution(timeSlotsGroups, convertChromoBandwidth(chromosome[permutation[i] + 1]),
                               positionChromosome, timeSlotsGroupsFULL, totalThroughput);
    i++;
  }

  while (i < permutation.size()) {
    positionChromosome = permutation[i] / 2;//Inserção pela ordem definida pelo permutation
    insertBestChannelSolution(timeSlotsGroups, convertChromoBandwidth(chromosome[permutation[i] + 1]),
                              positionChromosome, connectionInterference, totalThroughput);
    i++;
  }

  fitness = solutionOjectiveFunction(timeSlotsGroups);
  fitness = totalThroughput;
  writeSolution(timeSlotsGroups, nameFile, fitness);


  clearTimeSlotsGroups(timeSlotsGroups);
  connectionInterference.clear();
}


void SampleDecoder::solutionDecode(const std::vector<double> &solutionChromosome, ofstream &nameFile) {

  std::vector<std::pair<double, unsigned> > ranking;//ranking(chromosome.size());

  // Here we compute the fitness as the inner product of the random-key vector,
  // i.e. the chromosome, and the vector of indices [0,1,2,...,n-1]
  for (unsigned i = 0; i < solutionChromosome.size(); i = i + 2) {
    ranking.push_back(std::pair<double, unsigned>(solutionChromosome[i], i));
    //myFitness += (double(i+1) * chromosome[i]);
  }

  sort(ranking.begin(), ranking.end());
  // Here we obtain the permutation of [n]:

  vector<unsigned> permutation;
  for (vector<pair<double, unsigned> >::const_iterator i = ranking.begin(); i != ranking.end(); i++) {
    permutation.push_back(i->second);
  }

  // Here we decode the chromosome:
  buildSolution(solutionChromosome, permutation, nameFile);

}


void SampleDecoder::writeSolution(const vector<SetSlotTimes> &_timeSlotGroups, ofstream &solucoes, double FO_Star) {

  int quantCanais = 0, contador;
  Utility util;

  //***Aqui!!! util.printChannelsSlotTimeGroups(_timeSlotGroups);

  solucoes << (FO_Star * (-1)) / 10.0 << "\n";
  solucoes << _timeSlotGroups.size() << "\n";
  solucoes << "\n";

  for (int i = 0; i < _timeSlotGroups.size(); i++) {

    contador = 0;

    quantCanais = _timeSlotGroups[i].listSpectrum[0].listChannels.size() +
                  _timeSlotGroups[i].listSpectrum[1].listChannels.size() +
                  _timeSlotGroups[i].listSpectrum[2].listChannels.size();

    solucoes << quantCanais << "\n";//Aqui é a linha que indica a quantidade de canais

    solucoes << i << "\n";//Aqui é o índice do intervalo de tempo

    for (int j = 0; j <
                    3; j++)//Aqui é 3 pois o espectro extra que está no primeiro intervalo de tempo é para armazenas as conexões que não transmitiram
    {

      if (_timeSlotGroups[i].listSpectrum[j].listChannels.size() > 0) {


        for (int l = 0; l < _timeSlotGroups[i].listSpectrum[j].listChannels.size(); l++) {

          //Aqui é o índice do canal, largura de banda, throughput e quantidade de conexões dele
          //solucoes<<contador<<" "<<_timeSlotGroups[i].listSpectrum[j].listChannels[l].bandwidth<<" "<<_timeSlotGroups[i].listSpectrum[j].listChannels[l].throughput<<" "<<_timeSlotGroups[i].listSpectrum[j].listChannels[l].listConnections.size()<<"\n";
          solucoes << contador << " " << _timeSlotGroups[i].listSpectrum[j].listChannels[l].bandwidth << " " << 0.0
                   << " " << _timeSlotGroups[i].listSpectrum[j].listChannels[l].listConnections.size() << "\n";


          solucoes << "   ";
          for (int c = 0; c < _timeSlotGroups[i].listSpectrum[j].listChannels[l].listConnections.size(); c++) {
            solucoes << _timeSlotGroups[i].listSpectrum[j].listChannels[l].listConnections[c] << "  ";
          }

          solucoes << "\n\n";

          contador++;

        }

      }

    }


    solucoes << "\n\n\n";
  }


}


void
SampleDecoder::insertFreeChannelsSolution(vector<SetSlotTimes> &_timeSlotsGroups, int _bandwidth, int _idConnection,
                                          int &totalSpectrumUsed, double &totalThroughput) {
  Utility util;
  Channel auxChannel, auxChannelCase2;
  int auxfreeSpectrum, auxfreeSpectrumCase2, auxTime, auxSpectrum, auxChannelListPosition;
  bool channelCreated = false;


  double dataRate;
  int id_DataRate, id_Bandwidth;
  unsigned positionNewChannel;

  auxChannelCase2.frequency = -1;
  auxfreeSpectrumCase2 = 0;
  auxfreeSpectrum = INT_MIN;

  for (unsigned i = 0; i < _timeSlotsGroups.size(); i++) {

    for (unsigned j = 0; j < _timeSlotsGroups[i].listSpectrum.size(); j++) {

      //Verifica se em tal espectro há algum canal livre se há um espaço (mesmo que não seja CONTÍGUO) do tamanho do canal solicitado pela conexão
      if (_timeSlotsGroups[i].listSpectrum[j].spectrumUsed < _timeSlotsGroups[i].listSpectrum[j].maxSpectrum) {


        auxfreeSpectrum =
                _timeSlotsGroups[i].listSpectrum[j].maxSpectrum - _timeSlotsGroups[i].listSpectrum[j].spectrumUsed;

        if (_bandwidth <= auxfreeSpectrum) {

          if (_timeSlotsGroups[i].listSpectrum[j].spectrumUsed == 0) {
            auxChannel.frequency = 0;
          } else {
            auxChannel.frequency = _timeSlotsGroups[i].listSpectrum[j].listChannels[
                                           _timeSlotsGroups[i].listSpectrum[j].listChannels.size() - 1].frequency +
                                   _timeSlotsGroups[i].listSpectrum[j].listChannels[
                                           _timeSlotsGroups[i].listSpectrum[j].listChannels.size() - 1].bandwidth;
          }
          //Update throughput vector
          dataRate = SINR[9][3];

          id_Bandwidth = setBandwidthID(_bandwidth);
          id_DataRate = defineMaxDataRate(dataRate, id_Bandwidth);

          totalThroughput = totalThroughput + dataRates[id_DataRate][id_Bandwidth];


          _timeSlotsGroups[i].listSpectrum[j].listChannels.push_back(auxChannel);

          positionNewChannel = _timeSlotsGroups[i].listSpectrum[j].listChannels.size() - 1;


          _timeSlotsGroups[i].listSpectrum[j].listChannels[positionNewChannel].bandwidth = _bandwidth;
          _timeSlotsGroups[i].listSpectrum[j].listChannels[positionNewChannel].throughput = dataRates[id_DataRate][id_Bandwidth];

          _timeSlotsGroups[i].listSpectrum[j].listChannels[positionNewChannel].listConnections.push_back(_idConnection);

          _timeSlotsGroups[i].listSpectrum[j].spectrumUsed =
                  _timeSlotsGroups[i].listSpectrum[j].spectrumUsed + _bandwidth;

          totalSpectrumUsed = totalSpectrumUsed + _bandwidth;

          //É atribuído o valor -1 novamente para evitar que possíveis canais encontrados, de tamanho menor, sejam criados desnecessariamente no CASO 2
          auxChannelCase2.frequency = -1;
          j = _timeSlotsGroups[i].listSpectrum.size();
          channelCreated = true;
          break;


        }//if(_bandwidth > auxfreeSpectrum)
        else if (auxfreeSpectrum > auxfreeSpectrumCase2) {


          if (_timeSlotsGroups[i].listSpectrum[j].spectrumUsed == 0) {
            //***Aqui!!!  cout<<"\n---> === Primeiro Canal no Espectro "<<j<<" ==="<<endl;
            auxChannelCase2.frequency = 0;
          } else {
            auxChannelCase2.frequency = (_timeSlotsGroups[i].listSpectrum[j].listChannels[
                                                 _timeSlotsGroups[i].listSpectrum[j].listChannels.size() -
                                                 1].frequency + _timeSlotsGroups[i].listSpectrum[j].listChannels[
                                                 _timeSlotsGroups[i].listSpectrum[j].listChannels.size() -
                                                 1].bandwidth);
          }
          auxfreeSpectrumCase2 = auxfreeSpectrum;

          if (auxfreeSpectrumCase2 >= 80) {
            auxChannelCase2.bandwidth = 80;
          } else if (auxfreeSpectrumCase2 >= 40) {
            auxChannelCase2.bandwidth = 40;
          } else {
            auxChannelCase2.bandwidth = 20;
          }

          auxTime = i;
          auxSpectrum = j;
        }


      }
    }

    if (channelCreated == true) {
      i = _timeSlotsGroups.size();
    }

  }


  if (auxChannelCase2.frequency >= 0) {

    //Update throughput vector
    dataRate = SINR[9][3];

    id_Bandwidth = setBandwidthID(auxChannelCase2.bandwidth);
    id_DataRate = defineMaxDataRate(dataRate, id_Bandwidth);

    totalThroughput = totalThroughput + dataRates[id_DataRate][id_Bandwidth];


    _timeSlotsGroups[auxTime].listSpectrum[auxSpectrum].listChannels.push_back(auxChannelCase2);

    positionNewChannel = _timeSlotsGroups[auxTime].listSpectrum[auxSpectrum].listChannels.size() - 1;


    _timeSlotsGroups[auxTime].listSpectrum[auxSpectrum].listChannels[positionNewChannel].throughput = dataRates[id_DataRate][id_Bandwidth];

    _timeSlotsGroups[auxTime].listSpectrum[auxSpectrum].listChannels[positionNewChannel].listConnections.push_back(
            _idConnection);

    _timeSlotsGroups[auxTime].listSpectrum[auxSpectrum].spectrumUsed =
            _timeSlotsGroups[auxTime].listSpectrum[auxSpectrum].spectrumUsed + auxChannelCase2.bandwidth;
    totalSpectrumUsed = totalSpectrumUsed + auxChannelCase2.bandwidth;
  }


}


void SampleDecoder::insertBestChannelSolution(vector<SetSlotTimes> &_timeSlotsGroups, int _bandwidth, int _idConnection,
                                              vector<double> &connectionInterference, double &totalThroughput) {

  FreeChannel viableChannel, auxChannel;
  double interference = 0.0, sumSINR = 0.0, valueSINR = 0.0, bestViableSINR, bestSINR;
  int posInterference = -1, id_Bandwidth, id_DataRate;
  bool viableTransmissions, channelTransmission, viableTransmission;


  Utility util;

  viableTransmission = false;

  bestViableSINR = DBL_MIN;
  bestSINR = DBL_MIN;


  double interferenceID_Connection, interferenceChannelInsertion, channelThroughput, auxThroughput, bestChannelThroughput = DBL_MIN;
  double newTotalThroughput = totalThroughput;
  unsigned sizeChannel;
  int numberConnection;


  for (unsigned i = 0; i < _timeSlotsGroups.size(); i++) {

    for (unsigned j = 0; j < _timeSlotsGroups[i].listSpectrum.size(); j++) {
      for (unsigned s = 0; s < _timeSlotsGroups[i].listSpectrum[j].listChannels.size(); s++) {

        channelTransmission = true;


        channelThroughput = 0.0;
        interferenceID_Connection = 0.0;


        for (unsigned l = 0; l < _timeSlotsGroups[i].listSpectrum[j].listChannels[s].listConnections.size(); l++) {

          interference = connectionInterference[_timeSlotsGroups[i].listSpectrum[j].listChannels[s].listConnections[l]];

          if (_timeSlotsGroups[i].listSpectrum[j].listChannels[s].listConnections[l] < _idConnection) {
            posInterference = _idConnection - 1;
          } else if (_timeSlotsGroups[i].listSpectrum[j].listChannels[s].listConnections[l] > _idConnection) {
            posInterference = _idConnection;
          }

          interference = interference +
                         interferenceGraph[_timeSlotsGroups[i].listSpectrum[j].listChannels[s].listConnections[l]][posInterference].valueInterference;

          if (interference + noise > 0) {
            valueSINR = connections[_timeSlotsGroups[i].listSpectrum[j].listChannels[s].listConnections[l]].powerSR /
                        (interference + noise);
          } else if (interference + noise == 0) {
            valueSINR = SINR[9][3];
          }


          id_Bandwidth = setBandwidthID(_timeSlotsGroups[i].listSpectrum[j].listChannels[s].bandwidth);
          id_DataRate = defineMaxDataRate(valueSINR, id_Bandwidth);

          if (id_DataRate >= 0) {

            if (_idConnection < _timeSlotsGroups[i].listSpectrum[j].listChannels[s].listConnections[l]) {
              posInterference = _timeSlotsGroups[i].listSpectrum[j].listChannels[s].listConnections[l] - 1;
            } else if (_idConnection > _timeSlotsGroups[i].listSpectrum[j].listChannels[s].listConnections[l]) {
              posInterference = _timeSlotsGroups[i].listSpectrum[j].listChannels[s].listConnections[l];
            }
            channelThroughput = channelThroughput + dataRates[id_DataRate][id_Bandwidth];

            interferenceID_Connection =
                    interferenceID_Connection + interferenceGraph[_idConnection][posInterference].valueInterference;
          } else {

            channelTransmission = false;
            l = _timeSlotsGroups[i].listSpectrum[j].listChannels[s].listConnections.size();
            break;
          }
        }


        if (channelTransmission == true) {


          if (interferenceID_Connection + noise > 0) {
            valueSINR = connections[_idConnection].powerSR / (interferenceID_Connection + noise);
          } else if (interferenceID_Connection + noise == 0) {
            valueSINR = SINR[9][3];
          }


          id_Bandwidth = setBandwidthID(_timeSlotsGroups[i].listSpectrum[j].listChannels[s].bandwidth);
          id_DataRate = defineMaxDataRate(valueSINR, id_Bandwidth);

          if (id_DataRate >= 0) {
            channelThroughput = channelThroughput + dataRates[id_DataRate][id_Bandwidth];
            auxThroughput = totalThroughput - _timeSlotsGroups[i].listSpectrum[j].listChannels[s].throughput;

            auxThroughput = auxThroughput + channelThroughput;
            if (auxThroughput > newTotalThroughput) {

              newTotalThroughput = auxThroughput;

              bestChannelThroughput = channelThroughput;
              interferenceChannelInsertion = interferenceID_Connection;
              viableChannel.slotTime = i;
              viableChannel.spectrum = j;
              viableChannel.listChannelsPosition = s;

              viableTransmission = true;
            }

          }
        }
      }

    }

  }

  if (viableTransmission == true) {

    //***Aqui!!! cout<<"\n ======= Insertion Viable ======="<<endl;



    sizeChannel = _timeSlotsGroups[viableChannel.slotTime].listSpectrum[viableChannel.spectrum].listChannels[viableChannel.listChannelsPosition].listConnections.size();
    for (unsigned i = 0; i < sizeChannel; i++) {

      numberConnection = _timeSlotsGroups[viableChannel.slotTime].listSpectrum[viableChannel.spectrum].listChannels[viableChannel.listChannelsPosition].listConnections[i];

      if (numberConnection < _idConnection) {
        posInterference = _idConnection - 1;
      } else if (numberConnection > _idConnection) {
        posInterference = _idConnection;
      }


      connectionInterference[numberConnection] = connectionInterference[numberConnection] +
                                                 interferenceGraph[numberConnection][posInterference].valueInterference;
    }


    connectionInterference[_idConnection] = interferenceChannelInsertion;

    totalThroughput = totalThroughput -
                      _timeSlotsGroups[viableChannel.slotTime].listSpectrum[viableChannel.spectrum].listChannels[viableChannel.listChannelsPosition].throughput;

    totalThroughput = totalThroughput + bestChannelThroughput;
    _timeSlotsGroups[viableChannel.slotTime].listSpectrum[viableChannel.spectrum].listChannels[viableChannel.listChannelsPosition].throughput = bestChannelThroughput;

    _timeSlotsGroups[viableChannel.slotTime].listSpectrum[viableChannel.spectrum].listChannels[viableChannel.listChannelsPosition].listConnections.push_back(
            _idConnection);
  }
}
