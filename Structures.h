/*
 * Structures.h
 *
 *  Created on: 23 de abr de 2016
 *      Author: mauricio
 */

#ifndef STRUCTURES_H_
#define STRUCTURES_H_

#include <iostream>
#include <vector>

using namespace std;

struct Connection{
	 int time;
	 int bandwidth;
	 int frequency;
	 int id_dataRate;
	 double powerSR;//Valor da potência enviada de Si até Ri
	 double valueSINR;
	 double distanceSenderReceptor;
};

struct Interference{
	int idSender;
	int idReceptor;
	double distanceConnections;
	double valueInterference;
};

struct Coordenate{
	double x, y;
};

struct FreeChannel{
	int slotTime;
	int spectrum;
	int listChannelsPosition;
	int frequency;
	int bandwidth;
	int quantConnections;
};

struct Channel{
	int bandwidth;
	int frequency;
	double throughput;
	vector<int> listConnections;
};

struct Spectrum{
	int maxSpectrum;
	int spectrumUsed;
	vector<Channel> listChannels;
};

struct SetSlotTimes{
	int id_time;
	vector<Spectrum> listSpectrum;
};

struct TimeTableBandwidthRate{
 int timeTable[4][11];
};

#endif /* STRUCTURES_H_ */
