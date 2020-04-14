/*
 * Utility.h
 *
 *  Created on: Apr 15, 2016
 *      Author: jmauricio
 */

#ifndef UTIL_UTILITY_H_
#define UTIL_UTILITY_H_

#include <iostream>
#include <vector>
#include "Structures.h"

using namespace std;

class Utility {
public:
	Utility();
	virtual ~Utility();

	double **createMatrixDouble(int _lines,int _rows);
	void printMatrix(double **matrix,int lines,int rows);
	void printMatrixSINR(double **matrix,int lines,int rows);
	void convertTableToMw(double **SINR,double **SINR_Mw,int lines,int rows);
	double convertDBMtoMW(double _value);
	void printTableMwToDBm(double **SINR_Mw,int lines,int rows);
	void printGraph(vector< vector<Interference> > _graph);
	void printConnections(vector<Connection> _connections);
	void printSetChannels(const vector<Channel> &setChannels);

	void printChannelsSlotTimeGroups(const vector<SetSlotTimes> &_timeSlotsGroups);
	void printSlotTimeGroups(const vector<SetSlotTimes> &_timeSlotsGroups);
	void printSlotTimeGroupsChannels(const vector<SetSlotTimes> &_timeSlotsGroups,int _slotTime,int _spectrum);

	void writeSlotTimeGroupsChannels(const vector<SetSlotTimes> &_timeSlotsGroups,string nameFile,double fo_Star);
	void writeHistogramTableBandwidthRate(const vector<SetSlotTimes> &_timeSlotsGroups,string nameFile,int **histogramTable,TimeTableBandwidthRate histogramTimeTable[],int _timeSlots);
};

#endif /* UTIL_UTILITY_H_ */
