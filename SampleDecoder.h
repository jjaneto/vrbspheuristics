/*
 * SampleDecoder.h
 *
 * Any decoder must have the format below, i.e., implement the method decode(std::vector< double >&)
 * returning a double corresponding to the fitness of that vector. If parallel decoding is to be
 * used in the BRKGA framework, then the decode() method _must_ be thread-safe; the best way to
 * guarantee this is by adding 'const' to the end of decode() so that the property will be checked
 * at compile time.
 *
 * The chromosome inside the BRKGA framework can be changed if desired. To do so, just use the
 * first signature of decode() which allows for modification. Please use double values in the
 * interval [0,1) when updating, thus obeying the BRKGA guidelines.
 *
 *  Created on: Jan 14, 2011
 *      Author: rtoso
 */

#ifndef SAMPLEDECODER_H
#define SAMPLEDECODER_H

#include <list>
#include <vector>
#include <algorithm>
#include <fstream>
#include "SampleDecoder.h"
#include "Structures.h"
#include "Utility.h"

//#include "Population.h"

using namespace std;

class SampleDecoder {
private:
    clock_t TempoFO_StarInic;
    fstream dataFileRead;
    string fileName;
    vector<Coordenate> coordSender;
    vector<Coordenate> coordReceptor;

    vector<Connection> connections;


    vector<vector<Interference> > interferenceGraph;


    double **dataRates;
    double **SINR;
    int bandwidth[4], channels[24];
    int quantBandWidth, quantChannels, quantConnections, timeSlots;
    double alfa, noise, powerSender, seed;
    double penality;//O valor dessa variável é o da maior taxa de transmissão que uma conexão da rede pode assumir usando o padrão IEEE802.11ac


public:
    SampleDecoder(string _fileName);

    ~SampleDecoder();

    void cleanMemory();

    int getQuantConnections();

    double getSeed();

    void setInitialTime();

    void generateConnectionsAndInterfGraph();

    int defineMaxDataRate(double _receivedSINR, int _idBandwidth) const;

    int setBandwidthID(int _idBandwidth) const;

    double calcInterference(int _idConnection, int _time, int _spectrum, int _channel,
                            const vector<SetSlotTimes> &_timeSlotsGroups) const;

    double objectiveFunction(const vector<SetSlotTimes> &_timeSlotsGroups, std::vector<double> &chromosome) const;


    int convertChromoBandwidth(double _value) const;

    void initTimeSlotsGroups(vector<SetSlotTimes> &_timeSlotsGroups) const;

    void clearTimeSlotsGroups(vector<SetSlotTimes> &_timeSlotsGroups) const;

    void insertFreeChannels(vector<SetSlotTimes> &_timeSlotsGroups, int _bandwidth, int _idConnection,
                            int &totalSpectrumUsed, std::vector<double> &chromosome, double &totalThroughput) const;

    void insertBestChannel(vector<SetSlotTimes> &_timeSlotsGroups, int _bandwidth, int _idConnection,
                           vector<double> &connectionInterference, std::vector<double> &chromosome,
                           double &totalThroughput) const;

    double buildVRBSP(std::vector<double> &chromosome, vector<unsigned> &permutation) const;

    double decode(std::vector<double> &chromosome) const;

    void buildSolution(const std::vector<double> &chromosome, vector<unsigned> &permutation, ofstream &nameFile);

    void solutionDecode(const std::vector<double> &solutionChromosome, ofstream &nameFile);

    double solutionOjectiveFunction(const vector<SetSlotTimes> &_timeSlotsGroups);

    void writeSolution(const vector<SetSlotTimes> &_timeSlotGroups, ofstream &solucoes, double FO_Star);

    void insertFreeChannelsSolution(vector<SetSlotTimes> &_timeSlotsGroups, int _bandwidth, int _idConnection,
                                    int &totalSpectrumUsed, double &totalThroughput);

    void insertBestChannelSolution(vector<SetSlotTimes> &_timeSlotsGroups, int _bandwidth, int _idConnection,
                                   vector<double> &connectionInterference, double &totalThroughput);
};

#endif
