#! /usr/bin/env python3

import math
import sys


def bandIdx(channel):
    if channel < 25:
        return 0
    elif channel < 37:
        return 1
    elif channel < 43:
        return 2
    else:
        return 3


def computeThroughput(
    channel,
    arr,
    interferenceMatrix,
    SINR,
    dataRates,
    distanceMatrix,
    powerSender,
    alfa,
    y_set,
    interference_solution,
):
    to_ret = 0.0
    for conn1 in arr:
        interference = 0.0
        for conn2 in arr:
            if conn1 != conn2:
                interference += interferenceMatrix[conn1][conn2]

        # print(-10, conn1, interference, interference_solution[conn1])
        if interference != 0:
            if math.isclose(interference, interference_solution[conn1], rel_tol=1e-15):
                print("DIFFERENT INTERFERENCES IN conn " + str(conn1))
                print(interference, interference_solution[conn1])
                exit()

        conn_sinr = -1.0
        if interference != 0.0:
            conn_sinr = powerSender / math.pow(distanceMatrix[conn1][conn1], alfa)
            conn_sinr /= interference + noise

        # print(-20, conn1, conn_sinr)
        bw_idx = bandIdx(int(channel))

        MCS = 8 if bw_idx == 0 else 9
        if interference != 0.0:
            while MCS >= 0:
                if conn_sinr > SINR[MCS][bw_idx]:
                    break

                MCS -= 1

        if MCS < 0:
            print("ERROR")
            exit()

        conn_troughput = dataRates[MCS][bw_idx]
        # print(-30, conn1, conn_troughput)
        if [bw_idx, MCS] != y_set[conn1]:
            print(
                "erro em conn %d [%d, %d] != [%d, %d]"
                % (conn1, bw_idx, MCS, y_set[conn1][0], y_set[conn1][1])
            )
            exit()

        # print("connection %d bw %d MCS %d" % (conn1, bw_idx, MCS))
        to_ret += conn_troughput

    return to_ret


def converDBMToMW(noise):
    b = noise / 10.0
    result = math.pow(10.0, b)
    return result


def converTableToMW(SINR, SINR_):
    for i in range(len(SINR_)):
        for j in range(len(SINR_[i])):
            if SINR[i][j] != 0.0:
                b = SINR[i][j] / 10.0
                result = math.pow(10.0, b)

                SINR_[i][j] = result
            else:
                SINR_[i][j] = 0.0


def distance(a, b, c, d):
    return math.hypot((a - c), (b - d))


def distanceAndInterference(
    senders,
    receivers,
    interferenceMatrix,
    distanceMatrix,
    powerSender,
    nConnections,
    alfa,
):
    for i in range(nConnections):
        distanceMatrix.append([])
        interferenceMatrix.append([])
        X_si = receivers[i][0]
        Y_si = receivers[i][1]

        for j in range(nConnections):
            X_rj = senders[j][0]
            Y_rj = senders[j][1]

            dist = distance(X_si, Y_si, X_rj, Y_rj)
            distanceMatrix[i].append(dist)

            if i == j:
                interferenceMatrix[i].append(0.0)
            else:
                value = powerSender / math.pow(dist, alfa) if dist != 0.0 else 1e9
                interferenceMatrix[i].append(value)


def loadData(
    path, receivers, senders, dataRates, SINR, interferenceMatrix, distanceMatrix,
):
    with open(path, "r") as f:
        line = f.readline()
        line = f.readline()
        aux = line.split()
        nConnections = int(aux[0])
        ttm = float(aux[1])
        alfa = float(aux[2])
        noise = float(aux[3])
        powerSender = float(aux[4])
        nSpectrums = float(aux[5])

        if noise != 0:
            noise = converDBMToMW(noise)

        f.readline()
        for i in range(10):
            line = f.readline()
            aux = line.split()

            dataRates.append([])
            for j in range(4):
                dataRates[i].append(float(aux[j]))

        f.readline()
        for i in range(10):
            line = f.readline()
            aux = line.split()

            SINR.append([])
            for j in range(4):
                SINR[i].append(float(aux[j]))

        converTableToMW(SINR, SINR)

        f.readline()
        for i in range(nConnections):
            line = f.readline()
            aux = line.split()
            receivers.append([float(aux[0]), float(aux[1])])

        del receivers[0]

        f.readline()
        for i in range(nConnections):
            line = f.readline()
            aux = line.split()
            senders.append([float(aux[0]), float(aux[1])])

        del senders[0]

        distanceAndInterference(
            senders,
            receivers,
            interferenceMatrix,
            distanceMatrix,
            powerSender,
            nConnections,
            alfa,
        )

    return noise, powerSender, alfa, nConnections


def loadSolution(path):
    to_ret = {}
    y = {}
    interference = {}
    throughput = 0.0

    for i in range(45):
        to_ret[i] = []

    with open(path, "r") as f:
        throughput = float(f.readline())
        for line in f:
            aux = line.split(" ")
            to_ret[int(aux[1])].append(int(aux[0]))
            y[int(aux[0])] = [int(aux[2]), int(aux[3])]
            interference[int(aux[0])] = float(aux[4])

    return throughput, to_ret, y, interference


if __name__ == "__main__":
    receivers, senders, dataRates = [[]], [[]], [[]]
    SINR = []
    distanceMatrix, interferenceMatrix = [[]], [[]]

    path = sys.argv[1]
    noise, powerSender, alfa, nConnections = loadData(
        path, receivers, senders, dataRates, SINR, interferenceMatrix, distanceMatrix,
    )

    m_throughput, channels, y_set, interference = loadSolution(sys.argv[2])
    throughput = 0.0

    for k, v in channels.items():
        throughput += computeThroughput(
            k,
            v,
            interferenceMatrix,
            SINR,
            dataRates,
            distanceMatrix,
            powerSender,
            alfa,
            y_set,
            interference,
        )

    print(m_throughput, throughput)
    if m_throughput != throughput:
        print("ERROR IN THROUGHPUT")
        print(m_throughput, throughput)
        # exit()

    print(throughput)
