#!/usr/bin/env bash

instancias=("U_8" "U_16" "U_32" "U_64" "U_128")
mod="D250x250"
timeLimit=600
directoryResults="results"

for ((i = 1; i <= 30; i++)); do
  objectivesFile="objectivesFile.txt"
  for inst in "${instancias[@]}"; do
    pathResults="${directoryResults}/${inst}"

    mkdir -p ${pathResults}

    file="${mod}/${inst}/${inst}_${i}.txt"
    solutionsFile="solutionFile_"${inst}_${i}".txt"
    args="../Instancias/${file} ${pathResults}/"${solutionsFile}" ${pathResults}/"${objectivesFile}" ${timeLimit}"
    echo "Trying experiment ${file} with arguments ${args}"
    ./main ${args}
  done
done
