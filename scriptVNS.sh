#!/usr/bin/env bash

instancias=("U_16")
mod="D250x250"
timeLimit=10

for inst in "${instancias[@]}"; do
  mkdir -p results/${inst}
  pathResults="results/${inst}"
  objectivesFile="objectivesFile.txt"
  for ((i = 1; i <= 1; i++)); do
    file="${mod}/${inst}/${inst}_${i}.txt"
    solutionsFile="solutionFile_"${inst}_${i}".txt";
    echo "Trying experiment ${file}"
    ./main ../Instancias/${file} ${pathResults}/"${solutionsFile}" ${pathResults}/"${objectivesFile}" ${timeLimit}
  done
done