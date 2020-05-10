#!/usr/bin/env bash

instancias=("U_8")
mod="D250x250"
timeLimit=10

for inst in "${instancias[@]}"; do
  mkdir ${inst}
  objectivesFile="objectivesFile.txt"
  for ((i = 1; i <= 1; i++)); do
    file="${mod}/${inst}/${inst}_${i}.txt"
    solutionsFile="solutionFile_"${inst}_${i}".txt";
    echo "Trying experiment ${file}"
    ./main ../Instancias/${file} solutionsFile objectivesFile timeLimit
  done
done
