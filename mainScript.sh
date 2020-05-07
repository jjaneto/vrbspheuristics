#!/usr/bin/env bash
variant=("DE_BEST_1_BIN" "DE_BEST_1_EXP" "DE_BEST_2_BIN" "DE_BEST_2_EXP" "DE_RAND_1_BIN" "DE_RAND_1_EXP" "DE_RAND_2_BIN" "DE_RAND_2_EXP")
timelimit=10

#mkdir results;
#for var in "${variant[@]}"; do
#mkdir results/${var};
#for i in 16; do
#mkdir results/${var}/U_${i};
#done;
#done;

for var in "${variant[@]}"; do
  for i in 16 32; do
    outDirectory="./results/${var}/U_${i}/"
    echo ${outDirectory}
    for ((j = 1; j <= 3; j++)); do
      ./main ${var} ${timelimit} ${outDirectory} <../Instancias/D250x250/U_${i}/U_${i}_${j}.txt 2>error${variant}.txt
    done
  done
done
