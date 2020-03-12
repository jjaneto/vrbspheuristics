#!/usr/bin/env bash

for i in {1..30}; do
for j in 64 128; do
./mainDE < Instancias/D250x250/U_${j}/U_${j}_${i}.txt >> outDE${j};
done;
done;