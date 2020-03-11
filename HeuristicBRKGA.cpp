//
// Created by José Joaquim on 24/02/20.
//

#include <new>
#include <float.h>
#include <string>
#include "HeuristicDecoder.h"
#include "MTRand.h"
#include "BRKGA.h"
#include "SolutionVerifier.h"

using namespace std;

int main() {
  fprintf(stderr, "BRKGA WITH HEURISTIC\n");
  const unsigned p = 100;    // size of population //***Antes era 1000 o valor desse parâmetro
  const double pe = 0.25;     // fraction of population to be the elite-set
  const double pm = 0.05;     // fraction of population to be replaced by mutants
  const double rhoe = 0.70;   // probability that offspring inherit an allele from elite parent
  const unsigned K = 1;       // number of independent populations         //***Antes era 3 o valor desse parâmetro
  const unsigned MAXT = 1;    // number of threads for parallel decoding   //***Antes era 3 o valor desse parâmetro

  const unsigned X_INTVL = 100;   // exchange best individuals at every 100 generations
  const unsigned X_NUMBER = 2;    // exchange top 2 best
  const unsigned MAX_GENS = 1000; // run for 1000 gens

  string fileName; //TODO

  HeuristicDecoder decoder;

  const unsigned n = decoder.getQuantConnections();
  double seed = time(NULL);

  MTRand rng(seed);
  BRKGA<HeuristicDecoder, MTRand> algorithm(n, p, pe, pm, rhoe, decoder, rng, K, MAXT);

  double TempoExecTotal = 0.0, TempoFO_Star = 0.0, FO_Star = DBL_MAX, FO_Min = DBL_MIN;
  int bestGeneration = 0, minGeneration = 0;
  int iterSemMelhora, iterMax = 10, quantIteracoes = 0, bestIteration = 0;

  clock_t TempoFO_StarInic;//TempoInicial

  TempoFO_StarInic = clock();
  decoder.setInitialTime();

  FO_Star = DBL_MAX;
  FO_Min = DBL_MAX * -1;
  bestGeneration = 0;
  minGeneration = 0;

  iterSemMelhora = 0;

  iterMax = 1;

  quantIteracoes = 0;
  bestIteration = 0;

  unsigned generation = 0;        // current generation
  do {
    algorithm.evolve(); // evolve the population for one generation

    if ((++generation) % X_INTVL == 0) {
      algorithm.exchangeElite(X_NUMBER); // exchange top individuals
    }

    if (algorithm.getBestFitness() < FO_Star) {
      TempoFO_Star = (((double) (clock() - TempoFO_StarInic)) / CLOCKS_PER_SEC);
      FO_Star = algorithm.getBestFitness();
      bestGeneration = generation;
      bestIteration = quantIteracoes;
    }
    quantIteracoes++;

//    if (algorithm.getBestFitness() - FO_Star < 0.0001)
//      break;

  } while ((((double) (clock() - TempoFO_StarInic)) / CLOCKS_PER_SEC) < 100);

  TempoExecTotal = (((double) (clock() - TempoFO_StarInic)) / CLOCKS_PER_SEC);

  printf("%lf\n", algorithm.getBestFitness());

  vector<double> aux = algorithm.getBestChromosome();
  return 0;
}