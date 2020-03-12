//
// Created by José Joaquim on 04/03/20.
//

#include <iostream>
#include <vector>
#include <deque>
#include <algorithm>
#include <cstdio>
#include <unordered_set>
#include <chrono>
#include <assert.h>
#include "MTRand.h"
#include "HeuristicDecoder.h"

using namespace std;

const double LOWER_BOUND = 0.0;
const double UPPER_BOUND = 1.0;

const double f__ = 0.5;
const double cr__ = 0.5;
const int numberOfDifferenceVectors = 1;

class Individual {
    double objective;
    vector<double> variables;

public:

    Individual() {
      objective = 0.0;
    }

    Individual(const int _objective, const vector<double> &_variables) : objective(_objective), variables(_variables) {}

    Individual(vector<double> &_variables) : variables(_variables) {
      Individual();
    }

    double getObjective() const {
      return objective;
    }

    void setObjective(double objective) {
      this->objective = objective;
    }

    const vector<double> &getVariables() const {
      return variables;
    }

    void setVariables(const vector<double> &variables) {
      this->variables = variables;
    }

    void setVariableValue(const int index, const double value) {
      this->variables[index] = value;
    }

    double getVariableValue(const int index) {
      return this->variables[index];
    }

    bool operator<(const Individual &o1) const {
      return objective < o1.objective;
    }
};

enum STOPPING {
    timeLimit,
    iteration
};

int populationSize;
int numberVariables;
int evaluations, maxEvaluations;
int qtdSelected;
vector<Individual> _population;
MTRand rng(time(NULL));
HeuristicDecoder decoder;
enum STOPPING criteria;
double maximumTime;
clock_t startTime;

inline void printPopulation(const vector<Individual> aux) {
  for (int i = 0; i < aux.size(); i++) {

    vector<double> variables = aux[i].getVariables();
    double obj = decoder.decode(variables);

    fprintf(stderr, "%d:", i);
    for (int i = 0; i < variables.size(); i++) {
      fprintf(stderr, " %.4lf ", variables[i]);
    }
    fprintf(stderr, "  ==> %.5lf\n", obj);
  }
}

bool isStoppingCriteriaReached() {
  if (criteria == timeLimit) {
    return (((double) (clock() - startTime)) / CLOCKS_PER_SEC) >= maximumTime;
  } else if (criteria == iteration) {
    return evaluations > maxEvaluations;
  }

  return false; //TODO
}

vector<Individual> createInitialPopulation() {
  vector<Individual> initialPop;

  for (int i = 0; i < populationSize; i++) {
    vector<double> indVariables;
    for (int i = 0; i < numberVariables; i++) {
      indVariables.emplace_back(rng.rand());
    }

    initialPop.emplace_back(Individual(indVariables));
  }

  return initialPop;
};

void evaluatePopulation(vector<Individual> &population) {
  for (int i = 0; i < populationSize; i++) {
    vector<double> variables(population[i].getVariables());

    double obj = decoder.decode(variables);
//    printf("%.4lf\n", obj);

    population[i].setObjective(obj);
  }
}

void initProgress() {
  if (criteria == timeLimit) {
    startTime = clock();
  } else if (criteria == iteration) {
    evaluations = populationSize;
  }
}

vector<Individual> selection(const vector<Individual> &solutionList) {
  unordered_set<int> toReturn;

  assert(solutionList.size() > 0);

  do {
    int randomIdx = rng.randInt(solutionList.size() - 1);

    while (toReturn.count(randomIdx))
      randomIdx = rng.randInt(solutionList.size() - 1);

    toReturn.insert(randomIdx);
  } while (toReturn.size() < qtdSelected);

  vector<Individual> inds;
  for (const int el : toReturn)
    inds.emplace_back(solutionList[el]);

  return inds;
}

double mutate(vector<vector<double> > parent, int index) {
  double value = -1000000007;


  if (numberOfDifferenceVectors == 1) {
    return parent[2][index] + f__ * (parent[0][index] - parent[1][index]);
  } else if (numberOfDifferenceVectors == 2) {
    return parent[4][index] + f__ * (parent[0][index] - parent[1][index]) + f__ * (parent[2][index] - parent[3][index]);
  }

  return value;
}

vector<Individual>
crossover(const vector<Individual> &parentSolutions, Individual child) { //TODO: Maybe there is a bug here
  vector<vector<double> > parent;

  double requiredParents = 1 + numberOfDifferenceVectors * 2;
  assert(requiredParents == qtdSelected);

  for (int i = 0; i < requiredParents; i++)
    parent.emplace_back(parentSolutions[i].getVariables());

  int jrand = rng.randInt(numberVariables - 1);

  for (int j = 0; j < numberVariables; j++) {
    if (rng.rand() || j == jrand) {
      double value = mutate(parent, j);
      child.setVariableValue(j, value);
    }
  }

//  printPopulation({Individual(child)});

  //REPAIR BOUNDS
  vector<double> childVariables = child.getVariables();
  for (int i = 0; i < childVariables.size(); i++) {
    double varValue = childVariables[i];
    if (varValue < LOWER_BOUND) {
      varValue = LOWER_BOUND;
    }

    if (varValue > UPPER_BOUND) {
      varValue = UPPER_BOUND;
    }

    childVariables[i] = varValue;
  }

  child.setVariables(childVariables);
  return {child};
}

vector<Individual> reproduction(vector<Individual> matingPopulation) {
  vector<Individual> offspringPopulation;

  for (int i = 0; i < populationSize; i++) {
    vector<Individual> parents = selection(matingPopulation);
    vector<Individual> children = crossover(parents, matingPopulation[i]);
//    printPopulation(children);

    offspringPopulation.emplace_back(children[0]);
  }

  return offspringPopulation;
}

vector<Individual> replacement(vector<Individual> &population, vector<Individual> &offspringPopulation) {
  assert(population.size() == offspringPopulation.size());
  vector<Individual> pop;

  for (int i = 0; i < populationSize; i++) {
    if (population[i] < offspringPopulation[i]) {
      pop.emplace_back(population[i]);
    } else {
      pop.emplace_back(offspringPopulation[i]);
    }
  }

  sort(pop.begin(), pop.end());
  return pop;
}

inline void updateProgress() {
  evaluations += populationSize;
}

inline void init() {
  populationSize = 100;
  numberVariables = decoder.getQuantConnections();
  criteria = timeLimit;
  maximumTime = 100;
  qtdSelected = 3;
}


int main() { //DE_RAND_1_BIN
  fprintf(stderr, "DIFFERENTIAL EVOLUTION WITH HEURISTIC\n");
  init();
  vector<Individual> offspringPopulation;

  _population = createInitialPopulation();
  evaluatePopulation(_population);
  initProgress();

  while (!isStoppingCriteriaReached()) {
    offspringPopulation = reproduction(_population);
    evaluatePopulation(offspringPopulation);
    replacement(_population, offspringPopulation);
    updateProgress();
  }

  sort(_population.begin(), _population.end());
  printf("%lf\n", _population[0].getObjective());
  return 0;
}