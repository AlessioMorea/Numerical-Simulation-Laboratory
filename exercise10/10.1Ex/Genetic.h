#include <vector>
#include "random.h"
#include <string>

#ifndef __Genetic__
#define __Genetic__



using namespace std;

class Genetic {
private:
  int m_rank;
  int m_size;
  int m_N_Generation;
  int m_NumberCities;  
  int m_NumberParents;
  vector<double> Positions_x;
  vector<double> Positions_y;
  vector<int> CityIndex;
  Random rnd;
  vector<vector<int>> Chromosomes;//column = single chromosome
  vector<vector<int>> NewGeneration;//column = single chromosome
  vector<double> Fitness;
  //vector<int> AbsoluteBest;
  double m_p;
  double m_p_Pair;
  double m_p_Inversion;
  double m_p_Shift;
  double m_p_Permutation;
  double m_p_Selection;
  int m_Generation;
  string m_method;
public:
  
  // constructors
  Genetic(string method,int rank,int size);

  // destructor
  ~Genetic();
  // methods
  void Set(string method,int rank,int size);
  void Check();
  void PrintChromosomes();
  double FitnessCalculation(int i);

  void Selection();
  void MutationPair(int k);
  void MutationInversion(int k);
  void MutationShift(int k);
  void MutationPermutation(int k);
  void CrossOver(int k,int j);
  void Simulation();
  void Print();
  void Sort();
  void BestChromosome();
  void CompleteSim();
  void HalfAverage();
  };

#endif // __Genetic__




/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
