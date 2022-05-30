/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include "Genetic.h"
#include "random.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

Genetic ::Genetic(string method)
{

  int seed[4];
  int p1, p2;
  ifstream Primes("Primes");
  if (Primes.is_open())
  {
    Primes >> p1 >> p2;
  }
  else
    cerr << "PROBLEM: Unable to open Primes" << endl;
  Primes.close();
  ifstream input("seed.in");
  string property;
  if (input.is_open())
  {
    while (!input.eof())
    {
      input >> property;
      if (property == "RANDOMSEED")
      {
        input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
        rnd.SetRandom(seed, p1, p2);
      }
    }
    input.close();
  }
  else
    cerr << "PROBLEM: Unable to open seed.in" << endl;

  Set(method);
  m_Generation = 0;
  for (int i = 0; i < m_NumberCities; i++)
  {
    if (method == "circle")
    {
      double theta = rnd.Rannyu(0.0, 2.0 * M_PI);
      Positions_x.push_back(cos(theta));
      Positions_y.push_back(sin(theta));
    }
    else if (method == "square")
    {

      Positions_x.push_back(rnd.Rannyu(0.0, 1.0));
      Positions_y.push_back(rnd.Rannyu(0.0, 1.0));
    }
    else
    {
      cout << "Error" << endl;
    }
    // ADD CHECK IF cities are in the same position
    //    if(i!=0){ }
  }

  for (int i = 0; i < m_NumberCities; i++)
  {
    CityIndex.push_back(i);
  }

  ofstream City;
  City.open(m_method + "/Positions.dat");
  for (int i = 0; i < m_NumberCities; i++)
  {
    City << i << " " << Positions_x[i] << " " << Positions_y[i] << endl;
  }

  City.close();

  // Create the first population
  // shuffle the index of the cities

  for (int i = 0; i < m_NumberParents; i++)
  {
    random_shuffle(begin(CityIndex) + 1, end(CityIndex));
    Chromosomes.push_back(CityIndex);
  }
  // AbsoluteBest = Chromosomes[0];
}

Genetic ::~Genetic() {}

void Genetic::PrintChromosomes()
{

  for (int i = 0; i < m_NumberParents; i++)
  {
    for (int j = 0; j < m_NumberCities; j++)
    {
      cout << Chromosomes[i][j] << " ";
    }
    cout << endl;
  }
}

void Genetic ::Check()
{
  int count = 0;
  for (int i = 0; i < m_NumberParents; i++)
  {
    if (Chromosomes[i][0] != 0)
    {
      count++;
      cout << "error first gene is not 0" << endl;
    }
    for (int j = 0; j < m_NumberCities; j++)
    {
      for (int k = j + 1; k < m_NumberCities; k++)
      {
        if (Chromosomes[i][j] == Chromosomes[i][k])
        {
          count++;
        }
      }
    }
    if (count != 0)
    {
      cout << "Not Ok. Error in column: " << i << endl;
    }
    count = 0;
  }
}

double Genetic::FitnessCalculation(int i)
{
  double Fit = 0;
  int ind1, ind2;
  for (int j = 0; j < m_NumberCities; j++)
  {
    ind1 = (Chromosomes[i][j]);
    if (j == m_NumberCities - 1)
    {
      ind2 = Chromosomes[i][0];
    }
    else
    {
      ind2 = Chromosomes[i][j + 1];
    }
    Fit += (pow(Positions_x[ind1] - Positions_x[ind2], 2) +
            pow(Positions_y[ind1] - Positions_y[ind2], 2));
  }
  return Fit;
}

void Genetic::Print()
{
  ofstream OutputBest, OutputBestHalf;
  m_Generation++;
  OutputBest.open(m_method + "/output_best.dat", ios::app);
  OutputBest << m_Generation << " " << Fitness.back() << endl;
  OutputBest.close();

  double sum = 0;
  for (int i = 0; i < m_NumberParents / 2; i++)
  {
    sum += Fitness[(int)m_NumberParents / 2 + i];
  }
  sum = sum / ((int)m_NumberParents / 2);
  OutputBestHalf.open(m_method + "/HalfAverage.dat", ios::app);
  OutputBestHalf << m_Generation << " " << sum << endl;
  OutputBestHalf.close();
}
void Genetic::Sort()
{
  // sort the Chromosome matrix lower the fitness, higher the position
  for (int i = 0; i < m_NumberParents; i++)
  {
    for (int j = i + 1; j < m_NumberParents; j++)
    {
      if (FitnessCalculation(i) < FitnessCalculation(j))
      {
        swap(Chromosomes[i], Chromosomes[j]);
      }
    }
    Fitness.push_back(FitnessCalculation(i));

    // cout << i << " " << Fitness[i] << endl;
  }
  // Matrix problems with sort()
}

void Genetic::Selection()
{

  int r = (int)(m_NumberParents * pow(rnd.Rannyu(), m_p));
  int s;
  do
  {
    s = (int)(m_NumberParents * pow(rnd.Rannyu(), m_p));
  } while (s == r);
  Check();
  CrossOver(r, s);
}

// k is the Chromosome
void Genetic::MutationPair(int k)
{
  int ind1, ind2;
  // swap on k_th chromosome
  ind1 = (int)rnd.Rannyu(1.0, m_NumberCities);
  do
  {
    ind2 = (int)rnd.Rannyu(1.0, m_NumberCities);
  } while (ind2 == ind1);
  if (rnd.Rannyu() < m_p_Pair)
  {
    swap(Chromosomes[k][ind1], Chromosomes[k][ind2]);
  }
}

void Genetic::Set(string method)
{
  ifstream Input;
  Input.open(method + "/input.in");
  Input >> m_N_Generation;
  Input >> m_NumberCities;
  Input >> m_NumberParents;
  Input >> m_p;
  Input >> m_p_Pair;
  Input >> m_p_Inversion;
  Input >> m_p_Shift;
  Input >> m_p_Permutation;
  Input >> m_p_Selection;
  m_method = method;

  Input.close();

  cout << "Method used : " << method << endl;
  cout << "Number of Cities = " << m_NumberCities << endl;
  cout << "Number of Parents = " << m_NumberParents << endl;
  cout << "Probability used in selection = " << m_p << endl;
  cout << "Probability Pair Mutation = " << m_p_Pair << endl;
  cout << "Probability Inversion Mutation = " << m_p_Inversion << endl;
  cout << "Probability Permutation Mutation = " << m_p_Permutation << endl;
  cout << "Probability Selection = " << m_p_Selection << endl;
}

void Genetic::MutationInversion(int k)
{
  int ind1, ind2;
  // Inversion on k_th chromosome
  // ind1 is the starting position of inversion
  ind1 = (int)rnd.Rannyu(1.0, m_NumberCities - 1);
  // ind2 is the ending position of inversion
  // if ind2 ==ind1+1 reverse doesn't invert
  ind2 = (int)rnd.Rannyu(ind1 + 1, m_NumberCities);
  if (rnd.Rannyu() < m_p_Inversion)
  {
    reverse(Chromosomes[k].begin() + ind1, Chromosomes[k].begin() + ind2);
  }
}

void Genetic::MutationShift(int k)
{
  // if ind1==m_NumerCities it does nothing
  int ind1 = (int)rnd.Rannyu(1.0, m_NumberCities - 1);

  if (rnd.Rannyu() < m_p_Shift)
  {
    rotate(Chromosomes[k].begin() + 1, Chromosomes[k].end() - ind1,
           Chromosomes[k].end());
  }
}

void Genetic::MutationPermutation(int k)
{
  int ind1, ind2, ind3;
  // use a cycle as a check
  do
  {
    ind1 = (int)rnd.Rannyu(1.0,
                           m_NumberCities / 2);             // number of cities in the block
    ind2 = (int)rnd.Rannyu(1.0, m_NumberCities - 2 * ind1); // start of the
                                                            // block
    ind3 = ind1 + ind2 +
           rnd.Rannyu(0.0,
                      m_NumberCities - ind2 - 2 * ind1); // swap with this block
  } while (ind3 + ind1 > m_NumberCities);

  if (rnd.Rannyu() < m_p_Permutation)
  {
    for (int i = 0; i < ind1; i++)
    {
      swap(Chromosomes[k][ind2 + i], Chromosomes[k][ind3 + i]);
    }
  }
}

void Genetic::CrossOver(int k, int j)
{
  int ind1 = rnd.Rannyu(1.0, m_NumberCities);
  NewGeneration.push_back(Chromosomes[k]);
  double alpha = rnd.Rannyu();
  if (alpha < m_p_Selection)
  {
    sort(NewGeneration[NewGeneration.size() - 1].begin() + ind1,
         NewGeneration[NewGeneration.size() - 1].end(), [=](int a, int b)
         { return (find(Chromosomes[j].begin(), Chromosomes[j].end(), a) <
                   find(Chromosomes[j].begin(), Chromosomes[j].end(), b)); });
  }
  NewGeneration.push_back(Chromosomes[j]);

  if (alpha < m_p_Selection)
  {
    sort(NewGeneration[NewGeneration.size() - 1].begin() + ind1,
         NewGeneration[NewGeneration.size() - 1].end(), [=](int a, int b)
         { return (find(Chromosomes[k].begin(), Chromosomes[k].end(), a) <
                   find(Chromosomes[k].begin(), Chromosomes[k].end(), b)); });
  }
}

void Genetic::Simulation()
{
  Fitness.clear();
  Sort();
  BestChromosome();
  Print();
  for (int i = 0; i < m_NumberParents / 2; i++)
  {
    Selection();
  }

  Chromosomes = NewGeneration;
  NewGeneration.clear();
  Check();

  for (int k = 0; k < m_NumberParents; k++)
  {
    MutationInversion(k);
    MutationPair(k);
    MutationPermutation(k);
    MutationShift(k);
  }
  Check();
}

void Genetic::BestChromosome()
{
  ofstream Best;
  Best.open(m_method + "/BestChr.dat", ios::app);
  for (int i = 0; i < m_NumberCities; i++)
  {
    Best << Chromosomes[m_NumberParents - 1][i] << " ";
  }
  Best << endl;
  Best.close();

  /*if (AbsoluteBest > Chromosomes[m_NumberParents - 1]) {
    AbsoluteBest = Chromosomes[m_NumberParents - 1];

  }*/
}

void Genetic::CompleteSim()
{
  cout << "Start" << endl;

  string progressbar = "";

  vector<int> Test;

  for (int it = 1; it < m_N_Generation + 1; it++)
  {

    Simulation();
    if (it == m_N_Generation)
    {
      Fitness.clear();
      Sort();
      BestChromosome();
      Print();
    }

    if (it % 10 == 0)
    {

      progressbar += "=";

      cout << progressbar << " > " << (double)it * 100 / m_N_Generation << " % "
           << endl;
    }
  }
}