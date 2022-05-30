#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <cmath>
#include <vector>

using namespace std;

double error(vector<double> AV, vector<double> AV2, int n)
{
  if (n == 0)
  {
    return 0;
  }
  else
  {
    return sqrt((AV2[n] - AV[n] * AV[n]) / n);
  }
}
void Print(string a, double N, vector<double> sum_prog, vector<double> err_prog)
{
  ofstream output(a);
  for (int i = 0; i < N; i++)
  {
    output << i << " " << sum_prog[i] << " " << err_prog[i] << endl;
  }
  output.close();
}

void Data_Blocking_Method(int N, int L, vector<double> f, string a)
{

  // N is the number of blocks

  // L is the number of throws in each block

  // This algorithm calculates the average of the array f and its error

  vector<double> ave(N);
  vector<double> av2(N);
  vector<double> sum_prog(N);
  vector<double> su2_prog(N);
  vector<double> err_prog(N);

  for (int i = 0; i < N; i++)
  {

    double sum = 0;

    for (int j = 0; j < L; j++)
    {

      int k = j + i * L;
      sum += f[k];
    }

    ave[i] = sum / L;
    av2[i] = (ave[i] * ave[i]);
  }

  for (int i = 0; i < N; i++)
  {

    for (int j = 0; j < i + 1; j++)
    {

      sum_prog[i] += ave[j];
      su2_prog[i] += av2[j];
    }

    sum_prog[i] = sum_prog[i] / (i + 1);
    su2_prog[i] /= (i + 1);
    err_prog[i] = error(sum_prog, su2_prog, i);
  }
  Print(a, N, sum_prog, err_prog);
}

void putAndCall() {}

int main(int argc, char *argv[])
{

  Random rnd;
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

  // start exercise 3.1

  // double t = 0;                  // time
  double T = 1;                  // delivery time
  double strikePrice = 100;      // strike price
  double riskFreeInterest = 0.1; // risk free interest
  double volatility = 0.25;      // volatility
  double startingAssetPrice = 100;
  double finalAssetPrice;

  int nData = 1000000; // total number of throws
  int nBlocks = 100;   // number of blocks

  vector<double> call;
  vector<double> put;

  for (int i = 0; i < nData; i++)
  {
    finalAssetPrice = startingAssetPrice *
                      exp((riskFreeInterest - 0.5 * pow(volatility, 2)) * T +
                          volatility * rnd.Gauss(0, T));
    call.push_back(max(finalAssetPrice - strikePrice, 0.) *
                   exp(-riskFreeInterest * T));
    put.push_back(max(strikePrice - finalAssetPrice, 0.) *
                  exp(-riskFreeInterest * T));
  }

  Data_Blocking_Method(nBlocks, int(nData / nBlocks), call, "callDiscrete.dat");
  Data_Blocking_Method(nBlocks, int(nData / nBlocks), put, "putDiscrete.dat");

  call.clear();
  put.clear();
  int nSteps = 100;

  for (int i = 0; i < nData; i++)
  {
    finalAssetPrice = startingAssetPrice;

    for (int j = 0; j < nSteps; j++)
    {

      finalAssetPrice =
          finalAssetPrice *
          exp((riskFreeInterest - 0.5 * pow(volatility, 2)) * T / nSteps +
              volatility * rnd.Gauss(0., 1.) * sqrt(T / nSteps));
    }
    call.push_back(max(finalAssetPrice - strikePrice, 0.) *
                   exp(-riskFreeInterest * T));
    put.push_back(max(strikePrice - finalAssetPrice, 0.) *
                  exp(-riskFreeInterest * T));
  }

  Data_Blocking_Method(nBlocks, nData / nBlocks, call, "callContinue.dat");
  Data_Blocking_Method(nBlocks, nData / nBlocks, put, "putContinue.dat");

  rnd.SaveSeed();

  return 0;
}
