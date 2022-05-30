#include "random.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;
double error(vector<double> AV, vector<double> AV2, int n) {

  if (n == 0) {
    return 0;
  } else {
    return sqrt((AV2[n] - AV[n] * AV[n]) / n);
  }
}

void Print(string a, double N, vector<double> sum_prog,
           vector<double> err_prog) {
  ofstream output(a);
  for (int i = 0; i < N; i++) {
    output << i << " " << sum_prog[i] << " " << err_prog[i] << endl;
  }

  output.close();
}

void Data_Blocking_Method(int N, int L, vector<double> f, string a) {
  // N is the number of blocks
  // L is the number of throws in each block
  // This algorithm calculates the average of the array f and its error
  vector<double> ave(N);
  vector<double> av2(N);
  vector<double> sum_prog(N);
  vector<double> su2_prog(N);
  vector<double> err_prog(N);

  for (int i = 0; i < N; i++) {
    double sum = 0;

    for (int j = 0; j < L; j++) {

      int k = j + i * L;
      sum += f[k];
    }

    ave[i] = sum / L;
    av2[i] = (ave[i] * ave[i]);
  }

  for (int i = 0; i < N; i++) {

    for (int j = 0; j < i + 1; j++) {

      sum_prog[i] += ave[j];
      su2_prog[i] += av2[j];
    }

    sum_prog[i] = sum_prog[i] / (i + 1);
    su2_prog[i] /= (i + 1);
    err_prog[i] = error(sum_prog, su2_prog, i);
  }
  Print(a, N, sum_prog, err_prog);
}

int main(int argc, char *argv[]) {

  Random rnd;
  int seed[4];
  int p1, p2;
  ifstream Primes("Primes");
  if (Primes.is_open()) {
    Primes >> p1 >> p2;
  } else
    cerr << "PROBLEM: Unable to open Primes" << endl;
  Primes.close();

  ifstream input("seed.in");
  string property;
  if (input.is_open()) {
    while (!input.eof()) {
      input >> property;
      if (property == "RANDOMSEED") {
        input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
        rnd.SetRandom(seed, p1, p2);
      }
    }
    input.close();
  } else
    cerr << "PROBLEM: Unable to open seed.in" << endl;
  

  // start exercises 1.1.1 and 1.1.2

  int M = 1000000;    // total throws
  int N = 100;        // number of blocks
  int L = int(M / N); // throws inside each block

  vector<double> r(M);
  vector<double> r_variance(M);

  for (int i = 0; i < M; i++) {

    r[i] = rnd.Rannyu();
    r_variance[i] = pow(r[i] - 0.5, 2);
  }

  Data_Blocking_Method(N, L, r, "average.dat");
  Data_Blocking_Method(N, L, r_variance, "variance.dat");

  // 1.1.3

  int n_bins = 100;
  int n_throws = 10000;

  double expected_value = n_throws / n_bins;

  vector<double> chi_squared(n_bins);
  vector<int> bin_count(n_bins); // observed number of throws in each bin

  double sum = 0;
  for (int i = 0; i < n_bins; i++) {
    sum = 0;
    for (int j = 0; j < n_bins; j++) {
      bin_count[j] = 0;
    }

    for (int j = 0; j < n_throws; j++) {
      bin_count[int(floor(r[j + n_throws * i] * n_bins))]++;
    }

    for (int j = 0; j < n_bins; j++) {
      sum += (bin_count[j] - expected_value) * (bin_count[j] - expected_value);
    }
    chi_squared[i] = sum / expected_value;
  }

  ofstream chi_output("chi_squared.dat");
  for (int i = 0; i < n_bins; i++) {
    chi_output << i << " " << chi_squared[i] << endl;
  }

  chi_output.close();

  rnd.SaveSeed();
  return 0;
}



