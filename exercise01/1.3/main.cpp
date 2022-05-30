
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

  // start exercise 1.3

  int M = 1000000;    // total throws of the needle
  int N = 100;        // number of blocks
  int L = int(M / N); // throws inside each block

  double Lenght = 0.6;
  double d = 1;

  vector<double> r(M);

  vector<vector<double>> s(M, vector<double>(2));

  vector<double> theta(M);

  for (int i = 0; i < M; i++) {

    r[i] = rnd.Rannyu(0, d); // x of the CM
    // y of the CM is not necessary

    // the coordinates that are used to determine the angle are taken inside a
    // circle of radius 1
    do {
      s[i][0] = rnd.Rannyu(-1, 1); // x
      s[i][1] = rnd.Rannyu(-1, 1); // y
    } while (hypot(s[i][1], s[i][0]) > 1);

    theta[i] = atan(s[i][1] / s[i][0]);
  }

  vector<double> ave(N);
  vector<double> av2(N);
  vector<double> sum_prog(N);
  vector<double> su2_prog(N);
  vector<double> err_prog(N);

  for (int i = 0; i < N; i++) {

    double N_hit = 0;

    for (int j = 0; j < L; j++) {
      int k = j + i * L;

      if ((r[k] - Lenght * cos(theta[k]) * 0.5 <= 0) ||
          (r[k] + Lenght * cos(theta[k]) * 0.5 >= d)) {
        N_hit++;
      }
    }

    ave[i] = 2 * Lenght * L / (N_hit * d);
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

  ofstream pi_average_output("pi_average1.dat");
  for (int i = 0; i < N; i++) {
    pi_average_output << i << " " << sum_prog[i] << " " << err_prog[i] << endl;
  }
  pi_average_output.close();

  rnd.SaveSeed();
  return 0;
}
