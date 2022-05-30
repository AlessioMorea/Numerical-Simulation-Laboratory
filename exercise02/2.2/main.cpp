
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

double func(double x) { return M_PI * cos(M_PI * x / 2) / 2; }

double InvCDF(double y) {

  if (y <= 0) {
    cout << "ERROR: random variable is not positive";
    return -1;
  }
  return 1. - sqrt(1 - y);
}

double pdf(double x) { return 2. * (1 - x); }

void Print(string a, double N, vector<vector<double>> &b) {

  ofstream output(a);
  for (int i = 0; i < N; i++) {
    output << i << " " << b[i][0] << " " << b[i][1] << endl;
  }
  output.close();
}

void Data_Blocking_Method(vector<vector<double>> &results, int count, int N,
                          int L, vector<double> f) {
						
  //count is the given step
  vector<double> ave(N);
  vector<double> av2(N);
  vector<double> sum_prog(N);
  vector<double> su2_prog(N);
  vector<double> err_prog(N);

  for (int i = 0; i < N; i++) {

    ave[i] = 0;
    av2[i] = 0;
    sum_prog[i] = 0;
    su2_prog[i] = 0;
    err_prog[i] = 0;
  }
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
	//square root of the last value of sum_prog

  results[count][0] = sqrt(sum_prog[N - 1]);
	//uncertainty error on the square root calculated with propagation of errors
  results[count][1] = err_prog[N - 1] / (2. * results[count][0]);
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
  rnd.SaveSeed();

  // start exercise 2.2

  int n_steps = 101;               // number of steps per simulation
  int n_simulations = 10000;       // number of simulations
  int M = n_steps * n_simulations; // total steps

  vector<vector<double>> r(M, vector<double>(3));
  // step lenght
  double a = 1.0;

  // if p<1/3 ->x
  // if 1/3<p<2/3 -> y
  // else z

  vector<double> s(M);

  for (int j = 0; j < n_simulations; j++) {
    int k = j * n_steps;
    // coordinates
    r[k][0] = 0;
    r[k][1] = 0;
    r[k][2] = 0;
    // distance from the origin squared
    s[k] = 0;

    for (int t = 1; t < n_steps; t++) {
      int i = t + k;
      r[i][0] = r[i - 1][0];
      r[i][1] = r[i - 1][1];
      r[i][2] = r[i - 1][2];
      double t1 = rnd.Rannyu();
      double t2 = rnd.Rannyu();

      if (t1 <= 1. / 3.) {
        if (t2 <= 1. / 2.) {
          r[i][0] += a;
        } else {
          r[i][0] -= a;
        }
      } else if (t1 > 1. / 3. && t1 <= 2. / 3.) {
        if (t2 <= 1. / 2.) {
          r[i][1] += a;
        } else {
          r[i][1] -= a;
        }
      } else {

        if (t2 <= 1. / 2.) {
          r[i][2] += a;
        } else {
          r[i][2] -= a;
        }
      }
      s[i] = pow(r[i][0], 2) + pow(r[i][1], 2) + pow(r[i][2], 2);
    }
  }

  vector<double> dist(n_simulations); // the n_simulations squared distances for a given step
  vector<vector<double>> results(n_steps, vector<double>(2));//last value of the square root of the progressive average for each step and the respective uncertainty 

  for (int i = 0; i < n_steps; i++) {
    for (int j = 0; j < n_simulations; j++) {
      dist[j] = s[i + j * n_steps];
    }

    double N = 100; // number of blocks (10000 numbers in 100 blocks)
    double L = n_simulations / N; // numbers in each block
    Data_Blocking_Method(results, i, N, L, dist);
    results[0][1] = 0;
    Print("discrete.dat", n_steps, results);
  }

  vector<vector<double>> rc(M, vector<double>(3));

  vector<double> sc(M);
  for (int j = 0; j < n_simulations; j++) {
    int k = j * n_steps;
    rc[k][0] = 0;
    rc[k][1] = 0;
    rc[k][2] = 0;
    sc[k] = 0;

    for (int t = 1; t < n_steps; t++) {
      int i = t + k;
      rc[i][0] = rc[i - 1][0];
      rc[i][1] = rc[i - 1][1];
      rc[i][2] = rc[i - 1][2];

      double theta = acos(1 - 2 * rnd.Rannyu());
      double phi = rnd.Rannyu(0, M_PI * 2.0);

      rc[i][0] += a * sin(theta) * cos(phi);
      rc[i][1] += a * sin(theta) * sin(phi);
      rc[i][2] += a * cos(theta);

      sc[i] = pow(rc[i][0], 2) + pow(rc[i][1], 2) + pow(rc[i][2], 2);
    }
  }
  // the root of the average of n_simulations squared distances and its
  // uncertainty for a given step
  vector<vector<double>> resultsc(n_steps, vector<double>(2));

  for (int i = 0; i < n_steps; i++) {
    for (int j = 0; j < n_simulations; j++) {
      dist[j] = sc[i + j * n_steps];
    }

    double N = 100; // number of blocks (10000 numbers in 100 blocks)
    double L = n_simulations / N; // numbers in each block
    Data_Blocking_Method(resultsc, i, N, L, dist);
    resultsc[0][1] = 0;
    Print("continue.dat", n_steps, resultsc);
  }

  return 0;
}
