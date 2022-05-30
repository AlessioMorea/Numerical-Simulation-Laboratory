#include "random.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <math.h>
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
double f100(vector<double> x) {
  return exp(-2 * sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2])) / M_PI;
}
double f210(vector<double> x) {
  return x[2]*x[2]*exp(- sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2])) / (32.0*M_PI);
}



vector<double> MetropolisUniform(vector<double> x, int &Accepted, double min,
                                 double max, Random &rnd,
                                 double (&f)(vector<double>)) {
  vector<double> xproposed(3);
  vector<double> xcurrent(3);
  xcurrent = x;

  for (int i = 0; i < 3; i++) {

    xproposed[i] = xcurrent[i] + rnd.Rannyu(min, max);
  }

  double ratio = f(xproposed) / f(xcurrent);
  if (ratio >= 1) {
    xcurrent = xproposed;
    Accepted++;
  } else {
    if (rnd.Rannyu() <= ratio) {
      xcurrent = xproposed;
      Accepted++;
    }
  }

  return xcurrent;
}

vector<double> MetropolisNormal(vector<double> x, int &Accepted, double mean,double sigma,
                                Random &rnd, double (&f)(vector<double>)) {

  vector<double> xproposed(3);
  vector<double> xcurrent(3);
  xcurrent = x;

  for (int i = 0; i < 3; i++) {
    xproposed[i] = xcurrent[i] + rnd.Gauss(mean, sigma);
  }

  double ratio = f(xproposed) / f(xcurrent);

  if (ratio >= 1) {
    xcurrent = xproposed;
    Accepted++;
  } else {
    if (rnd.Rannyu() <= ratio) {
      xcurrent = xproposed;
      Accepted++;
    }
  }
  return xcurrent;
}

void Simulation(string method,int M,int N,int L,vector<double> InputVariables,Random & rnd,string folder,double inp1,double inp2,double (&f)(vector<double>),vector<double>(*function)(vector<double> , int &, double ,double ,Random &, double (&f)(vector<double>))){
  vector<vector<double>> r(M, vector<double>(3));
  r[0][0] = InputVariables[0];
  r[0][1] = InputVariables[1];
  r[0][2] = InputVariables[2];
  
  int Accepted = 0;

// Equilibration
  for (int i = 0; i < InputVariables[3]; i++) {
    r[0] = function(r[0], Accepted, inp1, inp2, rnd, f);
  }

  Accepted = 0;
  
  vector<double> distance(M);
  distance[0] = sqrt(r[0][0] * r[0][0] + r[0][1] * r[0][1] + r[0][2] * r[0][2]);
  for (int i = 0; i < M - 1; i++) {
    r[i + 1] = function(r[i], Accepted, inp1, inp2, rnd, f);
    distance[i + 1] =sqrt(r[i + 1][0] * r[i + 1][0] + r[i + 1][1] * r[i + 1][1] +r[i + 1][2] * r[i + 1][2]);
    }
  cout<< folder+"  "+method+" "<<endl;
  cout << "Percentage of accepted: " << ( double(Accepted) / double(M))*100.0 << endl;
  Data_Blocking_Method(N, L, distance, folder + method+"output.dat");  
  ofstream output(folder + "distance/"+method+"distance.dat");
  ofstream outputxyz(folder + "positions/"+method+"conf.dat");




  for (int i = 0; i < M; i++) {
    output << i << " " << distance[i] << endl;
    if(i%10==0){
    outputxyz<< i << " "<< r[i][0]<<" "<< r[i][1]<<" "<< r[i][2]<<endl;
    }
    }
    output.close();
    outputxyz.close();



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

  // start exercise 5.1
  if (argc != 2) {
    cout << "Error" << endl;
    cout << "Syntax: " << argv[0] << " <input_file_folder>" << endl;

    return -1;
  }
  string folder = argv[1];
  folder = folder + "/";
  vector<double> InputVariables(4);

  ifstream Input(folder + "input.dat");
  if (!Input.is_open()) {
    cerr << "Unable to open " << folder << "input.dat" << endl;
    return -10;
  } else {
    Input >> InputVariables[0] >> InputVariables[1] >> InputVariables[2] >>
        InputVariables[3];
  } // starting x          starting y          starting z           number of
    // equilibration steps

  int M = 1000000;    // total
  int N = 100;        // number of blocks
  int L = int(M / N); // number of throws in each block
  
  Simulation("Uniform100", M,  N,  L, InputVariables, rnd, folder, -1.223, +1.223, f100, MetropolisUniform);
  Simulation("Normal100", M, N, L, InputVariables, rnd, folder, 0, 0.758, f100, MetropolisNormal);
  Simulation("Uniform210", M, N, L, InputVariables, rnd, folder, -2.95, +2.95, f210, MetropolisUniform);
  Simulation("Normal210", M, N, L, InputVariables, rnd, folder, 0, 1.878, f210, MetropolisNormal);
//SEI QUI


  rnd.SaveSeed();

  return 0;
}
