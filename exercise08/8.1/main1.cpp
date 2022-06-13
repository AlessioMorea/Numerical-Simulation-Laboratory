/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include "MD_MC.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <string>

using namespace std;

double fT(double x) {
  return exp(-(x - mu) * (x - mu) / (2.0 * sigma * sigma)) +
         exp(-(x + mu) * (x + mu) / (2.0 * sigma * sigma));
}

double Potential(double x) { return pow(x, 4) - 5.0 * pow(x, 2) / 2.0; }

double SecondDerivative(double x) {

  double dx = exp(-(x - mu) * (x - mu) / (2.0 * sigma * sigma));
  double sx = exp(-(x + mu) * (x + mu) / (2.0 * sigma * sigma));
  double D2 =
      +pow(sigma, -2) * ((pow(x - mu, 2) / pow(sigma, 2)) * dx +
                         (pow(x + mu, 2) / pow(sigma, 2)) * sx - dx - sx);

  return D2;
}

int main() {
  folder =   "1method/";
  double sigma_best, sigma_old;
  double mu_best, mu_old;
  double energy_old, energy = 15.0, energy_best;
  double delta_sigma = 0.1, delta_mu = 0.1;
  double alpha=0.995;
  // int attempted_mu = 0, accepted_mu = 0;
  // int attempted_sigma = 0, accepted_sigma = 0;
  double P;
  ofstream data,best,data_old;


  Input(); // Inizialization

  // Equilibration
  /*for (int i = 0; i < n_equilibration; i++) {
    Move();
  }*/
  sigma_old = sigma;
  mu_old = mu;
  energy_old = energy;
  energy_best = energy_old;
  sigma_best = sigma_old;
  temp=1.0;
  int count=0;
  err_Ham=0;


data_old.open(folder + "output_OldData.dat", ios::app);
data_old <<count<<" "<<temp << " " << mu_old << " " << sigma_old <<" "<<energy_old<<" "<<err_Ham<< endl;
data_old.close();




  while(temp>=0.0001) {
    temp *= alpha;
    count++;
    beta = 1.0 / temp;
    mu = abs(mu_old + rnd.Gauss(0.0, delta_mu));
    sigma = abs(sigma_old + rnd.Gauss(0.0, delta_sigma));
    if (sigma <= 0.01) {
      sigma = sigma + 0.1;
    }

    for (int iblk = 1; iblk <= nblk; iblk++) // Simulation
    {
      Reset(iblk); // Reset block averages
      for (int istep = 1; istep <= nstep; istep++) {
        Move();
        Measure();
        Accumulate(); // Update block averages
      }
      Averages(iblk); // Print results for current block
    }
    energy = glob_av[iv] / (double)nblk;

    data.open(folder + "output_PossibleE.dat", ios::app);
    
    data <<count<<" "<<temp << " " << mu << " " << sigma <<" "<<energy<<" "<<err_Ham<< endl;
    


    data.close();

    P = exp(-beta * (energy - energy_old));
    if (P >= rnd.Rannyu()) {
      // Update
      energy_old = energy;
      sigma_old = sigma;
      mu_old =mu;
      // accepted_sigma++;
      
      data_old.open(folder + "output_OldData.dat", ios::app);
      data_old <<count<<" "<<temp << " " << mu_old << " " << sigma_old <<" "<<energy_old<<" "<<err_Ham<< endl;
      data_old.close();

      if (energy_old < energy_best) {
        sigma_best = sigma;
        energy_best = energy_old;
        mu_best = mu;
      }
    }
    // attempted_sigma++;
  }

  // cout << "Ac Mu rate:" << (double)accepted_mu / attempted_mu << endl;
  // cout << "Ac sigma rate:" << (double)accepted_sigma / attempted_sigma <<
  // endl;
  temp=0.0;
  sigma = sigma_best;
  mu = mu_best;
  nstep=nstep*10;
  ofstream Gofr1;
  Gofr1.open(folder + "output_gofr1.dat", ios::app);
  for (int iblk = 1; iblk <= nblk; iblk++) // Simulation
  {
    Reset(iblk); // Reset block averages
    for (int istep = 1; istep <= nstep; istep++) {
      Move();
      Measure();
      Accumulate(); // Update block averages
      Gofr1<<x<<endl;
    }
    Averages(iblk); // Print results for current block
  }
  Gofr1.close();
cout << "sigma old = " << sigma_old << " | mu old = " << mu_old<< " | energy = " << energy_old << endl;
cout << "sigma best = " << sigma_best << " | mu best = " << mu_best<< " | energy best = " << energy_best << endl;

best.open(folder+"best.dat");
best << "sigma old = " << sigma_old << " | mu old = " << mu_old<< " | energy = " << energy_old << endl;
best << "sigma best = " << sigma_best << " | mu best = " << mu_best<< " | energy best = " << energy_best << endl;







  //cout<<glob_av[iv]/(double)nblk<<endl;;
  return 0;
}

void Input(void) {
  ifstream ReadInput, ReadConf, ReadVelocity, Primes, Seed;

  // Read seed for random numbers
  int p1, p2;
  Primes.open("Primes");
  Primes >> p1 >> p2;
  Primes.close();

  // Read input informations
  ReadInput.open(folder + "input.in");
  ReadInput >> restart;

  if (restart)
    { 
      Seed.open("seed.out");}
  else{
    
    Seed.open("seed.in");
    }

  Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  
  
  rnd.SetRandom(seed, p1, p2);
  Seed.close();
  x = 0.0;

  ReadInput >> delta;

  ReadInput >> nblk;

  ReadInput >> nstep;
  ReadInput >> n_equilibration;
  ReadInput >> mu;
  ReadInput >> sigma;

  cout << "The program perform Metropolis moves with uniform translations"
       << endl;
  cout << "Moves parameter = " << delta << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of equilibration steps = " << n_equilibration << endl;
  cout << "Number of steps in one block = " << nstep << endl;
  cout << "Initial mu = " << mu << endl;
  cout << "Initial sigma = " << sigma << endl << endl;
  cout << "Initial temp = " << temp << endl;

  ReadInput.close();

  // Prepare arrays for measurements
  iv = 0; // Hamiltonian/psi
 
  n_props = 1; // Number of observables
  //igofr = 1;
  //nbins = 100;
  //n_props += nbins;
  //bin_size = (hb - lb) / (double)n_bins;
  // from x=-6 to x=+6

  // Read initial configuration

  // Evaluate properties of the initial configuration
  Measure();

  return;
}

void Move() {
  int o;
  double p, xold;
  double xnew;
  // Metropolis move

  // Old
  xold = x;

  // New
  xnew = x + delta * (rnd.Rannyu(-1.0, 1.0));

  // Metropolis test
  p = pow((fT(xnew) / fT(xold)), 2);
  if (p >= rnd.Rannyu()) {
    // Update
    x = xnew;
    accepted = accepted + 1.0;
  }
  attempted = attempted + 1.0;

  return;
}

void Measure() // Properties measurement
{
  int bin = 0;
  double v = 0.0;
  int r = 0;
  // reset the hystogram
  //for (int k = igofr; k < igofr + nbins; ++k)
  //  walker[k] = 0.0;

  // Bin

  //r = igofr + int((x - lb) / bin_size);
  //walker[r]++;

  // Hamiltonian
  v = (-0.5 * SecondDerivative(x)) / fT(x) + Potential(x);
  walker[iv] = v;

  return;
}

void Reset(int iblk) // Reset block averages
{

  if (iblk == 1) {
    for (int i = 0; i < n_props; ++i) {
      glob_av[i] = 0;
      glob_av2[i] = 0;
    }
  }

  for (int i = 0; i < n_props; ++i) {
    blk_av[i] = 0;
  }
  blk_norm = 0;
  attempted = 0;
  accepted = 0;
}

void Accumulate(void) // Update block averages
{

  for (int i = 0; i < n_props; ++i) {
    blk_av[i] = blk_av[i] + walker[i];
  }
  blk_norm = blk_norm + 1.0;
}

void Averages(int iblk) // Print results for current block
{
  double r, gdir;
  ofstream Gofr, Energy;
  double stima_g[igofr + nbins];
  double err_g[igofr + nbins];
  cout << "Block number " << iblk << endl;
  cout << "Acceptance rate " << accepted / attempted << endl << endl;

  Energy.open(folder + "output_energy.dat", ios::app);
  //Gofr.open(folder + "output_gofr.dat", ios::app);

  stima_Ham = blk_av[iv] / blk_norm;
  glob_av[iv] += stima_Ham;
  glob_av2[iv] += stima_Ham * stima_Ham;
  err_Ham = Error(glob_av[iv], glob_av2[iv], iblk);
/*not needed, i do it in python
  for (int i = igofr; i < igofr + nbins; i++) {
    stima_g[i] = blk_av[i] / (bin_size * nstep);
    glob_av[i] += stima_g[i];
    glob_av2[i] += stima_g[i] * stima_g[i];
    err_g[i] = Error(glob_av[i], glob_av2[i], iblk);
  }
*/
  if (temp == 0.0) {
    /*if (iblk == nblk) {
      for (int i = igofr; i < igofr + nbins; i++) {
        r = bin_size * (i - igofr) + lb;
        Gofr << r << " " << stima_g[i] << " " << glob_av[i] / iblk << " "
             << err_g[i] << endl;
      }
    }
    */
    Energy << iblk << " " << stima_Ham << " " << glob_av[iv] / (double)iblk
           << " " << err_Ham << endl;
  }

  cout << "----------------------------" << endl << endl;

  Energy.close();
 // Gofr.close();
}

double Error(double sum, double sum2, int iblk) {
  if (iblk == 1) {
    return 0;
  } else {
    return sqrt(fabs(sum2 / (double)iblk - pow(sum / (double)iblk, 2)) /
                (double)iblk);
  }
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
