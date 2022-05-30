/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include "Monte_Carlo_ISING_1D.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <string>

using namespace std;
string folder;
int main(int argc, char *argv[])
{
  if (argc != 2)
  {
    cout << "Error" << endl;
    cout << "Syntax: " << argv[0] << " <input_file_folder>" << endl;
    return -1;
  }

  folder = argv[1];
  folder = "input/" + folder + "/";

  Input(); // Inizialization

  while (temp > 0.49)
  {
    beta = 1.0 / temp;

    cout << "Starting Simulation for T=" << temp << endl;

    for (int i = 0; i < nequilibration; i++)
    {
      Move(metro);
    }
    for (int iblk = 1; iblk <= nblk; ++iblk) // Simulation
    {
      Reset(iblk); // Reset block averages
      for (int istep = 1; istep <= nstep; ++istep)
      {
        Move(metro);
        Measure();
        Accumulate(); // Update block averages
      }
      Averages(iblk); // Print results for current block
    }
    ConfFinal(); // Write final configuration

    temp -= 0.05;
  }
  return 0;
}

void Input(void)
{
  ifstream ReadInput, ReadConf;

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl
       << endl;
  cout << "Nearest neighbour interaction      " << endl
       << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl
       << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

  // Read seed for random numbers
  int p1, p2;
  ifstream Primes("Primes");
  Primes >> p1 >> p2;
  Primes.close();

  ifstream input("seed.in");
  input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  rnd.SetRandom(seed, p1, p2);
  input.close();

  // Read input informations
  ReadInput.open(folder + "input.in");
  ReadInput >> restart;
  ReadInput >> temp;
  beta = 1.0 / temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl
       << endl;

  ReadInput >> metro; // if=1 Metropolis else Gibbs

  ReadInput >> nblk;

  ReadInput >> nstep;

  ReadInput >> nequilibration;

  if (metro == 1)
    cout << "The program perform Metropolis moves" << endl;
  else
    cout << "The program perform Gibbs moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl;
  cout << "Number of equilibration steps  = " << nequilibration << endl
       << endl;

  ReadInput.close();

  // Prepare arrays for measurements
  iu = 0; // Energy
  ic = 1; // Heat capacity
  im = 2; // Magnetization
  ix = 3; // Magnetic susceptibility

  n_props = 4; // Number of observables
  if (restart)
  {
    ReadConf.open("config.final");
    for (int i = 0; i < nspin; ++i)
    {
      ReadConf >> s[i];
    }
    ReadConf.close();
  }
  else
  {
    // initial configuration
    for (int i = 0; i < nspin; ++i)
    {
      if (rnd.Rannyu() >= 0.5)
        s[i] = 1;
      else
        s[i] = -1;
    }
  }
  // Evaluate energy etc. of the initial configuration
  Measure();

  // Print initial values for the potential energy and virial
  cout << "Initial energy = " << walker[iu] / (double)nspin << endl;
}

void Move(int metro)
{
  int o;
  double p, energy_diff; // energy_old, energy_new, sm,;
  // double energy_up, energy_down;

  for (int i = 0; i < nspin; ++i)
  {
    // Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    o = (int)(rnd.Rannyu() * nspin);

    if (metro == 1) // Metropolis
    {
      // INCLUDE YOUR CODE HERE
      // Difference
      energy_diff = 2. * J * s[o] * (s[Pbc(o - 1)] + s[Pbc(o + 1)]) + 2. * h * s[o];

      // Metropolis test
      if (energy_diff > 0)
      {
        p = exp(-beta * (energy_diff));
        if (p >= rnd.Rannyu())
        {
          // Update
          s[o] = -s[o];
          accepted = accepted + 1.0;
        }
        attempted = attempted + 1.0;
      }
      else
      {
        s[o] = -s[o];
        accepted = accepted + 1.0;
        attempted = attempted + 1.0;
      }
    }
    else // Gibbs sampling
    {
      accepted = accepted + 1.0;
      attempted = attempted + 1.0;

      // INCLUDE YOUR CODE HERE
      s[o] = 1;
      energy_diff = 2. * J * s[o] * (s[Pbc(o - 1)] + s[Pbc(o + 1)]) + 2. * h * s[o];
      p = exp(-beta * (energy_diff));
      if (1.0 / (1.0 + p) < rnd.Rannyu())
      {
        s[o] = -1;
      }
    }
  }
}

double Boltzmann(int sm, int ip)
{
  double ene = -J * sm * (s[Pbc(ip - 1)] + s[Pbc(ip + 1)]) - h * sm;
  return ene;
}

void Measure()
{
  // int bin;
  double u = 0.0, m = 0.0;

  // cycle over spins
  for (int i = 0; i < nspin; ++i)
  {
    u += -J * s[i] * s[Pbc(i + 1)] - 0.5 * h * (s[i] + s[Pbc(i + 1)]);
    // INCLUDE YOUR CODE HERE
    m += s[i];
  }
  walker[iu] = u;
  // INCLUDE YOUR CODE HERE
  walker[ic] = u * u;
  walker[im] = m;
  walker[ix] = m * m;
}

void Reset(int iblk) // Reset block averages
{

  if (iblk == 1)
  {
    for (int i = 0; i < n_props; ++i)
    {
      glob_av[i] = 0;
      glob_av2[i] = 0;
    }
  }

  for (int i = 0; i < n_props; ++i)
  {
    blk_av[i] = 0;
  }
  blk_norm = 0;
  attempted = 0;
  accepted = 0;
}

void Accumulate(void) // Update block averages
{

  for (int i = 0; i < n_props; ++i)
  {
    blk_av[i] = blk_av[i] + walker[i];
  }
  blk_norm = blk_norm + 1.0;
}

void Averages(int iblk) // Print results for current block
{

  ofstream Ene, Heat, Mag, Chi;
  const int wd = 12;

  cout << "Block number " << iblk << endl;
  cout << "Acceptance rate " << accepted / attempted << endl
       << endl;

  stima_u = blk_av[iu] / blk_norm / (double)nspin; // Energy
  glob_av[iu] += stima_u;
  glob_av2[iu] += stima_u * stima_u;
  err_u = Error(glob_av[iu], glob_av2[iu], iblk);
  if (temp == 2 && metro == 1 && h == 0.0)
  {
    Ene.open("output.ene.0", ios::app);
    Ene << setw(wd) << iblk << setw(wd) << stima_u << setw(wd) << glob_av[iu] / (double)iblk << setw(wd) << err_u << endl;
    Ene.close();
  }
  // INCLUDE YOUR CODE HERE

  // Heat.open("output.heat.0", ios::app);

  stima_c = beta * beta * (blk_av[ic] / blk_norm / (double)nspin - (double)nspin * stima_u * stima_u); // Heat Capacity
  glob_av[ic] += stima_c;
  glob_av2[ic] += stima_c * stima_c;
  err_c = Error(glob_av[ic], glob_av2[ic], iblk);

  // Heat << setw(wd) << iblk << setw(wd) << stima_c << setw(wd) << glob_av[ic] / (double)iblk << setw(wd) << err_c << endl;

  // Heat.close();

  // Chi.open("output.chi.0", ios::app);
  stima_x = beta * (blk_av[ix] / blk_norm / (double)nspin); // susceptibility
  glob_av[ix] += stima_x;
  glob_av2[ix] += stima_x * stima_x;

  err_x = Error(glob_av[ix], glob_av2[ix], iblk);
  // Chi << setw(wd) << iblk << setw(wd) << stima_x << setw(wd)<< glob_av[ix] / (double)iblk << setw(wd) << err_x << endl;

  // Chi.close();

  // Mag.open("output.mag.0", ios::app);

  stima_m = (blk_av[im] / blk_norm / (double)nspin); // Magnetization
  glob_av[im] += stima_m;
  glob_av2[im] += stima_m * stima_m;
  err_m = Error(glob_av[im], glob_av2[im], iblk);

  // Mag << setw(wd) << iblk << setw(wd) << stima_m << setw(wd)<< glob_av[im] / (double)iblk << setw(wd) << err_m << endl;
  // Mag.close();

  if (iblk == nblk)
  {
    string type;
    if (metro == 1)
    {
      type = "metro";
    }
    else
    {
      type = "gibbs";
    }
    if (h == 0.0)
    {
      Heat.open("results/" + type + "/output.heat.0", ios::app);
      Heat << setw(wd) << temp << setw(wd) << stima_c << setw(wd) << glob_av[ic] / (double)iblk << setw(wd) << err_c << endl;
      Heat.close();
      Chi.open("results/" + type + "/output.chi.0", ios::app);
      Chi << setw(wd) << temp << setw(wd) << stima_x << setw(wd) << glob_av[ix] / (double)iblk << setw(wd) << err_x << endl;
      Chi.close();
      Ene.open("results/" + type + "/output.ene.0", ios::app);
      Ene << setw(wd) << temp << setw(wd) << stima_u << setw(wd) << glob_av[iu] / (double)iblk << setw(wd) << err_u << endl;
      Ene.close();
    }
    if (h == 0.02)
    {
      Mag.open("results/" + type + "/output.mag.0", ios::app);
      Mag << setw(wd) << temp << setw(wd) << stima_m << setw(wd) << glob_av[im] / (double)iblk << setw(wd) << err_m << endl;
      Mag.close();
    }
  }

  cout << "----------------------------" << endl
       << endl;
}

void ConfFinal(void)
{
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl
       << endl;
  WriteConf.open("config.final");
  for (int i = 0; i < nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

int Pbc(int i) // Algorithm for periodic boundary conditions
{
  if (i >= nspin)
    i = i - nspin;
  else if (i < 0)
    i = i + nspin;
  return i;
}

double Error(double sum, double sum2, int iblk)
{
  if (iblk == 1)
    return 0.0;
  else
    return sqrt((sum2 / (double)iblk - pow(sum / (double)iblk, 2)) /
                (double)(iblk - 1));
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
