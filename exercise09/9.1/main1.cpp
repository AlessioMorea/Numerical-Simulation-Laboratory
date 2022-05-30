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

#include "random.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <string>

using namespace std;
int main(int argc, char *argv[]) {
  if (argc != 2) {
    cout << "Error" << endl;
    cout << "Syntax: " << argv[0] << " <input_file_folder>" << endl;
    return -1;
  }

  string method = argv[1];

  

  Genetic gen(method);
  gen.CompleteSim(); 
}