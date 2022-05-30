#include <mpi.h>
#include "Genetic.h"
#include "random.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <string>

using namespace std;
int main(int argc, char *argv[]) {
    int size, rank;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    
  string method = "America";
  

  Genetic gen(method,rank,size);
  gen.CompleteSim(); 
  MPI_Finalize();
}
