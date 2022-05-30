

#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <cmath>
#include <vector>

using namespace std;

double error(vector <double> AV,vector <double> AV2,int n){
	
    if (n==0){
        return 0;
	    }
    else
    	{
        return sqrt((AV2[n] - AV[n]*AV[n])/n);}
		}
      //L is the number of throws
void Print(string a,int L,vector<vector<double>> distribution){

   ofstream output(a);
   for(int i = 0; i <L;i++){
      for(int t=0;t<4;t++){
         output <<distribution[i][t]<<" ";
         }
   output<<endl;
   }
   output.close();
}

int main (int argc, char *argv[]){

   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;


//start exercise 1.2
	
	
          
//int S=100;        //number of blocks       
int L=10000;    // throws inside each block        
//int M=L*S;    //total throws
	
vector <int> N {1,2,10,100};

//iterazioni 10000 di S_N
//j from 0 to 1 or 2 or 10 or 100
//
vector<vector<double>> S_N_uniform(L,vector<double>(4));
vector<vector<double>> S_N_Exponential(L,vector<double>(4));
vector<vector<double>> S_N_Lorentzian(L,vector<double>(4));

for(int t=0;t<4;t++){
	for(int i=0;i<L;i++){
		double sum_U=0;
		double sum_E=0;
		double sum_L=0;
		for(int j=0;j<N[t];j++){
			
			sum_U += rnd.Rannyu(0,7);		
			sum_E += rnd.Exponential(1.);		
			sum_L += rnd.Cauchy_Lorentz(1.,0.);		

		}

	sum_U=sum_U/N[t];
	S_N_uniform[i][t]=sum_U;

	sum_E=sum_E/N[t];
	S_N_Exponential[i][t]=sum_E;
		
	sum_L=sum_L/N[t];
	S_N_Lorentzian[i][t]=sum_L;
	
}

	}

Print("uniform.dat", L, S_N_uniform);
Print("exponential.dat", L, S_N_Exponential);
Print("Lorentzian.dat", L, S_N_Lorentzian);



   rnd.SaveSeed();
   return 0;
}
