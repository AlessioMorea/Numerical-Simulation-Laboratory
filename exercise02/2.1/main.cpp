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
	    	}else{
        	return sqrt((AV2[n] - AV[n]*AV[n])/n);}
	}


double func(double x){
	return M_PI*cos(M_PI*x/2)/2;
	}
//p(x)=pi/2*(1-x)
//F(x)=pi/2*(x-x^2/2)
//Inverse F(x): x=1-sqrt(1-y)

double InvCDF(double y) {
	
	if(y<=0){
		cout<<"ERROR: random variable is not positive";
		return -1;
		
	}
	return 	1.-sqrt(1-y);
}

double pdf(double x){
	return 2.*(1-x);
}




void Print(string a, double N,vector <double> sum_prog,vector <double> err_prog){
	
	ofstream output(a);
    	for(int i = 0; i <N;i++){
        	output << i << " " << sum_prog[i] << " " << err_prog[i] <<endl;
		}
    output.close();
}

void Data_Blocking_Method(int N,int L,vector<double> f,string a){
	
vector <double> ave(N);
vector <double> av2(N);
vector <double> sum_prog(N);
vector <double> su2_prog(N);
vector <double> err_prog(N);
	

for (int i=0;i<N;i++){
	
	double sum=0;
	for(int j=0;j<L;j++){
        	int k = j+i*L;
        	sum += f[k];
		}

	ave[i]=sum/L;
	av2[i]=(ave[i]*ave[i]);
	}
	
for (int i=0;i<N;i++){
	for(int j=0;j<i+1;j++){
        	sum_prog[i] += ave[j]; 
        	su2_prog[i] += av2[j];
	}
	
	sum_prog[i]=sum_prog[i]/(i+1);	
    	su2_prog[i]/=(i+1);
    	err_prog[i] = error(sum_prog,su2_prog,i);
	
	}

Print(a,N,sum_prog,err_prog);
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
rnd.SaveSeed();

//start exercise 2.1
	
int M=1000000;    //total throws          
int N=100;        //number of blocks       
int L=int(M/N);    // throws inside each block        
    
vector <double> r(M);
vector <double> f(M);
vector <double> g(M);  //importance sampling

	
for(int i=0;i<M;i++){
	
	r[i]=rnd.Rannyu();
	f[i]=func(r[i]);
	g[i]=func(InvCDF(r[i]))/pdf(InvCDF(r[i]));
	
}
	
Data_Blocking_Method(N,L,f,"standard.dat");
Data_Blocking_Method(N,L,g,"sampling.dat");

return 0;
}
