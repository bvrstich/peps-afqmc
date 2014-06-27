#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <complex>
#include <vector>

using std::cout;
using std::endl;
using std::ofstream;
using std::ifstream;
using std::complex;
using std::vector;

#include "include.h"

using namespace global;

int main(int argc,char *argv[]){

   cout.precision(15);
   srand(time(NULL));

   int L = atoi(argv[1]);
   int d = atoi(argv[2]);
   int D = atoi(argv[3]);
   int D_aux = atoi(argv[4]);

   //initialize the dimensions of the problem, set the trial
   global::init(D,d,L,L);

   //and some static objects
   Environment::init(D,D_aux);

   //set trotter terms on Heisenberg model
   double dtau = 0.001;

   Trotter::heisenberg(dtau);

   int Nw = 1024;

   AFQMC afqmc(Nw);
   afqmc.walk(1000000);

   return 0;

}
