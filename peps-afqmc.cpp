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

int main(int argc,char *argv[]){

   cout.precision(15);
   srand(time(NULL));

   int L = atoi(argv[1]);
   int d = atoi(argv[2]);
   int D = atoi(argv[3]);
   int D_aux = atoi(argv[4]);

   //initialize the dimensions of the problem
   global::init(d,L,L);

   //and some static objects
   Environment::init(D,D_aux);

   //read in the trial state
   char filename[200];
   sprintf(filename,"input/%dx%d/D=%d",L,L,D);

   PEPS< complex<double> > peps;
   peps.load(filename);

   Walker walker(10);

   Environment::calc_env('A',peps,walker);
   Environment::test_env();
/*
   walker.calc_properties('H',peps);

   double dtau = 0.01;
   int Nw = 100;

   AFQMC afqmc(mps,dtau,Nw);
   afqmc.walk(100);
*/
   return 0;

}
