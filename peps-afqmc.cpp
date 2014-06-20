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

   //set trotter terms on Heisenberg model
   double dtau = 0.001;

   Trotter::heisenberg(dtau);

   //read in the trial state
   char filename[200];
   sprintf(filename,"input/%dx%d/D=%d",L,L,D);

   PEPS< complex<double> > peps;
   peps.load(filename);

   Walker walker;

   walker.calc_properties('V',peps);
/*
   for(int k = 0;k < Trotter::n_trot;++k)
      for(int r = 0;r < 3;++r)
         cout << "(" << k << "," << r << ")\t|\t" << walker.gVL()[r*Trotter::n_trot + k] << endl;
*/
   /*
   int Nw = 100;

   AFQMC afqmc(mps,dtau,Nw);
   afqmc.walk(100);
*/
   return 0;

}
