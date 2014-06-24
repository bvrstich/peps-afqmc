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


   //initialize the dimensions of the problem
   global::init(d,L,L);

   //and some static objects
   Environment::init(D,D_aux);

   //set trotter terms on Heisenberg model
   double dtau = 0.01;

   Trotter::heisenberg(dtau);

   //read in the trial state
   char filename[200];
   sprintf(filename,"input/%dx%d/D=%d",L,L,D);

   PEPS< complex<double> > peps;
   peps.load(filename);

   Walker walker;

   //set the values
   Propagator P;

   //now loop over the auxiliary fields:
 //  for(int k = 0;k < Trotter::n_trot;++k)
 //     for(int r = 0;r < 3;++r){

         double x = 0.5;//RN.normal();

         complex<double> shift(0.0,0.0);// = walker[i].gVL(k,r);

         //set the values
         P.set(x + shift,1,1);

         //and fill the propagator
         P.fill();

         //and apply it to the walker:
         walker.propagate(P);

      //}

   walker.normalize();

   Environment::U.fill('H',peps,walker);
   Environment::calc_env('H',peps,walker);

   Environment::test_env();
/*
   walker.calc_properties('H',peps);

   cout << walker.gEL() << endl;
   cout << walker.gOverlap() << endl;

      int Nw = 1;

      AFQMC afqmc(peps,Nw);
      afqmc.walk(1);
    */
   return 0;

}
