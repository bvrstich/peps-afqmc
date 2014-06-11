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

   //read in the trial state
   char filename[200];
   sprintf(filename,"input/%dx%d/D=%d.mps",L,L,D);

   //intialize some storage

   /*
   double dtau = 0.01;
   int Nw = 100;

   AFQMC afqmc(mps,dtau,Nw);
   afqmc.walk(100);
   */

   return 0;

}
