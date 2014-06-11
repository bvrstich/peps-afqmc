#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <complex>

using std::cout;
using std::endl;
using std::vector;
using std::complex;
using std::ofstream;

#include "include.h"

Lattice Global::lat;

Random Global::RN;

//!function which generates random complex numbers uniformly on a square of side 2
template<>
complex<double> Global::rgen(){ 

   return complex<double>(2.0*RN() - 1.0,2.0*RN() - 1.0); 

}

//!function which generates uniform random numbers between [-1:1]
template<>
double Global::rgen(){ 

   return 2.0*RN() - 1.0;

}
