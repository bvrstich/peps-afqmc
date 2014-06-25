#ifndef GLOBAL_H
#define GLOBAL_H

#include <iostream>
#include <fstream>
#include <vector>
#include <complex>

using std::ostream;
using std::vector;
using std::complex;

using namespace btas;

namespace global {

   extern Random RN;

   //!x dimension of the lattice, nr of cols
   extern int Lx;

   //!y dimension of the lattice, nr of rows
   extern int Ly;

   //!physical dimension of sites
   extern int d;

   //!Sx matrix
   extern TArray<complex<double>,2> Sx;

   //!Sy matrix
   extern TArray<complex<double>,2> Sy;

   //!Sz matrix
   extern TArray<complex<double>,2> Sz;

   //!number of threads in the program: for openmp
   extern int omp_num_threads;

   //!intermediate storage for calculation of auxiliary operator expectation values
   extern vector< vector< vector< complex<double> > > > auxvec;

   void init(int,int,int);

   template<typename T>
      T rgen();
   
   //!eigenvector of Sx |Sx><Sx|
   extern std::vector< TArray<complex<double>,2> > Mx;

   //!eigenvector of Sy |Sy><Sy|
   extern std::vector< TArray<complex<double>,2> > My;

   //!eigenvector of Sz |Sz><Sz|
   extern std::vector< TArray<complex<double>,2> > Mz;

};

#endif

/* vim: set ts=3 sw=3 expandtab :*/
