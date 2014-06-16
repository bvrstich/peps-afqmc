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

/**
 * @author Brecht Verstichel
 * @date 29-04-2014\n\n
 * This class contains some global variables used throughout the program
 */
class Global {

   public:

      static void init(int,int,int);

      template<typename T>
         static T rgen();

      //!static Lattice object containing the info about the lattice
      static Random RN;

      //!x dimension of the lattice, nr of cols
      static int Lx;

      //!y dimension of the lattice, nr of rows
      static int Ly;

      //!physical dimension of sites
      static int d;

      //!Sx matrix
      static TArray<complex<double>,2> Sx;

      //!Sy matrix
      static TArray<complex<double>,2> Sy;

      //!Sz matrix
      static TArray<complex<double>,2> Sz;

};

#endif

/* vim: set ts=3 sw=3 expandtab :*/
