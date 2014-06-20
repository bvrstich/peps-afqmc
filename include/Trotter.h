#ifndef TROTTER_H
#define TROTTER_H

#include <iostream>
#include <iomanip>
#include <complex>
#include <vector>

class Trotter {

   public:
   
      //Constructor
      static void heisenberg(double);
      
      static double gdtau();
      
      static const TArray<complex<double>,2> &gV();
      
      //!the number of non-zero eigenvalues of the coupling matrix: equals the number of trotter product terms
      static int n_trot;

   private:
   
      //!the transformation matrix
      static TArray<complex<double>,2> V;


      //!The time step
      static double dtau;

};

#endif
