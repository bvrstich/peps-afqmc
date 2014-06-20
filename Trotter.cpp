#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <complex>
#include <vector>

using std::cout;
using std::endl;

#include "include.h"

using namespace global;

//!the transformation matrix
TArray<complex<double>,2> Trotter::V;

//!the number of non-zero eigenvalues of the coupling matrix: equals the number of trotter product terms
int Trotter::n_trot;

//!The time step
double Trotter::dtau;

/**
 * constructor, diagonalizes the coupling matrix J to make the auxiliary field operators V(r,k)
 * @param dtau time step
 */
void Trotter::heisenberg(double dtau_in){

   dtau = dtau_in;

   //set the heisenberg coupling matrix
   TArray<double,2> J(Lx*Ly,Lx*Ly);
   J = 0.0;

   //horizontal interaction
   for(int r = 0;r < Ly;++r)
      for(int c = 0;c < Lx - 1;++c)
         J(r*Lx + c,r*Lx + c + 1) = 1.0;

   //vertical interaction
   for(int r = 0;r < Ly-1;++r)
      for(int c = 0;c < Lx;++c)
         J(r*Lx + c,(r + 1)*Lx + c) = 1.0;

   //diagonalize the coupling matrix
   TArray<double,1> eig(Lx*Ly);

   lapack::syev(CblasRowMajor, 'V', 'U', Lx*Ly, J.data(), Lx*Ly, eig.data());

   n_trot = 0;

   for(int i = 0;i < Lx*Ly;++i){

      if(fabs(eig(i)) > 1.0e-14)
         n_trot++;

   }

   //now transform the elements with dtau and the eigenvalues for the propagator to form the transformation V
   V.resize(n_trot,Lx*Ly);

   int k = 0;

   for(int i = 0;i < Lx*Ly;++i){

      if(fabs(eig(i)) > 1.0e-14){

         complex<double> tmp =  std::sqrt( (complex<double>)-eig(i) * dtau);

         for(int j = 0;j < Lx*Ly;++j)
            V(k,j) = tmp * J(j,i);

         ++k;

      }

   }

}

/**
 * @return the timestep
 */
double Trotter::gdtau() {

   return dtau;

}

/**
 * @return the auxiliary field matrix
 */
const TArray<complex<double>,2> &Trotter::gV() {

   return V;

}
