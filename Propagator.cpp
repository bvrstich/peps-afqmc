#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include <fstream>

using std::ostream;
using std::ofstream;
using std::ifstream;
using std::cout;
using std::endl;
using std::ios;
using std::vector;
using std::complex;

#include "include.h"

using namespace global;

/**
 * standard constructor: sets the size of vector, and initializes the ZArray objects to the correct quantumnumbers and dimensions
 * @param L size
 */
Propagator::Propagator() : vector< TArray<complex<double>,2> > (Lx*Ly) {

   x = 0.0;
   k = 0;
   r = 0;
  
   for(int i = 0;i < Lx*Ly;++i)
      (*this)[i].resize(d,d);
 
}

/**
 * copy constructor
 * @param prop_copy input Propagator object
 */
Propagator::Propagator(const Propagator &prop_copy) : vector< TArray<complex<double>,2> > (prop_copy) {

   x = prop_copy.gx();
   k = prop_copy.gk();
   r = prop_copy.gr();
   
}

/**
 * destructor
 */
Propagator::~Propagator() { }

/**
 * @return the auxiliary field variable x
 */
complex<double> Propagator::gx() const {

   return x;

}

/**
 * @return the trotter index k
 */
int Propagator::gk() const {

   return k;

}

/**
 * @return the type of operator r (0=x,1=y,2=z)
 */
int Propagator::gr() const {

   return r;

}

/**
 * set the variables to fill the propagator on
 */
void Propagator::set(complex<double> x_in,int k_in,int r_in){

   x = x_in;
   k = k_in;
   r = r_in;

}

/**
 * fill the propagator using a Trotter input object
 */
void Propagator::fill() {

   if(r == 0){//x

      for(int row = 0;row < Ly;++row)
         for(int col = 0;col < Lx;++col){

            blas::copy(d*d, Mx[0].data(), 1, (*this)[row*Lx + col].data(), 1);

            double m = -0.5 * ( d - 1.0 );

            blas::scal(d*d, exp(x * Trotter::gV()(k,row*Lx + col) * m ) , (*this)[row*Lx + col].data(), 1);

            for(int i = 1;i < d;++i){

               m++;

               blas::axpy(d*d, exp(  x * Trotter::gV()(k,row*Lx + col) * m ) , Mx[i].data(), 1, (*this)[row*Lx + col].data(), 1);

            }

         }

   }
   else if(r == 1){//y

      for(int row = 0;row < Ly;++row)
         for(int col = 0;col < Lx;++col){

            blas::copy(d*d, My[0].data(), 1, (*this)[row*Lx + col].data(), 1);

            double m = -0.5 * ( d - 1.0 );

            blas::scal(d*d, exp(x * Trotter::gV()(k,row*Lx + col) * m ) , (*this)[row*Lx + col].data(), 1);

            for(int i = 1;i < d;++i){

               m++;

               blas::axpy(d*d, exp(  x * Trotter::gV()(k,row*Lx + col) * m ) , My[i].data(), 1, (*this)[row*Lx + col].data(), 1);

            }

         }

   }
   else{//z

      for(int row = 0;row < Ly;++row)
         for(int col = 0;col < Lx;++col){

            blas::copy(d*d, Mz[0].data(), 1, (*this)[row*Lx + col].data(), 1);

            double m = -0.5 * ( d - 1.0 );

            blas::scal(d*d, exp(x * Trotter::gV()(k,row*Lx + col) * m ) , (*this)[row*Lx + col].data(), 1);

            for(int i = 1;i < d;++i){

               m++;

               blas::axpy(d*d, exp(  x * Trotter::gV()(k,row*Lx + col) * m ) , Mz[i].data(), 1, (*this)[row*Lx + col].data(), 1);

            }

         }

   }

}
