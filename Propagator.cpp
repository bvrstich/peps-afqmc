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
 /*
void Propagator::fill(const Trotter &trotter) {

   if(r == 0){//x

      for(int site = 0;site < L;++site){

         blas::copy(size, trotter.gMx(0).data(), 1, (*this)[site].data(), 1);

         double m = -0.5 * ( d - 1.0 );

         blas::scal(size, exp(x * trotter.gV()(k,site) * m ) , (*this)[site].data(), 1);

         for(int i = 1;i < d;++i){

            m++;

            blas::axpy(size, exp(  x * trotter.gV()(k,site) * m ) , trotter.gMx(i).data(), 1, (*this)[site].data(), 1);

         }

      }

   }
   else if(r == 1){//y

      for(int site = 0;site < L;++site){

         blas::copy(size, trotter.gMy(0).data(), 1, (*this)[site].data(), 1);

         double m = -0.5 * ( d - 1.0 );

         blas::scal(size, exp(x * trotter.gV()(k,site) * m ) , (*this)[site].data(), 1);

         for(int i = 1;i < d;++i){

            m++;

            blas::axpy(size, exp(  x * trotter.gV()(k,site) * m ) , trotter.gMy(i).data(), 1, (*this)[site].data(), 1);

         }

      }

   }
   else{//z

      for(int site = 0;site < L;++site){

         blas::copy(size, trotter.gMz(0).data(), 1, (*this)[site].data(), 1);

         double m = -0.5 * ( d - 1.0 );

         blas::scal(size, exp(x * trotter.gV()(k,site) * m ) , (*this)[site].data(), 1);

         for(int i = 1;i < d;++i){

            m++;

            blas::axpy(size, exp(  x * trotter.gV()(k,site) * m ) , trotter.gMz(i).data(), 1, (*this)[site].data(), 1);

         }

      }

   }
}
*/
