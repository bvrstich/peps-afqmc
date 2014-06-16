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

/**
 * empty constructor
 */
MPO::MPO() : vector< TArray<complex<double>,4> > ( Global::lat.gLx() ) { }

/**
 * constructor of MPO: allocates the matrices
 * @param D_in bond input dimension
 */
MPO::MPO(int D_in) : vector< TArray<complex<double>,4> > ( Global::lat.gLx() ) {

   this->D = D_in;

   int Lx = Global::lat.gLx();

   (*this)[0].resize(1,D,D,D);

   for(int c = 1;c < Lx - 1;++c)
      (*this)[c].resize(D,D,D,D);

   (*this)[Lx-1].resize(D,D,D,1);

}

/**
 * constructs a standard MPO object, by creating a single layer from the constraction of a peps with a physical vector
 * @param option 'V' == vertical stripe, 'H' == horizontal stripe
 * @param rc row or column index
 */
void MPO::fill(char option,int rc,const PEPS< complex<double> > &peps,const Walker &walker) {

   int Lx = Global::lat.gLx();
   int Ly = Global::lat.gLy();

   int d = Global::lat.gd();

   complex<double> one(1.0,0.0);
   complex<double> zero(0.0,0.0);

   if(option == 'H'){//horizontal

      for(int c = 0;c < Lx;++c){

         int dim = (*this)[c].size();

         blas::gemv(CblasRowMajor, CblasTrans, d, dim , one, peps(rc,c).data(), dim , walker(rc,c).data(), 1, zero, (*this)[c].data(), 1);

      }

   }
   else{//vertical!

      TArray<complex<double>,4> tmp;

      for(int r = 0;r < Ly;++r){

         int dim = (*this)[r].size();
         
         tmp.resize(peps(r,rc).shape(1),peps(r,rc).shape(2),peps(r,rc).shape(3),peps(r,rc).shape(4));

         blas::gemv(CblasRowMajor, CblasTrans, d, dim , one, peps(r,rc).data(), dim , walker(r,rc).data(), 1, zero, tmp.data(), 1);

         //PERMUTE!!
         Permute(tmp,shape(2,3,0,1),(*this)[r]);

      }


   }

}

/**
 * copy constructor
 */
MPO::MPO(const MPO &mpo_copy) : vector< TArray< complex<double> ,4> >(mpo_copy) {

   D = mpo_copy.gD();

}

/**
 * empty destructor
 */
MPO::~MPO(){ }

/**
 * @return virtual dimension of the MPO
 */
int MPO::gD() const {

   return D;

}
