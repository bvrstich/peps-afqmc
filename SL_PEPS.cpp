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

using namespace global;

/** 
 * empty constructor: just sets the length of the vector
 */
SL_PEPS::SL_PEPS() : vector< TArray<complex<double>,4> >( Lx*Ly ) { }

/** 
 * standard constructor: just takes in
 * @param L_in length of the chain
 * @param D_in virtual max bond dimension
 * allocates the tensors and fills them randomly
 */
SL_PEPS::SL_PEPS(int D_in) : vector< TArray<complex<double>,4> >( Lx*Ly ) {

   D = D_in;

   //corners first

   //r == 0 : c == 0
   (*this)[ 0 ].resize(1,D,1,D);

   //r == 0 : c == L - 1
   (*this)[ Lx - 1 ].resize(D,D,1,1);

   //r == L - 1 : c == 0
   (*this)[ (Ly-1)*Lx ].resize(1,1,D,D);

   //r == L - 1 : c == L - 1
   (*this)[ (Ly-1)*Lx + Lx - 1 ].resize(D,1,D,1);

   //sides:

   //r == 0
   for(int c = 1;c < Lx - 1;++c)
      (*this)[ c ].resize(D,D,1,D);

   //r == Ly - 1
   for(int c = 1;c < Lx - 1;++c)
      (*this)[ (Ly-1)*Lx + c ].resize(D,1,D,D);

   //c == 0
   for(int r = 1;r < Ly - 1;++r)
      (*this)[ r*Lx ].resize(1,D,D,D);

   //c == Lx - 1
   for(int r = 1;r < Ly - 1;++r)
      (*this)[ r*Lx + Lx - 1 ].resize(D,D,D,1);

   //the rest is full
   for(int r = 1;r < Ly - 1;++r)
      for(int c = 1;c < Lx - 1;++c)
         (*this)[ r*Lx + c ].resize(D,D,D,D);

}

/**
 * copy constructor
 */
SL_PEPS::SL_PEPS(const SL_PEPS &slp_copy) : vector< TArray<complex<double>,4> >(slp_copy) {

   this->D = slp_copy.gD();

}

/**
 * empty destructor
 */
SL_PEPS::~SL_PEPS(){ }

/**
 * @return virtual dimension of the SL_PEPS
 */
int SL_PEPS::gD() const {

   return D;

}

/**
 * fill the SL_PEPS object by contracting a peps with a walker
 * @param option == 'H' keep regular order of indices
 * @param peps input PEPS<> object
 * @param walker the Walker object
 */
void SL_PEPS::fill(char option,const PEPS< complex<double> > &peps,const Walker &walker){

   complex<double> one(1.0,0.0);
   complex<double> zero(0.0,0.0);

   if(option == 'H'){

      for(int r = 0;r < Ly;++r)
         for(int c = 0;c < Lx;++c){

            int dim = (*this)[r*Lx +c].size();

            blas::gemv(CblasRowMajor, CblasTrans, d, dim , one, peps(r,c).data(), dim , walker(r,c).data(), 1, zero, (*this)[r*Lx + c].data(), 1);

         }

   }
   else{//VERTICAL!

      TArray<complex<double>,4> tmp;

      for(int r = 0;r < Ly;++r)
         for(int c = 0;c < Lx;++c){

            int dim = (*this)[r*Lx +c].size();

            tmp.resize(peps(r,c).shape(1),peps(r,c).shape(2),peps(r,c).shape(3),peps(r,c).shape(4));

            blas::gemv(CblasRowMajor, CblasTrans, d, dim , one, peps(r,c).data(), dim , walker(r,c).data(), 1, zero, tmp.data(), 1);

            //PERMUTE!!
            Permute(tmp,shape(2,3,0,1),(*this)[r*Lx + c]);

         }

   }

}

/**
 * fill the SL_PEPS object by contracting a peps with a walker
 * @param option == 'H' keep regular order of indices
 * @param peps input PEPS<> object
 * @param walker the Walker object
 */
void SL_PEPS::fill(char option,const PEPS< complex<double> > &peps,TArray<complex<double>,2> &O,const Walker &walker){

   complex<double> one(1.0,0.0);
   complex<double> zero(0.0,0.0);

   if(option == 'H'){

      TArray<complex<double>,1> vec(d);

      for(int r = 0;r < Ly;++r)
         for(int c = 0;c < Lx;++c){

            //first act with O on walker
            blas::gemv(CblasRowMajor, CblasTrans, d, d , one, O.data(), d , walker(r,c).data(), 1, zero, vec.data(), 1);

            int dim = (*this)[r*Lx +c].size();

            blas::gemv(CblasRowMajor, CblasTrans, d, dim , one, peps(r,c).data(), dim , vec.data(), 1, zero, (*this)[r*Lx + c].data(), 1);

         }

   }
   else{//VERTICAL!

      TArray<complex<double>,4> tmp;
      TArray<complex<double>,1> vec(d);

      for(int r = 0;r < Ly;++r)
         for(int c = 0;c < Lx;++c){

            //first act with O on walker
            blas::gemv(CblasRowMajor, CblasTrans, d, d , one, O.data(), d , walker(r,c).data(), 1, zero, vec.data(), 1);

            int dim = (*this)[r*Lx +c].size();

            tmp.resize(peps(r,c).shape(1),peps(r,c).shape(2),peps(r,c).shape(3),peps(r,c).shape(4));

            blas::gemv(CblasRowMajor, CblasTrans, d, dim , one, peps(r,c).data(), dim , vec.data(), 1, zero, tmp.data(), 1);

            //PERMUTE!!
            Permute(tmp,shape(2,3,0,1),(*this)[r*Lx + c]);

         }

   }

}

/**
 * access to the individual tensors: const version
 * @param r row index
 * @param c col index
 * @return the tensor on site (r,c)
 */
const TArray<complex<double>,4> &SL_PEPS::operator()(int r,int c) const {

   return (*this)[r*Lx + c];

}

/**
 * access to the individual tensors: const version
 * @param r row index
 * @param c col index
 * @return the tensor on site (r,c)
 */
TArray<complex<double>,4> &SL_PEPS::operator()(int r,int c) {

   return (*this)[r*Lx + c];

}
