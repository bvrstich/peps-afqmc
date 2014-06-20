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

namespace global{

   //global variables...
   TArray<complex<double>,2> Sx;

   TArray<complex<double>,2> Sy;

   TArray<complex<double>,2> Sz;

   std::vector< std::vector< complex<double> > > auxvec;

   std::vector< TArray<complex<double>,2> > Mx;

   std::vector< TArray<complex<double>,2> > My;

   std::vector< TArray<complex<double>,2> > Mz;

   int Lx;
   int Ly;

   int d;

   Random RN;

   /**
    * @param d_in physical dimension
    * @param Lx_in x dimension of the square lattice
    * @param Ly_in y dimension of the square lattice
    */
   void init(int d_in,int Lx_in,int Ly_in){

      Lx = Lx_in;
      Ly = Ly_in;

      d = d_in;

      //Sx
      Sx.resize(d,d);

      Sx(0,0) = complex<double>(0.0,0.0);
      Sx(0,1) = complex<double>(0.5,0.0);
      Sx(1,0) = complex<double>(0.5,0.0);
      Sx(1,1) = complex<double>(0.0,0.0);

      //Sy
      Sy.resize(d,d);

      Sy(0,0) = complex<double>(0.0,0.0);
      Sy(0,1) = complex<double>(0.0,0.5);
      Sy(1,0) = complex<double>(0.0,-0.5);
      Sy(1,1) = complex<double>(0.0,0.0);

      //Sz
      Sz.resize(d,d);

      Sz(0,0) = complex<double>(-0.5,0.0);
      Sz(0,1) = complex<double>(0.0,0.0);
      Sz(1,0) = complex<double>(0.0,0.0);
      Sz(1,1) = complex<double>(0.5,0.0);

      auxvec.resize(Lx*Ly);

      for(int i = 0;i < auxvec.size();++i)//for x,y and z components
         auxvec[i].resize(3);

      //now construct the propagating operators

      //Sx
      TArray<double,1> eig(d);
      TArray<complex<double>,2> U(Sx);

      lapack::heev(CblasRowMajor, 'V', 'U', d, U.data(), d, eig.data());

      Mx.resize(d);

      for(int i = 0;i < d;++i){

         Mx[i].resize(d,d);

         for(int j = 0;j < d;++j)
            for(int k = 0;k < d;++k)
               Mx[i](j,k) = U(j,i) * std::conj(U(k,i));

      }

      //Sy
      U = Sy;
      lapack::heev(CblasRowMajor, 'V', 'U', d, U.data(), d, eig.data());

      My.resize(d);

      for(int i = 0;i < d;++i){

         My[i].resize(d,d);

         for(int j = 0;j < d;++j)
            for(int k = 0;k < d;++k)
               My[i](j,k) = U(j,i) * std::conj(U(k,i));

      }

      //Sz
      Mz.resize(d);

      for(int i = 0;i < d;++i){

         Mz[i].resize(d,d);
         Mz[i] = 0.0;

         Mz[i](i,i) = 1.0;

      }

   }

   //!function which generates random complex numbers uniformly on a square of side 2
   template<>
      complex<double> rgen(){ 

         return complex<double>(2.0*RN() - 1.0,2.0*RN() - 1.0); 

      }

   //!function which generates uniform random numbers between [-1:1]
   template<>
      double rgen(){ 

         return 2.0*RN() - 1.0;

      }

}
