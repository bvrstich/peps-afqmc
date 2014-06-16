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
 * construct an empty PEPS object, note: be sure to initialize the Lattice object before calling the constructor
 */
template<typename T>
PEPS<T>::PEPS() : vector< TArray<T,5> >(Global::Lx * Global::Ly) { }

/**
 * construct constructs a standard PEPS object, note: be sure to initialize the Lattice object before calling the constructor
 * @param D_in cutoff virtual dimension
 */
template<typename T>
PEPS<T>::PEPS(int D_in) : vector< TArray<T,5> >(Global::Lx * Global::Ly) {

   D = D_in;

   int Lx = Global::Lx;
   int Ly = Global::Ly;
   int d = Global::d;

   //corners first

   //r == 0 : c == 0
   (*this)[ 0 ].resize(d,1,D,1,D);

   //r == 0 : c == L - 1
   (*this)[ Lx - 1 ].resize(d,D,D,1,1);

   //r == L - 1 : c == 0
   (*this)[ (Ly-1)*Lx ].resize(d,1,1,D,D);

   //r == L - 1 : c == L - 1
   (*this)[ (Ly-1)*Lx + Lx - 1 ].resize(d,D,1,D,1);

   //sides:

   //r == 0
   for(int c = 1;c < Lx - 1;++c)
      (*this)[ c ].resize(d,D,D,1,D);

   //r == Ly - 1
   for(int c = 1;c < Lx - 1;++c)
      (*this)[ (Ly-1)*Lx + c ].resize(d,D,1,D,D);

   //c == 0
   for(int r = 1;r < Ly - 1;++r)
      (*this)[ r*Lx ].resize(d,1,D,D,D);

   //c == Lx - 1
   for(int r = 1;r < Ly - 1;++r)
      (*this)[ r*Lx + Lx - 1 ].resize(d,D,D,D,1);

   //the rest is full
   for(int r = 1;r < Ly - 1;++r)
      for(int c = 1;c < Lx - 1;++c)
         (*this)[ r*Lx + c ].resize(d,D,D,D,D);

   //now initialize with random numbers
   for(int r = 0;r < Ly;++r)
      for(int c = 0;c < Lx;++c){

         (*this)[ r*Lx + c ].generate(Global::rgen<T>);

         Normalize((*this)[ r*Lx + c ]);
         Scal((T)D,(*this)[ r*Lx + c ]);

      }

}

/**
 * copy constructor
 */
template<typename T>
PEPS<T>::PEPS(const PEPS<T> &peps_copy) : vector< TArray<T,5> >(peps_copy) {

   D = peps_copy.gD();

}

/**
 * empty destructor
 */
template<typename T>
PEPS<T>::~PEPS(){ }

/**
 * access to the individual tensors: const version
 * @param r row index
 * @param c col index
 * @return the tensor on site (r,c)
 */
template<typename T>
const TArray<T,5> &PEPS<T>::operator()(int r,int c) const {

   return (*this)[r*Global::Lx + c];

}

/**
 * access to the individual tensors: const version
 * @param r row index
 * @param c col index
 * @return the tensor on site (r,c)
 */
template<typename T>
TArray<T,5> &PEPS<T>::operator()(int r,int c) {

   return (*this)[r*Global::Lx + c];

}

/**
 * @return the cutoff virutal dimension
 */
template<typename T>
int PEPS<T>::gD() const {

   return D;

}

/**
 * @param D_in value to the D to
 */
template<typename T>
void PEPS<T>::sD(int D_in){

   this->D = D_in;

}

/**
 * @param mpx will be constructed from file
 * @param filename name of the file
 * load the MPX object from a file in binary format.
 */
template<typename T>
void PEPS<T>::load(const char *filename){

   int Lx = Global::Lx;
   int Ly = Global::Ly;

   for(int row = 0;row < Ly;++row)
      for(int col = 0;col < Lx;++col){

         char name[200];

         sprintf(name,"%s/site_(%d,%d).peps",filename,row,col);

         std::ifstream fin(name);

         int Da,Db,Dc,Dd,De;

         fin >> Da >> Db >> Dc >> Dd >> De;

         //different order!
         (*this)(row,col).resize(Dc,Da,Db,Dd,De);

         for(int a = 0;a < Da;++a)
            for(int b = 0;b < Db;++b)
               for(int c = 0;c < Dc;++c)
                  for(int d = 0;d < Dd;++d)
                     for(int e = 0;e < De;++e)
                        fin >> a >> b >> c >> d >> e >> (*this)(row,col)(c,a,b,d,e);

   }

}

//forward declarations for types to be used!
template PEPS<double>::PEPS();
template PEPS< complex<double> >::PEPS();

template PEPS<double>::PEPS(int);
template PEPS< complex<double> >::PEPS(int);

template PEPS<double>::PEPS(const PEPS<double> &);
template PEPS< complex<double> >::PEPS(const PEPS< complex<double> > &);

template PEPS<double>::~PEPS();
template PEPS< complex<double> >::~PEPS();

template TArray<double,5> &PEPS<double>::operator()(int r,int c);
template TArray<complex<double>,5> &PEPS< complex<double> >::operator()(int r,int c);

template const TArray<double,5> &PEPS<double>::operator()(int r,int c) const;
template const TArray<complex<double>,5> &PEPS< complex<double> >::operator()(int r,int c) const;

template int PEPS<double>::gD() const;
template int PEPS< complex<double> >::gD() const;

template void PEPS<double>::sD(int);
template void PEPS< complex<double> >::sD(int);

template void PEPS<double>::load(const char *filename);
template void PEPS< complex<double> >::load(const char *filename);
