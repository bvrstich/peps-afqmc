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
Lattice::Lattice(){ }

/**
 * constructor of the lattice object
 * @param Lx x dimension, nr of cols
 * @param Ly y dimension, nr of rows
 * @param d physical dimension of the sites
 */
Lattice::Lattice(int Lx,int Ly,int d){

   this->Lx = Lx;
   this->Ly = Ly;

   this->d = d;

   rc2i.resize(Ly);

   for(int r = 0;r < Ly;++r)
      rc2i[r].resize(Lx);

   vector<int> v(2);

   int i = 0;

   for(int r = 0;r < Ly;++r)
      for(int c = 0;c < Lx;++c){

         v[0] = r;
         v[1] = c;

         i2rc.push_back(v);

         rc2i[r][c] = i;

         ++i;

      }

}

/**
 * copy constructor of the lattice object
 * @param lat_copy Lattice object to be copied
 */
Lattice::Lattice(const Lattice &lat_copy){

   this->Lx = lat_copy.gLx();
   this->Ly = lat_copy.gLy();
   this->d = lat_copy.gd();

   this->rc2i = lat_copy.rc2i;
   this->i2rc = lat_copy.i2rc;

}

/**
 * empty destructor
 */
Lattice::~Lattice(){ }

/**
 * @return the x dimension
 */
int Lattice::gLx() const {

   return Lx;

}

/**
 * @return the y dimension
 */
int Lattice::gLy() const {

   return Ly;

}

/**
 * @return the physical dimension
 */
int Lattice::gd() const {

   return d;

}

/**
 * access to the lists from outside the class
 * @param row index
 * @param col index
 * @return correpsonding lattice index
 */
int Lattice::grc2i(int r,int c) const {

   return rc2i[r][c];

}

/**
 * access to the lists from outside the class
 * @param i lattice index
 * @param option parameter
 * #return row if option == 0, if == 1 return col index 
 */
int Lattice::gi2rc(int i,int option) const {

   return i2rc[i][option];

}

/**
 * reset the parameters of the lattice object
 * @param Lx x dimension, nr of cols
 * @param Ly y dimension, nr of rows
 * @param d physical dimension of the sites
 */
void Lattice::set(int Lx,int Ly,int d){

   this->Lx = Lx;
   this->Ly = Ly;

   this->d = d;

   rc2i.resize(Ly);

   for(int r = 0;r < Ly;++r)
      rc2i[r].resize(Lx);

   vector<int> v(2);

   int i = 0;

   for(int r = 0;r < Ly;++r)
      for(int c = 0;c < Lx;++c){

         v[0] = r;
         v[1] = c;

         i2rc.push_back(v);

         rc2i[r][c] = i;

         ++i;

      }

}
