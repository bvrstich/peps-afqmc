#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <complex>
#include <vector>

#include "include.h"

using std::cout;
using std::ifstream;
using std::ofstream;
using std::endl;

/**
 * empty constructor:
 */
Walker::Walker() : std::vector< TArray<complex<double>,1> >( Global::lat.gLx() * Global::lat.gLy() ){ }

/**
 * construct a Walker object: initialize on AF state
 * @param n_trot_in number of trotter terms
 */
Walker::Walker(int n_trot_in) : std::vector< TArray<complex<double>,1> >( Global::lat.gLx() * Global::lat.gLy() ){

   this->n_trot = n_trot_in;

   int Lx = Global::lat.gLx();
   int Ly = Global::lat.gLy();
   int d = Global::lat.gd();

   for(int r = 0;r < Ly;++r)
      for(int c = 0;c < Lx;++c){

         (*this)[ Global::lat.grc2i(r,c) ].resize(d);

         if( (r + c)%2 == 0){

            (*this)[ Global::lat.grc2i(r,c) ](0) = 0.0;
            (*this)[ Global::lat.grc2i(r,c) ](1) = 1.0;

         }
         else{

            (*this)[ Global::lat.grc2i(r,c) ](0) = 1.0;
            (*this)[ Global::lat.grc2i(r,c) ](1) = 0.0;

         }

      }

   w_xyz.resize(3*Lx*Ly);

   for(int i = 0;i < Lx*Ly;++i)
      for(int r = 0;r < 3;++r)//regular vector
         w_xyz[3*i + r].resize(d);

   VL.resize(3*n_trot);

}

/**
 * copy constructor
 * @param walker input Walker object to be copied
 */
Walker::Walker(const Walker &walker) : std::vector< TArray<complex<double>,1> >(walker) {

   this->weight = walker.gWeight();
   this->n_trot = walker.gn_trot();
   this->overlap = walker.gOverlap();
   this->EL = walker.gEL();
   this->VL = walker.gVL();
   this->w_xyz = walker.w_xyz;

}

/**
 * destructor
 */
Walker::~Walker(){ }

/** 
 * @return the weight corresponding to the walker
 */
double Walker::gWeight() const{

   return weight; 

}

/**
 * muliply the weight by a factor
 */
void Walker::multWeight(double factor){

   weight *= factor; 

}

/**
 * set new weight
 */
void Walker::sWeight(double new_weight){

   weight = new_weight;

}


/** 
 * @return the number of trotter terms
 */
int Walker::gn_trot() const {

   return n_trot;

}

/** 
 * @return the overlap of the walker with the Trial
 */
complex<double> Walker::gOverlap() const{

   return overlap; 

}

/** 
 * @return the local energy
 */
complex<double> Walker::gEL() const{

   return EL; 

}

/** 
 * @return the shifts, i.e. the vector containing the projected expectation values of the auxiliary field operators
 */
const std::vector< complex<double> > &Walker::gVL() const{

   return VL; 

}

/**
 * access to the individual tensors: const version
 * @param r row index
 * @param c col index
 * @return the tensor on site (r,c)
 */
const TArray<complex<double>,1> &Walker::operator()(int r,int c) const {

   return (*this)[Global::lat.grc2i(r,c)];

}
