#ifndef WALKER_H
#define WALKER_H

#include <iostream>
#include <iomanip>
#include <complex>
#include <vector>

using std::complex;
using std::vector;

using namespace btas;

#include "PEPS.h"

class Propagator;

/**
 * class definition of Walker, made to describe the product state walkers. An array of L size-2 vector. Each site represents a rotation of the spin.
 */
class Walker : public vector< TArray<complex<double>,1> > {

   public:

      //empty contstructor
      Walker();
   
      //Constructor copying an entire Walker
      Walker(const Walker &walker);
      
      //Destructor
      virtual ~Walker();

      double gWeight() const;

      void sWeight(double);

      void multWeight(double);

      void calc_properties(char,const PEPS< complex<double> > &);
   
      complex<double> gOverlap() const;

      complex<double> gEL() const;

      void normalize();

      const std::vector< complex<double> > &gVL() const;

      complex<double> gVL(int,int) const;

      const TArray<complex<double>,1> &operator()(int r,int c) const;

      void propagate(const Propagator &);

      void copy_essential(const Walker &);

  private:

      //!The walker weight
      double weight;

      //!The walker overlap with the trial wfn
      complex<double> overlap;

      //!local energy
      complex<double> EL;

      //!local auxiliary operators: <PsiT|v|phi>/<PsiT|phi>
      std::vector< complex<double> > VL;

};

#endif
