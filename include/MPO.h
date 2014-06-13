#ifndef MPO_H
#define MPO_H

#include <iostream>
#include <fstream>
#include <vector>

using std::ostream;
using std::vector;

using namespace btas;

#include "PEPS.h"

/**
 * @author Brecht Verstichel
 * @date 26-03-2014\n\n
 * This class MPO is a class written for the construction of matrix product operators without symmetry.
 * More specifically it will be used for the contraction of PEPS networks. Where the reduction to a MPO-like form is done.
 */
class MPO : public vector< TArray<complex<double>,4> > {

   public:

      //constructors
      MPO();

      MPO(int);

      //copy constructor
      MPO(const MPO &);

      //destructor
      virtual ~MPO();

      void fill(char,int,const PEPS< complex<double> > &,const Walker &);

      int gD() const;

   private:

      //!dimension of the bonds
      int D;


};

#endif

/* vim: set ts=3 sw=3 expandtab :*/
