#ifndef SL_PEPS_H
#define SL_PEPS_H

#include <iostream>
#include <fstream>
#include <vector>

using std::ostream;
using std::vector;

using namespace btas;

/**
 * @author Brecht Verstichel
 * @date 16-06-2014\n\n
 * This class SL_PEPS is a class written for the construction of the overlap between a PEPS and a D=1 walker
 * SL stands for single layer, since the physical dimension is contracted over with the walker only a single layer remains
 */
class SL_PEPS : public vector< TArray<complex<double>,4> > {

   public:

      SL_PEPS();

      SL_PEPS(int D);

      //copy constructor
      SL_PEPS(const SL_PEPS &);

      //destructor
      virtual ~SL_PEPS();

      int gD() const;

      //fill by just overlap
      void fill(char,const PEPS< complex<double> > &,const Walker &);

      //fill with expectation value between walker, peps and operator
      void fill(char,const PEPS< complex<double> > &,TArray<complex<double>,2> &,const Walker &);

      const TArray<complex<double>,4> &operator()(int,int) const;

      TArray<complex<double>,4> &operator()(int,int);



   private:

      //!dimension of the bonds
      int D;

};

#endif

/* vim: set ts=3 sw=3 expandtab :*/
