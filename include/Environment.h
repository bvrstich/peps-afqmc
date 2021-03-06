#ifndef ENVIRONMENT_H
#define ENVIRONMENT_H

#include <iostream>
#include <fstream>
#include <vector>

using std::ostream;
using std::vector;

/**
 * @author Brecht Verstichel
 * @data 02-05-2014\n\n
 * Class used to calculate the enviroment of a peps. Needed for the calculation of expectation values and the update of tensors.
 */
class Environment {

   public:

      static void init(int,int);

      static void calc_env(char,const PEPS< complex<double> > &,const Walker &walker);

      static void calc_env(char,int,const PEPS< complex<double> > &,const Walker &walker);

      static void test_env();

      //!stores an array environment MPS's for l(eft) , r(ight), t(op) and b(ottom)
      static vector< vector< MPS > > l;
      static vector< vector< MPS > > r;
      static vector< vector< MPS > > t;
      static vector< vector< MPS > > b;

      //!overlap between walker and trial PEPS
      static vector< SL_PEPS > U;

      //!expectation between walker and trial, with operator O=Sn
      static vector< SL_PEPS > Sx;
      static vector< SL_PEPS > Sy;
      static vector< SL_PEPS > Sz;

      //!auxiliary bond dimension for the environment contractions
      static int D_aux;

   private:

};

#endif
