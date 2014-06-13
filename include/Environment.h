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
      static vector< MPS > l;
      static vector< MPS > r;
      static vector< MPS > t;
      static vector< MPS > b;

      //!for the contraction
      static MPO mpo;

      //!auxiliary bond dimension for the environment contractions
      static int D_aux;

   private:

};

#endif
