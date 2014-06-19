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

//statics
vector< MPS > Environment::l;
vector< MPS > Environment::r;
vector< MPS > Environment::t;
vector< MPS > Environment::b;

int Environment::D_aux;

SL_PEPS Environment::U;
SL_PEPS Environment::Sx;
SL_PEPS Environment::Sy;
SL_PEPS Environment::Sz;

/** 
 * initialize all the static variables
 * @param D bond dimension of the trial peps
 * @param D_aux_in auxiliary bond dimension for the contractions
 */
void Environment::init(int D,int D_aux_in){

   D_aux = D_aux_in;

   t.resize(Ly - 1);
   b.resize(Ly - 1);

   b[0] = MPS(D);
   t[0] = MPS(D_aux);

   for(int i = 1;i < Ly - 2;++i){

      t[i] = MPS(D_aux);
      b[i] = MPS(D_aux);

   }

   b[Ly-2] = MPS(D_aux);
   t[Ly-2] = MPS(D);

   r.resize(Lx - 1);
   l.resize(Lx - 1);

   l[0] = MPS(D);
   r[0] = MPS(D_aux);

   for(int i = 1;i < Lx - 2;++i){

      l[i] = MPS(D_aux);
      r[i] = MPS(D_aux);

   }

   l[Lx-2] = MPS(D_aux);
   r[Lx-2] = MPS(D);

   U = SL_PEPS(D);

   Sx = SL_PEPS(D);
   Sy = SL_PEPS(D);
   Sz = SL_PEPS(D);

}

/**
 * construct the enviroment mps's for the input PEPS: make sure the appropriate SL_PEPS's are filled
 * make sure that the SL_PEPS U is filled before with the right 'H' or 'V' option before calling this function!
 * @param option if 'H' construct top and bottom environment
 *               if 'V' construct left and right environment
 * @param peps input PEPS<double>
 */
void Environment::calc_env(char option,const PEPS< complex<double> > &peps,const Walker &walker){

   if(option == 'H'){

      //construct bottom layer
      b[0].fill('b',U);

      for(int r = 1;r < Ly - 1;++r){

         MPS tmp(b[r - 1]);

         //apply to form MPS with bond dimension D^2
         tmp.gemv('L','H',r,U);

         //reduce the dimensions of the edge states using thin svd
         tmp.cut_edges();

         //compress in sweeping fashion
         b[r].compress(D_aux,tmp,1);

      }

      //then construct top layer
      t[Ly - 2].fill('t',U);

      for(int r = Ly - 2;r > 0;--r){

         //apply to form MPS with bond dimension D^4
         MPS tmp(t[r]);

         tmp.gemv('U','H',r,U);

         //reduce the dimensions of the edge states using thin svd
         tmp.cut_edges();

         //compress in sweeping fashion
         t[r - 1].compress(D_aux,tmp,1);

      }

   }
   else{//Vertical

      //then left layer
      l[0].fill('l',U);

      for(int c = 1;c < Lx - 1;++c){

         //i'th col as MPO
         MPS tmp(l[c - 1]);

         //apply to form MPS with bond dimension D^4
         tmp.gemv('L','V',c,U);

         //reduce the dimensions of the edge states using thin svd
         tmp.cut_edges();

         //compress in sweeping fashion
         l[c].compress(D_aux,tmp,1);

      }

      //finally construct right layer
      r[Lx - 2].fill('r',U);

      for(int c = Lx - 2;c > 0;--c){

         //apply to form MPS with bond dimension D^4
         MPS tmp(r[c]);

         tmp.gemv('U','V',c,U);

         //reduce the dimensions of the edge states using thin svd
         tmp.cut_edges();

         //compress in sweeping fashion
         r[c - 1].compress(D_aux,tmp,1);

      }

   }

}

/**
 * test if the enviroment is correctly contracted
 */
void Environment::test_env(){

   cout << endl;
   cout << "FROM BOTTOM TO TOP" << endl;
   cout << endl;
   for(int i = 0;i < Ly - 1;++i)
      cout << i << "\t" << b[i].dot(t[i]) << endl;

   cout << endl;
   cout << "FROM LEFT TO RIGHT" << endl;
   cout << endl;
   for(int i = 0;i < Lx - 1;++i)
      cout << i << "\t" << r[i].dot(l[i]) << endl;
   cout << endl;

}
