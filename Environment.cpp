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

MPO Environment::mpo;

int Environment::D_aux;

SL_PEPS Environment::U;

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

   mpo = MPO(D);

   U = SL_PEPS(D);

}

/**
 * construct the enviroment mps's for the input PEPS: make sure the appropriate SL_PEPS's are filled
 * @param option if 'H' construct top and bottom environment
 *               if 'V' construct left and right environment
 * @param peps input PEPS<double>
 */
void Environment::calc_env(char option,const PEPS< complex<double> > &peps,const Walker &walker){

   if(option == 'H'){

      U.fill('H',peps,walker);

      //construct bottom layer
      b[0].fill(0,U);

      for(int r = 1;r < Ly - 1;++r){

         //i'th row as MPO
         mpo.fill('H',r,peps,walker);

         MPS tmp(b[r - 1]);

         //apply to form MPS with bond dimension D^2
         tmp.gemv('L',mpo);

         //reduce the dimensions of the edge states using thin svd
         tmp.cut_edges();

         //compress in sweeping fashion
         b[r].compress(D_aux,tmp,1);

      }

      //then construct top layer
      t[Ly - 2].fill('t',peps,walker);

      for(int i = Ly - 2;i > 0;--i){

         //i'th row as MPO
         mpo.fill('H',i,peps,walker);

         //apply to form MPS with bond dimension D^4
         MPS tmp(t[i]);

         tmp.gemv('U',mpo);

         //reduce the dimensions of the edge states using thin svd
         tmp.cut_edges();

         //compress in sweeping fashion
         t[i - 1].compress(D_aux,tmp,1);

      }

   }
   else{//Vertical

      //then left layer
      l[0].fill('l',peps,walker);

      for(int i = 1;i < Lx - 1;++i){

         //i'th col as MPO
         mpo.fill('V',i,peps,walker);

         MPS tmp(l[i - 1]);

         //apply to form MPS with bond dimension D^4
         tmp.gemv('L',mpo);

         //reduce the dimensions of the edge states using thin svd
         tmp.cut_edges();

         //compress in sweeping fashion
         l[i].compress(D_aux,tmp,1);

      }

      //finally construct right layer
      r[Lx - 2].fill('r',peps,walker);

      for(int i = Lx - 2;i > 0;--i){

         //i'th row as MPO
         mpo.fill('V',i,peps,walker);

         //apply to form MPS with bond dimension D^4
         MPS tmp(r[i]);

         tmp.gemv('U',mpo);

         //reduce the dimensions of the edge states using thin svd
         tmp.cut_edges();

         //compress in sweeping fashion
         r[i - 1].compress(D_aux,tmp,1);

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

/**
 * construct the enviroment for a specific row/column of the input PEPS, it is assumed that 
 * all prerequisites for the construction of the environment are there!
 * @param option if 'L' construct full left environment
 *               if 'R' construct full right environment
 *               if 'T' construct full top environment
 *               if 'B' construct full bottom environment
 * @param rc row or column index for the L,R,T or B environment to be constructed
 * @param peps input PEPS<double>
 * @param D_aux dimension to which environment will be compressed
 */
void Environment::calc_env(char option,int rc,const PEPS< complex<double> > &peps,const Walker &walker){

   if(option == 'B'){

      //construct bottom layer
      if(rc == 0)
         b[0].fill('b',peps,walker);
      else{

         //i'th row as MPO
         mpo.fill('H',rc,peps,walker);

         MPS tmp(b[rc - 1]);

         //apply to form MPS with bond dimension D^2
         tmp.gemv('L',mpo);

         //reduce the dimensions of the edge states using thin svd
         tmp.cut_edges();

         //compress in sweeping fashion
         b[rc].compress(D_aux,tmp,1);

      }

   }
   else if(option == 'T'){

      //then construct top layer
      if(rc == Ly-1)
         t[Ly - 2].fill('t',peps,walker);
      else{

         //i'th row as MPO
         mpo.fill('H',rc,peps,walker);

         //apply to form MPS with bond dimension D^4
         MPS tmp(t[rc]);

         tmp.gemv('U',mpo);

         //reduce the dimensions of the edge states using thin svd
         tmp.cut_edges();

         //compress in sweeping fashion
         t[rc - 1].compress(D_aux,tmp,1);

      }

   }
   else if(option == 'L'){

      //then left layer
      if(rc == 0)
         l[0].fill('l',peps,walker);
      else{

         //rc'th col as MPO
         mpo.fill('V',rc,peps,walker);

         MPS tmp(l[rc - 1]);

         //apply to form MPS with bond dimension D^4
         tmp.gemv('L',mpo);

         //reduce the dimensions of the edge states using thin svd
         tmp.cut_edges();

         //compress in sweeping fashion
         l[rc].compress(D_aux,tmp,1);

      }

   }
   else{//option == R

      //finally construct right layer
      if(rc == Lx - 1)
         r[Lx - 2].fill('r',peps,walker);
      else{

         //i'th row as MPO
         mpo.fill('V',rc,peps,walker);

         //apply to form MPS with bond dimension D^4
         MPS tmp(r[rc]);

         tmp.gemv('U',mpo);

         //reduce the dimensions of the edge states using thin svd
         tmp.cut_edges();

         //compress in sweeping fashion
         r[rc - 1].compress(D_aux,tmp,1);

      }

   }

}
