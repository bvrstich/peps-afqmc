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

//statics
vector< MPS > Environment::l;
vector< MPS > Environment::r;
vector< MPS > Environment::t;
vector< MPS > Environment::b;

MPO Environment::mpo;

int Environment::D_aux;

/** 
 * initialize all the static variables
 * @param D bond dimension of the trial peps
 * @param D_aux_in auxiliary bond dimension for the contractions
 */
void Environment::init(int D,int D_aux_in){

   D_aux = D_aux_in;

   int Lx = Global::Lx;
   int Ly = Global::Ly;

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

}

/**
 * construct the enviroment mps's for the input PEPS
 * @param option if 'L' construct full left environment
 *               if 'R' construct full right environment
 *               if 'T' construct full top environment
 *               if 'B' construct full bottom environment
 *               if 'A' construct all environments
 * @param peps input PEPS<double>
 * @param D_aux dimension to which environment will be compressed
 */
void Environment::calc_env(char option,const PEPS< complex<double> > &peps,const Walker &walker){

   int Lx = Global::Lx;
   int Ly = Global::Ly;

   int d = Global::d;

   if(option == 'B' || option == 'A'){

      //construct bottom layer
      b[0].fill('b',peps,walker);

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

   }

   if(option == 'T' || option == 'A'){

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

   if(option == 'L' || option == 'A'){

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

   }

   if(option == 'R' || option == 'A'){

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

   int Lx = Global::Lx;
   int Ly = Global::Ly;

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
 /*
void Environment::calc_env(char option,int rc,const PEPS<double> &peps,int D_aux){

   int Lx = Global::lat.gLx();
   int Ly = Global::lat.gLy();

   int d = Global::lat.gd();

   if(option == 'B'){

      //construct bottom layer
      if(rc == 0)
         b[0] = MPS<double>('b',peps,peps);
      else{

         //i'th row as MPO
         MPO<double> mpo('H',rc,peps,peps);

         MPS<double> tmp(b[rc - 1]);

         //apply to form MPS with bond dimension D^4
         tmp.gemv('L',mpo);

         //reduce the dimensions of the edge states using thin svd
         tmp.cut_edges();

         //compress in sweeping fashion
         b[rc].resize(Lx);
         b[rc].compress(D_aux,tmp,5);

      }

   }
   else if(option == 'T'){

      //then construct top layer
      if(rc == Ly-1)
         t[Ly - 2] = MPS<double>('t',peps,peps);
      else{

         //i'th row as MPO
         MPO<double> mpo('H',rc,peps,peps);

         //apply to form MPS with bond dimension D^4
         MPS<double> tmp(t[rc]);

         tmp.gemv('U',mpo);

         //reduce the dimensions of the edge states using thin svd
         tmp.cut_edges();

         //compress in sweeping fashion
         t[rc - 1].resize(Lx);
         t[rc - 1].compress(D_aux,tmp,5);

      }

   }
   else if(option == 'L'){

      //then left layer
      if(rc == 0)
         l[0] = MPS<double>('l',peps,peps);
      else{

         //i'th col as MPO
         MPO<double> mpo('V',rc,peps,peps);

         MPS<double> tmp(l[rc - 1]);

         //apply to form MPS with bond dimension D^4
         tmp.gemv('L',mpo);

         //reduce the dimensions of the edge states using thin svd
         tmp.cut_edges();

         //compress in sweeping fashion
         l[rc].resize(Ly);
         l[rc].compress(D_aux,tmp,5);

      }

   }
   else{//option == R

      //finally construct right layer
      if(rc == Lx - 1)
         r[Lx - 2] = MPS<double>('r',peps,peps);
      else{

         //i'th row as MPO
         MPO<double> mpo('V',rc,peps,peps);

         //apply to form MPS with bond dimension D^4
         MPS<double> tmp(r[rc]);

         tmp.gemv('U',mpo);

         //reduce the dimensions of the edge states using thin svd
         tmp.cut_edges();

         //compress in sweeping fashion
         r[rc - 1].resize(Ly);
         r[rc - 1].compress(D_aux,tmp,5);

      }

   }

}
*/
