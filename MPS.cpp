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

/** 
 * empty constructor: just sets the length of the vector
 */
MPS::MPS() : vector< TArray<complex<double>,3> >( Lx ) { }

/** 
 * standard constructor: just takes in
 * @param L_in length of the chain
 * @param D_in virtual max bond dimension
 * allocates the tensors and fills them randomly
 */
MPS::MPS(int D_in) : vector< TArray<complex<double>,3> >( Lx ) {

   this->D = D_in;

   (*this)[0].resize(1,D,D);

   for(int c = 1;c < Lx-1;++c)
      (*this)[c].resize(D,D,D);

   (*this)[Lx-1].resize(D,D,1);

}

/**
 * copy constructor
 */
MPS::MPS(const MPS &mps_copy) : vector< TArray<complex<double>,3> >(mps_copy) {

   this->D = mps_copy.gD();

}

/**
 * empty destructor
 */
MPS::~MPS(){ }

/**
 * fill a MPS object, by creating a single layer from contracting a peps with a physical vector
 * @param option construct mps from bottom layer (r == 0) if option == 'b', from top layer (r = Ly-1) if option == 't'
 * from left 'l' c = 0, or right 'r' c = Lx-1
 */
void MPS::fill(char option,const PEPS< complex<double> > &peps,const Walker &walker) {

   complex<double> one(1.0,0.0);
   complex<double> zero(0.0,0.0);

   if(option == 'b'){

      for(int c = 0;c < Lx;++c){

         int dim = (*this)[c].size();

         blas::gemv(CblasRowMajor, CblasTrans, d, dim , one, peps(0,c).data(), dim , walker(0,c).data(), 1, zero, (*this)[c].data(), 1);

      }

   }
   else if(option == 't'){//top

      for(int c = 0;c < Lx;++c){

         int dim = (*this)[c].size();

         blas::gemv(CblasRowMajor, CblasTrans, d, dim , one, peps(Ly-1,c).data(), dim , walker(Ly-1,c).data(), 1, zero, (*this)[c].data(), 1);

      }

   }
   else if(option == 'l'){//left

      TArray<complex<double>,3> tmp;

      for(int r = 0;r < Ly;++r){

         int dim = (*this)[r].size();
         tmp.resize(peps(r,0).shape(2),peps(r,0).shape(3),peps(r,0).shape(4));

         blas::gemv(CblasRowMajor, CblasTrans, d, dim , one, peps(r,0).data(), dim , walker(r,0).data(), 1, zero, tmp.data(), 1);

         Permute(tmp,shape(1,2,0),(*this)[r]);

      }

   }
   else{//finally right

      TArray<complex<double>,3> tmp;

      for(int r = 0;r < Ly;++r){

         int dim = (*this)[r].size();
         tmp.resize(peps(r,Lx-1).shape(1),peps(r,Lx-1).shape(2),peps(r,Lx-1).shape(3));

         blas::gemv(CblasRowMajor, CblasTrans, d, dim , one, peps(r,Lx-1).data(), dim , walker(r,Lx-1).data(), 1, zero, tmp.data(), 1);

         Permute(tmp,shape(2,0,1),(*this)[r]);

      }

   }

}

/**
 * @return virtual dimension of the MPS
 */
int MPS::gD() const {

   return D;

}

/**
 * act with an MPO on this MPS, resulting MPS is returned as *this object
 * @param uplo if == 'U' contract with the upper physical index of the MPO, if == 'L', contract with the lower
 * @param mpo the MPO
 */
void MPS::gemv(char uplo,const MPO &mpo){

   int DO = mpo.gD();

   int L = this->size();

   complex<double> one(1.0,0.0);
   complex<double> zero(0.0,0.0);

   if(uplo == 'U'){

      //first site
      TArray<complex<double>,5> tmp;

      enum {i,j,k,l,m,n,o,p,q};

      //dimensions of the new MPS
      int d_phys = (*this)[0].shape(1);
      int DL = 1;
      int DR = (*this)[0].shape(2) * mpo[0].shape(3);

      Contract(one,mpo[0],shape(i,j,k,l),(*this)[0],shape(m,j,n),zero,tmp,shape(m,i,k,n,l));

      (*this)[0] = tmp.reshape_clear(shape(DL,d_phys,DR));

      //middle sites
      for(int c = 1;c < L - 1;++c){

         DL = DR;
         DR = (*this)[c].shape(2) * mpo[c].shape(3);

         Contract(one,mpo[c],shape(i,j,k,l),(*this)[c],shape(m,j,n),zero,tmp,shape(m,i,k,n,l));

         (*this)[c] = tmp.reshape_clear(shape(DL,d_phys,DR));

      }

      DL = DR;
      DR = 1;

      Contract(one,mpo[L - 1],shape(i,j,k,l),(*this)[L - 1],shape(m,j,n),zero,tmp,shape(m,i,k,n,l));

      (*this)[L - 1] = tmp.reshape_clear(shape(DL,d_phys,DR));

   }
   else{//L

      //first site
      TArray<complex<double>,5> tmp;

      enum {i,j,k,l,m,n,o,p,q};

      //dimensions of the new MPS
      int d_phys = (*this)[0].shape(1);
      int DL = 1;
      int DR = (*this)[0].shape(2) * mpo[0].shape(3);

      Contract(one,mpo[0],shape(i,j,k,l),(*this)[0],shape(m,k,n),zero,tmp,shape(m,i,j,n,l));

      (*this)[0] = tmp.reshape_clear(shape(DL,d_phys,DR));

      //middle sites
      for(int c = 1;c < L - 1;++c){

         DL = DR;
         DR = (*this)[c].shape(2) * mpo[c].shape(3);

         Contract(one,mpo[c],shape(i,j,k,l),(*this)[c],shape(m,k,n),zero,tmp,shape(m,i,j,n,l));

         (*this)[c] = tmp.reshape_clear(shape(DL,d_phys,DR));

      }

      DL = DR;
      DR = 1;

      Contract(one,mpo[L - 1],shape(i,j,k,l),(*this)[L - 1],shape(m,k,n),zero,tmp,shape(m,i,j,n,l));

      (*this)[L - 1] = tmp.reshape_clear(shape(DL,d_phys,DR));

   }

   //vdim is increased
   D *= DO;

}

/**
 * canonicalize the mps
 * @param dir Left or Right canonicalization
 * @param norm if true: normalize, else not
 */
void MPS::canonicalize(const BTAS_SIDE &dir,bool norm){

   int length = this->size();

   complex<double> one(1.0,0.0);
   complex<double> zero(0.0,0.0);

   if(dir == Left){//QR

      TArray<complex<double>,2> R;
      TArray<complex<double>,3> tmp;

      for(int i = 0;i < length - 1;++i){

         R.clear();

         //do QR
         Geqrf((*this)[i],R);

         //paste to next matrix
         tmp.clear();

         Contract(one,R,shape(1),(*this)[i + 1],shape(0),zero,tmp);

         (*this)[i + 1] = std::move(tmp);

      }

      if(norm){

         complex<double> nrm = sqrt(Dotc((*this)[length-1],(*this)[length-1]));
         Scal(1.0/nrm,(*this)[length-1]);

      }

   }
   else{//LQ

      TArray<complex<double>,2> L;
      TArray<complex<double>,3> tmp;

      for(int i = length - 1;i > 0;--i){

         L.clear();

         //do QR
         Gelqf(L,(*this)[i]);

         //paste to previous matrix
         tmp.clear();

         Contract(one,(*this)[i - 1],shape(2),L,shape(0),zero,tmp);

         (*this)[i - 1] = std::move(tmp);

      }

      if(norm){

         complex<double> nrm = sqrt(Dotc((*this)[0],(*this)[0]));
         Scal(1.0/nrm,(*this)[0]);

      }

   }

}

/**
 * find an approximate form of the state 'mps' compressed to a bond dimension 'Dc' by performing an SVD on an non-canonical state.
 * @param dir Left or Right - going compression
 * @param Dc the compressed dimension
 * @param mps state to be compressed
 */
void MPS::guess(const BTAS_SIDE &dir,int Dc,const MPS &mps){

   int L = mps.size();

   complex<double> one(1.0,0.0);
   complex<double> zero(0.0,0.0);

   if(dir == Left){

      TArray<complex<double>,3> U;
      TArray<complex<double>,2> V;
      TArray<double,1> S;

      Gesvd('S','S',mps[0],S,U,V,Dc);

      (*this)[0] = std::move(U);

      //multiply S to V
      Dimm(S,V);

      //and contract V with mps on next site
      (*this)[1].clear();

      Contract(one,V,shape(1),mps[1],shape(0),zero,(*this)[1]);

      for(int i = 1;i < L - 1;++i){

         Gesvd('S','S',(*this)[i],S,U,V,Dc);

         (*this)[i] = std::move(U);

         //multiply S to V
         Dimm(S,V);

         //and contract V with mps on next site
         (*this)[i + 1].clear();

         Contract(one,V,shape(1),mps[i + 1],shape(0),zero,(*this)[i + 1]);

      }

   }
   else{

      TArray<complex<double>,2> U;
      TArray<complex<double>,3> V;
      TArray<double,1> S;

      Gesvd('S','S',mps[L - 1],S,U,V,Dc);

      (*this)[L - 1] = std::move(V);

      //multiply U and S
      Dimm(U,S);

      //and contract U with mps on previous site
      (*this)[L - 2].clear();

      Contract(one,mps[L - 2],shape(2),U,shape(0),zero,(*this)[L - 2]);

      for(int i = L - 2;i > 0;--i){

         Gesvd('S','S',(*this)[i],S,U,V,Dc);

         (*this)[i] = std::move(V);

         //multiply S to V
         Dimm(U,S);

         //and contract V with mps on next site
         (*this)[i - 1].clear();

         Contract(one,mps[i - 1],shape(2),U,shape(0),zero,(*this)[i - 1]);

      }

   }

   if(Dc < mps.gD())
      this->D = Dc;
   else
      Dc = mps.gD();

}

/**
 * find the best compression of the state 'mps' a bond dimension 'Dc' by optimizing the tensor in a sweeping fashion
 * @param Dc the compressed dimension
 * @param mps state to be compressed
 */
void MPS::compress(int Dc,const MPS &mps,int n_iter){

   int L = mps.size();

   complex<double> one(1.0,0.0);
   complex<double> zero(0.0,0.0);

   //initial guess by performing svd compression of uncanonicalized state: output is right-canonicalized state
   guess(Right,Dc,mps);

   //construct renormalized operators
   std::vector< TArray<complex<double>,2> > RO(L - 1);
   std::vector< TArray<complex<double>,2> > LO(L - 1);

   compress::init_ro(Right,RO,mps,*this);

   int iter = 0;

   while(iter < n_iter){

      //first site
      (*this)[0].clear();

      Contract(one,mps[0],shape(2),RO[0],shape(1),zero,(*this)[0]);

      //QR
      Geqrf((*this)[0],RO[0]);

      //paste to next matrix
      TArray<complex<double>,3> tmp;

      Contract(one,RO[0],shape(1),(*this)[1],shape(0),zero,tmp);

      (*this)[1] = std::move(tmp);

      compress::update_L(0,LO,mps,*this);

      for(int i = 1;i < L - 1;++i){

         TArray<complex<double>,3> I;

         Contract(one,mps[i],shape(2),RO[i],shape(1),zero,I);

         (*this)[i].clear();

         Contract(one,LO[i - 1],shape(1),I,shape(0),zero,(*this)[i]);

         Geqrf((*this)[i],RO[i]);

         //paste to next matrix
         tmp.clear();

         Contract(one,RO[i],shape(1),(*this)[i + 1],shape(0),zero,tmp);

         (*this)[i + 1] = std::move(tmp);

         compress::update_L(i,LO,mps,*this);

      }

      //and backward!
      (*this)[L - 1].clear();

      Contract(one,LO[L - 2],shape(1),mps[L - 1],shape(0),zero,(*this)[L - 1]);

      //LQ
      Gelqf(LO[L - 2],(*this)[L - 1]);

      //paste to next matrix
      tmp.clear();

      Contract(one,(*this)[L - 2],shape(2),LO[L -  2],shape(0),zero,tmp);

      (*this)[L - 2] = std::move(tmp);

      compress::update_R(L-1,RO,mps,*this);

      for(int i = L - 2;i > 0;--i){

         TArray<complex<double>,3> I;

         Contract(one,mps[i],shape(2),RO[i],shape(1),zero,I);

         (*this)[i].clear();

         Contract(one,LO[i - 1],shape(1),I,shape(0),zero,(*this)[i]);

         Gelqf(LO[i],(*this)[i]);

         //paste to previous matrix
         tmp.clear();

         Contract(one,(*this)[i - 1],shape(2),LO[i],shape(0),zero,tmp);

         (*this)[i - 1] = std::move(tmp);

         compress::update_R(i,RO,mps,*this);

      }

      ++iter;

   }

   this->D = Dc;

}

/**
 * @param bra the bra of the inner product
 * @return the inner product of two MPS's, with *this being the ket
 */
complex<double> MPS::dot(const MPS &bra) const {

   complex<double> one(1.0,0.0);
   complex<double> zero(0.0,0.0);

   TArray<complex<double>,2> E;

   Contract(one,bra[0],shape(0,1),(*this)[0],shape(0,1),zero,E);

   TArray<complex<double>,3> I;

   for(int i = 1;i < this->size();++i){

      I.clear();

      Contract(one,bra[i],shape(0),E,shape(0),zero,I);

      E.clear();

      Contract(one,I,shape(2,0),(*this)[i],shape(0,1),zero,E);

   }

   return E(0,0);

}

/**
 * reduce the dimension of the edge states after MPO action using thin svd.
 */
void MPS::cut_edges() {

   int L = this->size();

   //Left
   TArray<complex<double>,3> U;
   TArray<complex<double>,2> V;
   TArray<double,1> S;

   int i = 0;

   complex<double> one(1.0,0.0);
   complex<double> zero(0.0,0.0);

   //easy compression
   while( (*this)[i].shape(0)*(*this)[i].shape(1) < (*this)[i].shape(2) ){

      U.clear();
      S.clear();
      V.clear();

      Gesvd('S','S',(*this)[i],S,U,V);

      (*this)[i] = std::move(U);

      //multiply S to V
      Dimm(S,V);

      U.clear();

      //and contract V with mps on next site
      Contract(one,V,shape(1),(*this)[i+1],shape(0),zero,U);

      (*this)[i+1] = std::move(U);

      ++i;

   }

   i = L - 1;

   while( (*this)[i].shape(0) > (*this)[L - 1].shape(1)*(*this)[i].shape(2) ){

      //Right
      U.clear();
      V.clear();
      S.clear();

      Gesvd('S','S',(*this)[i],S,V,U);

      (*this)[i] = std::move(U);

      //multiply U and S
      Dimm(V,S);

      //and contract U with mps on previous site
      U.clear();

      Contract(one,(*this)[i-1],shape(2),V,shape(0),zero,U);
      (*this)[i-1] = std::move(U);

      --i;

   }

}
