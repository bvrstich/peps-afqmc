#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>

using std::cout;
using std::endl;
using std::ostream;
using std::complex;
using std::ofstream;

#include "include.h"

using namespace btas;

namespace compress {

   /**
    * initialize the renormalized operator
    * @param dir left or right renormalized operator
    * @param ro input empty, output will contain left or right ro
    * @param bra MPS containing bra state: this is the large-D state
    * @param ket MPS containing ket state: this is the state which is optimized
    */
   void init_ro(const BTAS_SIDE &dir,std::vector< TArray<complex<double>,2> > &ro,const MPS &bra,const MPS &ket){

      int L = bra.size();

      complex<double> one(1.0,0.0);
      complex<double> zero(0.0,0.0);

      ro.resize(L - 1);

      //for Left to right bra is on top!
      if(dir == Left){

         int M = bra[0].shape(2);
         int N = ket[0].shape(2);
         int K = ket[0].shape(0) * ket[0].shape(1);

         ro[0].resize(shape(M,N));

         blas::gemm(CblasRowMajor,CblasConjTrans,CblasNoTrans, M, N, K, one, bra[0].data(),M,ket[0].data(),N,zero,ro[0].data(),N);

         TArray<complex<double>,3> I;

         for(int i = 1;i < L - 1;++i){

            M = ro[i-1].shape(0);
            N = ket[i].shape(1) * ket[i].shape(2);
            K = ro[i-1].shape(1);

            I.resize(shape(M,ket[i].shape(1),ket[i].shape(2)));

            blas::gemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, M, N, K, one, ro[i-1].data(),K,ket[i].data(),N,zero,I.data(),N);

            M = bra[i].shape(2);
            N = ket[i].shape(2);
            K = bra[i].shape(0) * bra[i].shape(1);

            ro[i].resize(shape(M,N));

            blas::gemm(CblasRowMajor,CblasConjTrans,CblasNoTrans, M, N, K, one, bra[i].data(),M,I.data(),N,zero,ro[i].data(),N);

         }

      }
      else{

         int M = ket[L-1].shape(0);
         int N = bra[L-1].shape(0);
         int K = ket[L-1].shape(1) * ket[L-1].shape(2);

         ro[L-2].resize(shape(M,N));

         blas::gemm(CblasRowMajor,CblasNoTrans,CblasConjTrans, M, N, K, one, ket[L-1].data(),K,bra[L-1].data(),K,zero,ro[L-2].data(),N);

         TArray<complex<double>,3> I;

         for(int i = L - 2;i > 0;--i){

            M = ket[i].shape(0) * ket[i].shape(1);
            N = ro[i].shape(1);
            K = ket[i].shape(0);

            I.resize(shape(ket[i].shape(0),ket[i].shape(1),N));

            blas::gemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, M, N, K, one, ket[i].data(),K,ro[i].data(),N,zero,I.data(),N);

            M = ket[i].shape(0) ;
            N = bra[i].shape(0);
            K = bra[i].shape(1)*bra[i].shape(2);

            ro[i-1].resize(shape(M,N));

            blas::gemm(CblasRowMajor,CblasNoTrans,CblasConjTrans, M, N, K, one, I.data(),K,bra[i].data(),K,zero,ro[i-1].data(),N);

         }

      }

   }

   /**
    * update the left renormalized operator
    * @param site index of the site
    * @param LO left renormalized operator
    * @param bra MPS containing bra state: this is the large-D state
    * @param ket MPS containing ket state: this is the state which is optimized
    */
   void update_L(int site,std::vector< TArray<complex<double>,2> > &LO,const MPS &bra,const MPS &ket){

      complex<double> one(1.0,0.0);
      complex<double> zero(0.0,0.0);

      LO[site].clear();

      if(site == 0)
         Contract(one,ket[0],shape(0,1),bra[0],shape(0,1),zero,LO[0]);
      else{

         TArray<complex<double>,3> I;

         Contract(one,ket[site],shape(0),LO[site - 1],shape(0),zero,I);

         Contract(one,I,shape(2,0),bra[site],shape(0,1),zero,LO[site]);

      }

   }

   /**
    * update the right renormalized operator
    * @param site index of the site
    * @param RO Right renormalized operator
    * @param bra MPS containing bra state: this is the large-D state
    * @param ket MPS containing ket state: this is the state which is optimized
    */
   void update_R(int site,std::vector< TArray<complex<double>,2> > &RO,const MPS &bra,const MPS &ket){

      complex<double> one(1.0,0.0);
      complex<double> zero(0.0,0.0);

      RO[site - 1].clear();

      int L = bra.size();

      if(site == L - 1)
         Contract(one,ket[L - 1],shape(1,2),bra[L - 1],shape(1,2),zero,RO[L - 2]);
      else{

         TArray<complex<double>,3> I;

         Contract(one,ket[site],shape(2),RO[site],shape(0),zero,I);

         Contract(one,I,shape(1,2),bra[site],shape(1,2),zero,RO[site-1]);

      }

   }

}
