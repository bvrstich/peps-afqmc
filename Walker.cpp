#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <complex>
#include <vector>

using std::cout;
using std::ifstream;
using std::ofstream;
using std::endl;

#include "include.h"

using namespace global;

/**
 * construct a Walker object: initialize on AF state
 * @param n_trot_in number of trotter terms
 */
Walker::Walker() : std::vector< TArray<complex<double>,1> >( Lx * Ly ){

   weight = 1.0;

   for(int r = 0;r < Ly;++r)
      for(int c = 0;c < Lx;++c){

         (*this)[ r*Lx + c ].resize(d);

         if( (r + c)%2 == 0){

            (*this)[ r*Lx + c ](0) = 0.0;
            (*this)[ r*Lx + c ](1) = 1.0;

         }
         else{

            (*this)[ r*Lx + c ](0) = 1.0;
            (*this)[ r*Lx + c ](1) = 0.0;

         }

      }

   VL.resize(3*Trotter::n_trot);

}

/**
 * copy constructor
 * @param walker input Walker object to be copied
 */
Walker::Walker(const Walker &walker) : std::vector< TArray<complex<double>,1> >(walker) {

   this->weight = walker.gWeight();
   this->overlap = walker.gOverlap();
   this->EL = walker.gEL();
   this->VL = walker.gVL();

}

/**
 * destructor
 */
Walker::~Walker(){ }

/** 
 * @return the weight corresponding to the walker
 */
double Walker::gWeight() const{

   return weight; 

}

/**
 * muliply the weight by a factor
 */
void Walker::multWeight(double factor){

   weight *= factor; 

}

/**
 * set new weight
 */
void Walker::sWeight(double new_weight){

   weight = new_weight;

}


/** 
 * @return the overlap of the walker with the Trial
 */
complex<double> Walker::gOverlap() const{

   return overlap; 

}

/** 
 * @return the local energy
 */
complex<double> Walker::gEL() const{

   return EL; 

}

/** 
 * @return the shifts, i.e. the vector containing the projected expectation values of the auxiliary field operators
 */
const std::vector< complex<double> > &Walker::gVL() const{

   return VL; 

}

/**
 * access to the auxiliary field expectation values
 * @param k trotter index
 * @param r x,y or z
 * @return VL(k,r)
 */
complex<double> Walker::gVL(int k,int r) const {

   return VL[r*Trotter::n_trot + k];

}

/**
 * access to the individual tensors: const version
 * @param r row index
 * @param c col index
 * @return the tensor on site (r,c)
 */
const TArray<complex<double>,1> &Walker::operator()(int r,int c) const {

   return (*this)[r*Lx + c];

}

/**
 * calculate the overlap, the S(x,y,z) expectation values and the energy in one sweep
 * @param option 'H' = horizontal sweep, 'V' is vertical
 * @param peps trial wave function represented as peps
 */
void Walker::calc_properties(char option,const PEPS< complex<double> >& peps){

#ifdef _OPENMP
      int myID = omp_get_thread_num();
#else
      int myID = 0;
#endif

   // ---- || evaluate the expectation values in an MPO/MPS manner, first from bottom to top, then left to right || ----

   complex<double> energy(0.0,0.0);

   complex<double> one(1.0,0.0);
   complex<double> zero(0.0,0.0);

   if(option == 'H'){

      //calculate the single layer contractions first:
      Environment::U[myID].fill('H',peps,*this);

      Environment::Sx[myID].fill('H',peps,Sx,*this);
      Environment::Sy[myID].fill('H',peps,Sy,*this);
      Environment::Sz[myID].fill('H',peps,Sz,*this);

      //first construct the top and bottom (horizontal) environment layers
      Environment::calc_env('H',peps,*this);

      // #################################################################
      // ### ---- from bottom to top: contract in mps/mpo fashion ---- ### 
      // #################################################################

      // -- (1) -- || bottom row: similar to overlap calculation

      //first construct the right renormalized operators
      vector< ZArray<2> > R(Lx - 1);

      //first the rightmost operator
      ZArray<4> tmp4;
      ZArray<3> tmp3;

      //tmp comes out index (t,b)
      Contract(one,Environment::t[myID][0][Lx - 1],shape(1),Environment::b[myID][0][Lx - 1],shape(1),zero,tmp4);

      //reshape tmp to a 2-index array
      R[Lx - 2] = tmp4.reshape_clear(shape(Environment::t[myID][0][Lx - 1].shape(0),Environment::b[myID][0][Lx - 1].shape(0)));

      //now construct the rest
      for(int col = Lx - 2;col > 0;--col){

         tmp3.clear();
         Contract(one,Environment::t[myID][0][col],shape(2),R[col],shape(0),zero,tmp3);

         Contract(one,tmp3,shape(1,2),Environment::b[myID][0][col],shape(1,2),zero,R[col-1]);

      }

      //4 left going operators: Sx, Sy, Sz, and 1
      ZArray<2> LSx;
      ZArray<2> LSy;
      ZArray<2> LSz;
      ZArray<2> LU;

      TArray<complex<double>,5> tmp5;

      //tmp comes out index (t,b)
      Contract(one,Environment::t[myID][0][0],shape(1),Environment::Sx[myID](0,0),shape(1),zero,tmp5);

      LSx = tmp5.reshape_clear(shape(Environment::t[myID][0][0].shape(2),Environment::Sx[myID](0,0).shape(3)));

      //then Sy
      Contract(one,Environment::t[myID][0][0],shape(1),Environment::Sy[myID](0,0),shape(1),zero,tmp5);

      LSy = tmp5.reshape_clear(shape(Environment::t[myID][0][0].shape(2),Environment::Sy[myID](0,0).shape(3)));

      //then Sz
      Contract(one,Environment::t[myID][0][0],shape(1),Environment::Sz[myID](0,0),shape(1),zero,tmp5);

      LSz = tmp5.reshape_clear(shape(Environment::t[myID][0][0].shape(2),Environment::Sz[myID](0,0).shape(3)));

      //finally unity
      Contract(one,Environment::t[myID][0][0],shape(1),Environment::U[myID](0,0),shape(1),zero,tmp5);

      LU = tmp5.reshape_clear(shape(Environment::t[myID][0][0].shape(2),Environment::Sz[myID](0,0).shape(3)));

      int dim = R[0].size();

      //now contract x,y and z with R for local expectation values:
      auxvec[myID][0][0] = blas::dot(dim,LSx.data(),1,R[0].data(),1);
      auxvec[myID][0][1] = blas::dot(dim,LSy.data(),1,R[0].data(),1);
      auxvec[myID][0][2] = blas::dot(dim,LSz.data(),1,R[0].data(),1);

      //now for the middle terms
      for(int col = 1;col < Lx - 1;++col){

         //first close down the x,y and z terms from the previous site for the energy

         //construct the right intermediate contraction (paste top to right)
         tmp3.clear();
         Contract(one,Environment::t[myID][0][col],shape(2),R[col],shape(0),zero,tmp3);

         // 1) paste Sx to the right
         int M = tmp3.shape(0);
         int N = Environment::Sx[myID](0,col).shape(0);
         int K = tmp3.shape(1) * tmp3.shape(2);

         blas::gemm(CblasRowMajor, CblasNoTrans, CblasTrans, M, N, K, one, tmp3.data(),K,Environment::Sx[myID](0,col).data(),K,zero,R[col-1].data(),N);

         //contract with left Sx
         energy += Dot(LSx,R[col - 1]);

         // 2) then paste Sy to the right
         blas::gemm(CblasRowMajor, CblasNoTrans, CblasTrans, M, N, K, one, tmp3.data(),K,Environment::Sy[myID](0,col).data(),K,zero,R[col-1].data(),N);

         //contract with left Sy
         energy += Dot(LSy,R[col - 1]);

         // 3) finally Sz
         blas::gemm(CblasRowMajor, CblasNoTrans, CblasTrans, M, N, K, one, tmp3.data(),K,Environment::Sz[myID](0,col).data(),K,zero,R[col-1].data(),N);

         //contract with left Sz
         energy += Dot(LSz,R[col - 1]);

         //construct left renormalized operators for next site: first paste top to Left unity
         tmp3.clear();
         Contract(one,LU,shape(0),Environment::t[myID][0][col],shape(0),zero,tmp3);

         // 1) construct new Sx left operator
         LSx.resize(Environment::t[myID][0][col].shape(2),Environment::Sx[myID](0,col).shape(3));

         M = tmp3.shape(2);
         N = Environment::Sx[myID](0,col).shape(3);
         K = tmp3.shape(0) * tmp3.shape(1);

         blas::gemm(CblasRowMajor, CblasTrans, CblasNoTrans, M, N, K, one, tmp3.data(),M,Environment::Sx[myID](0,col).data(),N,zero,LSx.data(),N);

         // 2) construct new Sy left operator
         LSy.resize(Environment::t[myID][0][col].shape(2),Environment::Sy[myID](0,col).shape(3));
         blas::gemm(CblasRowMajor, CblasTrans, CblasNoTrans, M, N, K, one, tmp3.data(),M,Environment::Sy[myID](0,col).data(),N,zero,LSy.data(),N);

         // 3) construct new Sz left operator
         LSz.resize(Environment::t[myID][0][col].shape(2),Environment::Sz[myID](0,col).shape(3));
         blas::gemm(CblasRowMajor, CblasTrans, CblasNoTrans, M, N, K, one, tmp3.data(),M,Environment::Sz[myID](0,col).data(),N,zero,LSz.data(),N);

         //now contract x,y and z with R for local expectation values:
         dim = R[col].size();

         auxvec[myID][col][0] = blas::dot(dim,LSx.data(),1,R[col].data(),1);
         auxvec[myID][col][1] = blas::dot(dim,LSy.data(),1,R[col].data(),1);
         auxvec[myID][col][2] = blas::dot(dim,LSz.data(),1,R[col].data(),1);

         // 4) finally construct new unity on the left
         LU.resize(Environment::t[myID][0][col].shape(2),Environment::U[myID](0,col).shape(3));
         blas::gemm(CblasRowMajor, CblasTrans, CblasNoTrans, M, N, K, one, tmp3.data(),M,Environment::U[myID](0,col).data(),N,zero,LU.data(),N);

      }

      //last site of bottom row:close down the left x,y and z

      //1) Sx to close down LSx
      Contract(one,Environment::t[myID][0][Lx-1],shape(1),Environment::Sx[myID](0,Lx-1),shape(1),zero,tmp5);

      R[Lx-2] = tmp5.reshape_clear(shape(Environment::t[myID][0][Lx-1].shape(0),Environment::Sx[myID](0,Lx-1).shape(0)));

      energy += Dot(LSx,R[Lx-2]);

      dim = R[Lx-2].size();
      auxvec[myID][Lx-1][0] = blas::dot(dim,LU.data(),1,R[Lx-2].data(),1);

      //2) Sy to close down Ly
      Contract(one,Environment::t[myID][0][Lx-1],shape(1),Environment::Sy[myID](0,Lx-1),shape(1),zero,tmp5);

      R[Lx-2] = tmp5.reshape_clear(shape(Environment::t[myID][0][Lx-1].shape(0),Environment::Sy[myID](0,Lx-1).shape(0)));

      energy += Dot(LSy,R[Lx-2]);

      auxvec[myID][Lx-1][1] = blas::dot(dim,LU.data(),1,R[Lx-2].data(),1);

      //3) Sz to close down Lz
      Contract(one,Environment::t[myID][0][Lx-1],shape(1),Environment::Sz[myID](0,Lx-1),shape(1),zero,tmp5);

      R[Lx-2] = tmp5.reshape_clear(shape(Environment::t[myID][0][Lx-1].shape(0),Environment::Sz[myID](0,Lx-1).shape(0)));

      energy += Dot(LSz,R[Lx-2]);

      auxvec[myID][Lx-1][2] = blas::dot(dim,LU.data(),1,R[Lx-2].data(),1);

      // -- (2) -- now move from bottom to top calculating everything like an MPO/MPS expectation value

      //Right renormalized operators
      vector< TArray<complex<double>,3> > RO(Lx - 1);

      //4 left renormalized operators needed
      TArray<complex<double>,3> LOSx;
      TArray<complex<double>,3> LOSy;
      TArray<complex<double>,3> LOSz;
      TArray<complex<double>,3> LOU;

      for(int row = 1;row < Ly - 1;++row){

         //first create right renormalized operator

         //paste top environment on
         tmp5.clear();
         Contract(one,Environment::t[myID][row][Lx - 1],shape(1),Environment::U[myID](row,Lx-1),shape(1),zero,tmp5);

         //then bottom enviroment
         TArray<complex<double>,6> tmp6;
         Contract(one,tmp5,shape(3),Environment::b[myID][row-1][Lx-1],shape(1),zero,tmp6);

         //move to a ZArray<3> object
         RO[Lx - 2] = tmp6.reshape_clear(shape(Environment::t[myID][row][Lx - 1].shape(0),Environment::U[myID](row,Lx-1).shape(0),Environment::b[myID][row-1][Lx - 1].shape(0)));

         ZArray<4> I4;
         ZArray<4> I4bis;

         //now construct the middle operators
         for(int col = Lx-2;col > 0;--col){

            I4.clear();
            Contract(one,Environment::t[myID][row][col],shape(2),RO[col],shape(0),zero,I4);

            enum {i,j,k,o,m,n};

            I4bis.clear();
            Contract(one,I4,shape(i,j,k,o),Environment::U[myID](row,col),shape(m,j,n,k),zero,I4bis,shape(i,m,n,o));

            RO[col-1].clear();
            Contract(one,I4bis,shape(2,3),Environment::b[myID][row-1][col],shape(1,2),zero,RO[col-1]);

         }

         // --- now move from left to right to get the expecation value of the interactions ---
         // --- First construct the left going operators for the first site -----

         // 1) Sx

         //paste top environment on local Sx
         tmp5.clear();
         Contract(one,Environment::t[myID][row][0],shape(1),Environment::Sx[myID](row,0),shape(1),zero,tmp5);

         //then bottom enviroment on that
         tmp6.clear();
         Contract(one,tmp5,shape(3),Environment::b[myID][row-1][0],shape(1),zero,tmp6);

         //move to a ZArray<3> object: order (top-env,peps-row,bottom-env)
         LOSx = tmp6.reshape_clear(shape(Environment::t[myID][row][0].shape(2),Environment::Sx[myID](row,0).shape(3),Environment::b[myID][row-1][0].shape(2)));

         // 2) Sy

         //paste top environment on local Sy
         tmp5.clear();
         Contract(one,Environment::t[myID][row][0],shape(1),Environment::Sy[myID](row,0),shape(1),zero,tmp5);

         //then bottom enviroment on that
         tmp6.clear();
         Contract(one,tmp5,shape(3),Environment::b[myID][row-1][0],shape(1),zero,tmp6);

         //move to a ZArray<3> object: order (top-env,peps-row,bottom-env)
         LOSy = tmp6.reshape_clear(shape(Environment::t[myID][row][0].shape(2),Environment::Sy[myID](row,0).shape(3),Environment::b[myID][row-1][0].shape(2)));

         // 3) Sz

         //paste top environment on local Sz
         tmp5.clear();
         Contract(one,Environment::t[myID][row][0],shape(1),Environment::Sz[myID](row,0),shape(1),zero,tmp5);

         //then bottom enviroment on that
         tmp6.clear();
         Contract(one,tmp5,shape(3),Environment::b[myID][row-1][0],shape(1),zero,tmp6);

         //move to a ZArray<3> object: order (top-env,peps-row,bottom-env)
         LOSz = tmp6.reshape_clear(shape(Environment::t[myID][row][0].shape(2),Environment::Sz[myID](row,0).shape(3),Environment::b[myID][row-1][0].shape(2)));

         // 4) 1 -- finally construct left renormalized operator with unity

         //paste top environment on local Sz
         tmp5.clear();
         Contract(one,Environment::t[myID][row][0],shape(1),Environment::U[myID](row,0),shape(1),zero,tmp5);

         //then bottom enviroment on that
         tmp6.clear();
         Contract(one,tmp5,shape(3),Environment::b[myID][row-1][0],shape(1),zero,tmp6);

         //move to a ZArray<3> object: order (top-env,peps-row,bottom-env)
         LOU = tmp6.reshape_clear(shape(Environment::t[myID][row][0].shape(2),Environment::U[myID](row,0).shape(3),Environment::b[myID][row-1][0].shape(2)));

         //now contract x,y and z with R for local expectation values:
         dim = RO[0].size();

         auxvec[myID][row*Lx][0] = blas::dot(dim,LOSx.data(),1,RO[0].data(),1);
         auxvec[myID][row*Lx][1] = blas::dot(dim,LOSy.data(),1,RO[0].data(),1);
         auxvec[myID][row*Lx][2] = blas::dot(dim,LOSz.data(),1,RO[0].data(),1);

         // --- now for the middle sites, close down the operators on the left and construct new ones --- 
         for(int col = 1;col < Lx - 1;++col){

            //first add top to the right side, put it in I4
            I4.clear();
            Contract(one,Environment::t[myID][row][col],shape(2),RO[col],shape(0),zero,I4);

            enum {i,j,k,o,m,n};

            //1) close down LOSx with Sx
            I4bis.clear();
            Contract(one,I4,shape(i,j,k,o),Environment::Sx[myID](row,col),shape(m,j,n,k),zero,I4bis,shape(i,m,n,o));

            Contract(one,I4bis,shape(2,3),Environment::b[myID][row-1][col],shape(1,2),zero,RO[col-1]);

            //expectation value:
            energy += Dot(LOSx,RO[col-1]);

            //2) close down LOSy with Sy
            I4bis.clear();
            Contract(one,I4,shape(i,j,k,o),Environment::Sy[myID](row,col),shape(m,j,n,k),zero,I4bis,shape(i,m,n,o));

            Contract(one,I4bis,shape(2,3),Environment::b[myID][row-1][col],shape(1,2),zero,RO[col-1]);

            //expectation value:
            energy += Dot(LOSy,RO[col-1]);

            //3) finally close down LOSz with Sz
            I4bis.clear();
            Contract(one,I4,shape(i,j,k,o),Environment::Sz[myID](row,col),shape(m,j,n,k),zero,I4bis,shape(i,m,n,o));

            Contract(one,I4bis,shape(2,3),Environment::b[myID][row-1][col],shape(1,2),zero,RO[col-1]);

            //expectation value:
            energy += Dot(LOSz,RO[col-1]);

            // now construct the new left going renormalized operators

            //first attach top to left unity
            I4.clear();
            Contract(one,Environment::t[myID][row][col],shape(0),LOU,shape(0),zero,I4);

            // 1) construct new left Sx operator
            I4bis.clear();
            Contract(one,I4,shape(i,j,k,o),Environment::Sx[myID](row,col),shape(k,i,m,n),zero,I4bis,shape(j,n,o,m));

            LOSx.clear();
            Contract(one,I4bis,shape(2,3),Environment::b[myID][row-1][col],shape(0,1),zero,LOSx);

            // 2) construct new left Sy operator
            Contract(one,I4,shape(i,j,k,o),Environment::Sy[myID](row,col),shape(k,i,m,n),zero,I4bis,shape(j,n,o,m));

            LOSy.clear();
            Contract(one,I4bis,shape(2,3),Environment::b[myID][row-1][col],shape(0,1),zero,LOSy);

            // 3) construct new left Sz operator
            Contract(one,I4,shape(i,j,k,o),Environment::Sz[myID](row,col),shape(k,i,m,n),zero,I4bis,shape(j,n,o,m));

            LOSz.clear();
            Contract(one,I4bis,shape(2,3),Environment::b[myID][row-1][col],shape(0,1),zero,LOSz);

            // 4) finally construct new left unity
            Contract(one,I4,shape(i,j,k,o),Environment::U[myID](row,col),shape(k,i,m,n),zero,I4bis,shape(j,n,o,m));

            LOU.clear();
            Contract(one,I4bis,shape(2,3),Environment::b[myID][row-1][col],shape(0,1),zero,LOU);

            //now contract x,y and z with R for local expectation values:
            dim = RO[col].size();

            auxvec[myID][row*Lx + col][0] = blas::dot(dim,LOSx.data(),1,RO[col].data(),1);
            auxvec[myID][row*Lx + col][1] = blas::dot(dim,LOSy.data(),1,RO[col].data(),1);
            auxvec[myID][row*Lx + col][2] = blas::dot(dim,LOSz.data(),1,RO[col].data(),1);

         }

         //last site on the right: close down on the incomings

         //1) first LSx with Sx

         //paste top environment on
         tmp5.clear();
         Contract(one,Environment::t[myID][row][Lx - 1],shape(1),Environment::Sx[myID](row,Lx-1),shape(1),zero,tmp5);

         //then bottom enviroment
         Contract(one,tmp5,shape(3),Environment::b[myID][row-1][Lx-1],shape(1),zero,tmp6);

         //move to a ZArray<3> object
         RO[Lx - 2] = tmp6.reshape_clear(shape(Environment::t[myID][row][Lx - 1].shape(0),Environment::Sx[myID](row,Lx-1).shape(0),Environment::b[myID][row-1][Lx - 1].shape(0)));

         //add to energy
         energy += Dot(LOSx,RO[Lx - 2]);

         //local expectation value
         dim = RO[Lx-2].size();
         auxvec[myID][row*Lx + Lx-1][0] = blas::dot(dim,LOU.data(),1,RO[Lx-2].data(),1);

         //2) then LSy with Sy

         //paste top environment on
         tmp5.clear();
         Contract(one,Environment::t[myID][row][Lx - 1],shape(1),Environment::Sy[myID](row,Lx-1),shape(1),zero,tmp5);

         //then bottom enviroment
         Contract(one,tmp5,shape(3),Environment::b[myID][row-1][Lx-1],shape(1),zero,tmp6);

         //move to a DArray<3> object
         RO[Lx - 2] = tmp6.reshape_clear(shape(Environment::t[myID][row][Lx - 1].shape(0),Environment::Sy[myID](row,Lx-1).shape(0),Environment::b[myID][row-1][Lx - 1].shape(0)));

         //add to value
         energy += Dot(LOSy,RO[Lx - 2]);

         //local expectation value
         auxvec[myID][row*Lx + Lx-1][1] = blas::dot(dim,LOU.data(),1,RO[Lx-2].data(),1);

         //3) then LZz with Sz

         //paste top environment on
         tmp5.clear();
         Contract(one,Environment::t[myID][row][Lx - 1],shape(1),Environment::Sz[myID](row,Lx-1),shape(1),zero,tmp5);

         //then bottom enviroment
         Contract(one,tmp5,shape(3),Environment::b[myID][row-1][Lx-1],shape(1),zero,tmp6);

         //move to a DArray<3> object
         RO[Lx - 2] = tmp6.reshape_clear(shape(Environment::t[myID][row][Lx - 1].shape(0),Environment::Sz[myID](row,Lx-1).shape(0),Environment::b[myID][row-1][Lx - 1].shape(0)));

         //add to value
         energy += Dot(LOSz,RO[Lx - 2]);

         //local expectation value
         auxvec[myID][row*Lx + Lx-1][2] = blas::dot(dim,LOU.data(),1,RO[Lx-2].data(),1);

      }

      // -- (3) -- || top row = Ly-1: again similar to overlap calculation

      //first construct the right renormalized operators

      //tmp comes out index (t,b)
      tmp4.clear();
      Contract(one,Environment::t[myID][Ly-2][Lx - 1],shape(1),Environment::b[myID][Ly-2][Lx - 1],shape(1),zero,tmp4);

      //reshape tmp to a 2-index array
      R[Lx - 2] = tmp4.reshape_clear(shape(Environment::t[myID][Ly-2][Lx - 1].shape(0),Environment::b[myID][Ly-2][Lx - 1].shape(0)));

      //now construct the rest
      for(int col = Lx - 2;col > 0;--col){

         tmp3.clear();
         Contract(one,Environment::t[myID][Ly-2][col],shape(2),R[col],shape(0),zero,tmp3);

         R[col - 1].clear();
         Contract(one,tmp3,shape(1,2),Environment::b[myID][Ly-2][col],shape(1,2),zero,R[col-1]);

      }

      //construct the left going operators on the first top site

      //first Sx

      //tmp comes out index (t,b)
      tmp5.clear();
      Contract(one,Environment::Sx[myID](Ly-1,0),shape(2),Environment::b[myID][Ly-2][0],shape(1),zero,tmp5);

      LSx = tmp5.reshape_clear(shape(Environment::Sx[myID](Ly-1,0).shape(3),Environment::b[myID][Ly-2][0].shape(2)));

      //then Sy

      //tmp5 comes out index (t,b)
      Contract(one,Environment::Sy[myID](Ly-1,0),shape(2),Environment::b[myID][Ly-2][0],shape(1),zero,tmp5);

      LSy = tmp5.reshape_clear(shape(Environment::Sy[myID](Ly-1,0).shape(3),Environment::b[myID][Ly-2][0].shape(2)));

      //then Sz 

      //tmp comes out index (t,b)
      Contract(one,Environment::Sz[myID](Ly-1,0),shape(2),Environment::b[myID][Ly-2][0],shape(1),zero,tmp5);

      LSz = tmp5.reshape_clear(shape(Environment::Sz[myID](Ly-1,0).shape(3),Environment::b[myID][Ly-2][0].shape(2)));

      //and finally unity
      Contract(one,Environment::t[myID][Ly-2][0],shape(1),Environment::b[myID][Ly-2][0],shape(1),zero,tmp4);

      LU = tmp4.reshape_clear(shape(Environment::t[myID][Ly-2][0].shape(2),Environment::b[myID][Ly-2][0].shape(2)));

      dim = R[0].size();

      //now contract x,y and z with R for local expectation values:
      auxvec[myID][(Ly-1)*Lx][0] = blas::dot(dim,LSx.data(),1,R[0].data(),1);
      auxvec[myID][(Ly-1)*Lx][1] = blas::dot(dim,LSy.data(),1,R[0].data(),1);
      auxvec[myID][(Ly-1)*Lx][2] = blas::dot(dim,LSz.data(),1,R[0].data(),1);

      //middle of the chain:
      for(int col = 1;col < Lx-1;++col){

         //first close down the x,y and z terms from the previous site for the energy

         //construct the right intermediate contraction (paste bottom to right)
         tmp3.clear();
         Contract(one,Environment::b[myID][Ly-2][col],shape(2),R[col],shape(1),zero,tmp3);

         // 1) paste Sx to the right
         int M = Environment::Sx[myID](Ly-1,col).shape(0);
         int N = tmp3.shape(0);
         int K = tmp3.shape(1) * tmp3.shape(2);

         blas::gemm(CblasRowMajor, CblasNoTrans, CblasTrans, M, N, K, one, Environment::Sx[myID](Ly-1,col).data(),K,tmp3.data(),K,zero,R[col-1].data(),N);

         //contract with left Sx
         energy += Dot(LSx,R[col - 1]);

         // 2) paste Sy to the right
         M = Environment::Sy[myID](Ly-1,col).shape(0);
         N = tmp3.shape(0);
         K = tmp3.shape(1) * tmp3.shape(2);

         blas::gemm(CblasRowMajor, CblasNoTrans, CblasTrans, M, N, K, one, Environment::Sy[myID](Ly-1,col).data(),K,tmp3.data(),K,zero,R[col-1].data(),N);

         //contract with left Sy
         energy += Dot(LSy,R[col - 1]);

         // 3) paste Sz to the right
         M = Environment::Sz[myID](Ly-1,col).shape(0);
         N = tmp3.shape(0);
         K = tmp3.shape(1) * tmp3.shape(2);

         blas::gemm(CblasRowMajor, CblasNoTrans, CblasTrans, M, N, K, one, Environment::Sz[myID](Ly-1,col).data(),K,tmp3.data(),K,zero,R[col-1].data(),N);

         //contract with left Sy
         energy += Dot(LSz,R[col - 1]);

         //construct left renormalized operators for next site: first paste bottom to Left unity
         tmp3.clear();
         Contract(one,LU,shape(1),Environment::b[myID][Ly-2][col],shape(0),zero,tmp3);

         // 1) construct new Sx left operator
         LSx.resize(Environment::Sx[myID](Ly-1,col).shape(3),Environment::b[myID][Ly-2][col].shape(2));

         M = Environment::Sx[myID](Ly-1,col).shape(3);
         N = tmp3.shape(2);
         K = tmp3.shape(0) * tmp3.shape(1);

         blas::gemm(CblasRowMajor, CblasTrans, CblasNoTrans, M, N, K, one, Environment::Sx[myID](Ly-1,col).data(),M,tmp3.data(),N,zero,LSx.data(),N);

         // 2) construct new Sy left operator
         LSy.resize(Environment::Sy[myID](Ly-1,col).shape(3),Environment::b[myID][Ly-2][col].shape(2));
         blas::gemm(CblasRowMajor, CblasTrans, CblasNoTrans, M, N, K, one, Environment::Sy[myID](Ly-1,col).data(),M,tmp3.data(),N,zero,LSy.data(),N);

         // 3) construct new Sz left operator
         LSz.resize(Environment::Sz[myID](Ly-1,col).shape(3),Environment::b[myID][Ly-2][col].shape(2));
         blas::gemm(CblasRowMajor, CblasTrans, CblasNoTrans, M, N, K, one, Environment::Sz[myID](Ly-1,col).data(),M,tmp3.data(),N,zero,LSz.data(),N);

         // 4) finally construct new unity on the left
         LU.resize(Environment::U[myID](Ly-1,col).shape(3),Environment::b[myID][Ly-2][col].shape(2));
         blas::gemm(CblasRowMajor, CblasTrans, CblasNoTrans, M, N, K, one, Environment::t[myID][Ly-2][col].data(),M,tmp3.data(),N,zero,LU.data(),N);

         //now contract x,y and z with R for local expectation values:
         dim = R[col].size();

         auxvec[myID][(Ly-1)*Lx + col][0] = blas::dot(dim,LSx.data(),1,R[col].data(),1);
         auxvec[myID][(Ly-1)*Lx + col][1] = blas::dot(dim,LSy.data(),1,R[col].data(),1);
         auxvec[myID][(Ly-1)*Lx + col][2] = blas::dot(dim,LSz.data(),1,R[col].data(),1);

      }

      //finally close down on last top site

      //first calculate overlap
      overlap = Dot(LU,R[Lx-2]);

      //1) Sx to close down LSx

      //tmp comes out index (t,b)
      tmp5.clear();
      Contract(one,Environment::Sx[myID](Ly-1,Lx-1),shape(2),Environment::b[myID][Ly-2][Lx - 1],shape(1),zero,tmp5);

      //reshape tmp to a 2-index array
      R[Lx - 2] = tmp5.reshape_clear(shape(Environment::Sx[myID](Ly-1,Lx-1).shape(0),Environment::b[myID][Ly-2][Lx - 1].shape(0)));

      //energy
      energy += Dot(LSx,R[Lx-2]);

      //local expectation value
      auxvec[myID][(Ly-1)*Lx + Lx-1][0] = blas::dot(dim,LU.data(),1,R[Lx-2].data(),1);

      //2) Sy to close down Ly
      Contract(one,Environment::Sy[myID](Ly-1,Lx-1),shape(2),Environment::b[myID][Ly-2][Lx - 1],shape(1),zero,tmp5);

      //reshape tmp to a 2-index array
      R[Lx - 2] = tmp5.reshape_clear(shape(Environment::Sy[myID](Ly-1,Lx-1).shape(0),Environment::b[myID][Ly-2][Lx - 1].shape(0)));

      //energy
      energy += Dot(LSy,R[Lx-2]);

      //local expectation value
      auxvec[myID][(Ly-1)*Lx + Lx-1][1] = blas::dot(dim,LU.data(),1,R[Lx-2].data(),1);

      //3) Sz to close down Lz

      //tmp comes out index (t,b)
      Contract(one,Environment::Sz[myID](Ly-1,Lx-1),shape(2),Environment::b[myID][Ly-2][Lx - 1],shape(1),zero,tmp5);

      //reshape tmp to a 2-index array
      R[Lx - 2] = tmp5.reshape_clear(shape(Environment::Sz[myID](Ly-1,Lx-1).shape(0),Environment::b[myID][Ly-2][Lx - 1].shape(0)));

      //energy
      energy += Dot(LSz,R[Lx-2]);

      //local expectation value
      auxvec[myID][(Ly-1)*Lx + Lx-1][2] = blas::dot(dim,LU.data(),1,R[Lx-2].data(),1);

      // #################################################################
      // ### ----               SET THE PROPERTIES                ---- ### 
      // #################################################################

      EL = energy/overlap;

      //set VL: Sx, Sy and Sz on every site: ready to construct auxiliary expectation values quickly!
      for(int k = 0;k < Trotter::n_trot;++k)
         for(int r = 0;r < 3;++r){

            VL[r*Trotter::n_trot + k] = auxvec[myID][0][r] * Trotter::gV()(k,0);


            for(int row = 0;row < Ly;++row)
               for(int col = 0;col < Lx;++col)
                  VL[r*Trotter::n_trot + k] += auxvec[myID][row*Lx + col][r] * Trotter::gV()(k,row*Lx + col);

            VL[r*Trotter::n_trot + k] /= overlap;

         }

   }
   else{//VERTICAL: Left to right

      // #################################################################
      // ### ---- from left to right: contract in mps/mpo fashion ---- ### 
      // #################################################################

      //calculate the single layer contractions first:
      Environment::U[myID].fill('V',peps,*this);

      Environment::Sx[myID].fill('V',peps,Sx,*this);
      Environment::Sy[myID].fill('V',peps,Sy,*this);
      Environment::Sz[myID].fill('V',peps,Sz,*this);

      //first construct the top and bottom (horizontal) environment layers
      Environment::calc_env('V',peps,*this);

      // -- (1) -- || left column: similar to overlap calculation

      //first construct the right renormalized operators
      vector< ZArray<2> > R(Ly - 1);

      //first the rightmost operator
      ZArray<4> tmp4;
      ZArray<3> tmp3;

      //tmp comes out index (t,b)
      Contract(one,Environment::r[myID][0][Ly - 1],shape(1),Environment::l[myID][0][Ly - 1],shape(1),zero,tmp4);

      //reshape tmp to a 2-index array
      R[Ly - 2] = tmp4.reshape_clear(shape(Environment::r[myID][0][Ly - 1].shape(0),Environment::l[myID][0][Ly - 1].shape(0)));

      //now construct the rest
      for(int row = Ly - 2;row > 0;--row){

         tmp3.clear();
         Contract(one,Environment::r[myID][0][row],shape(2),R[row],shape(0),zero,tmp3);

         Contract(one,tmp3,shape(1,2),Environment::l[myID][0][row],shape(1,2),zero,R[row-1]);

      }

      //4 left going operators: Sx, Sy, Sz, and 1
      ZArray<2> LSx;
      ZArray<2> LSy;
      ZArray<2> LSz;
      ZArray<2> LU;

      TArray<complex<double>,5> tmp5;

      //tmp comes out index (r,l)
      Contract(one,Environment::r[myID][0][0],shape(1),Environment::Sx[myID](0,0),shape(1),zero,tmp5);

      LSx = tmp5.reshape_clear(shape(Environment::r[myID][0][0].shape(2),Environment::Sx[myID](0,0).shape(3)));

      //then Sy
      Contract(one,Environment::r[myID][0][0],shape(1),Environment::Sy[myID](0,0),shape(1),zero,tmp5);

      LSy = tmp5.reshape_clear(shape(Environment::r[myID][0][0].shape(2),Environment::Sy[myID](0,0).shape(3)));

      //then Sz
      Contract(one,Environment::r[myID][0][0],shape(1),Environment::Sz[myID](0,0),shape(1),zero,tmp5);

      LSz = tmp5.reshape_clear(shape(Environment::r[myID][0][0].shape(2),Environment::Sz[myID](0,0).shape(3)));

      //finally unity
      Contract(one,Environment::r[myID][0][0],shape(1),Environment::U[myID](0,0),shape(1),zero,tmp5);

      LU = tmp5.reshape_clear(shape(Environment::r[myID][0][0].shape(2),Environment::Sz[myID](0,0).shape(3)));

      int dim = R[0].size();

      //now contract x,y and z with R for local expectation values:
      auxvec[myID][0][0] = blas::dot(dim,LSx.data(),1,R[0].data(),1);
      auxvec[myID][0][1] = blas::dot(dim,LSy.data(),1,R[0].data(),1);
      auxvec[myID][0][2] = blas::dot(dim,LSz.data(),1,R[0].data(),1);

      //now for the middle terms
      for(int row = 1;row < Ly - 1;++row){

         //first close down the x,y and z terms from the previous site for the energy

         //construct the right intermediate contraction (paste top to right)
         tmp3.clear();
         Contract(one,Environment::r[myID][0][row],shape(2),R[row],shape(0),zero,tmp3);

         // 1) paste Sx to the right
         int M = tmp3.shape(0);
         int N = Environment::Sx[myID](row,0).shape(0);
         int K = tmp3.shape(1) * tmp3.shape(2);

         blas::gemm(CblasRowMajor, CblasNoTrans, CblasTrans, M, N, K, one, tmp3.data(),K,Environment::Sx[myID](row,0).data(),K,zero,R[row-1].data(),N);

         //contract with left Sx
         energy += Dot(LSx,R[row - 1]);

         // 2) then paste Sy to the right
         blas::gemm(CblasRowMajor, CblasNoTrans, CblasTrans, M, N, K, one, tmp3.data(),K,Environment::Sy[myID](row,0).data(),K,zero,R[row-1].data(),N);

         //contract with left Sy
         energy += Dot(LSy,R[row - 1]);

         // 3) finally Sz
         blas::gemm(CblasRowMajor, CblasNoTrans, CblasTrans, M, N, K, one, tmp3.data(),K,Environment::Sz[myID](row,0).data(),K,zero,R[row-1].data(),N);

         //contract with left Sz
         energy += Dot(LSz,R[row - 1]);

         //construct left renormalized operators for next site: first paste top to Left unity
         tmp3.clear();
         Contract(one,LU,shape(0),Environment::r[myID][0][row],shape(0),zero,tmp3);

         // 1) construct new Sx left operator
         LSx.resize(Environment::r[myID][0][row].shape(2),Environment::Sx[myID](row,0).shape(3));

         M = tmp3.shape(2);
         N = Environment::Sx[myID](row,0).shape(3);
         K = tmp3.shape(0) * tmp3.shape(1);

         blas::gemm(CblasRowMajor, CblasTrans, CblasNoTrans, M, N, K, one, tmp3.data(),M,Environment::Sx[myID](row,0).data(),N,zero,LSx.data(),N);

         // 2) construct new Sy left operator
         LSy.resize(Environment::r[myID][0][row].shape(2),Environment::Sy[myID](row,0).shape(3));
         blas::gemm(CblasRowMajor, CblasTrans, CblasNoTrans, M, N, K, one, tmp3.data(),M,Environment::Sy[myID](row,0).data(),N,zero,LSy.data(),N);

         // 3) construct new Sz left operator
         LSz.resize(Environment::r[myID][0][row].shape(2),Environment::Sz[myID](row,0).shape(3));
         blas::gemm(CblasRowMajor, CblasTrans, CblasNoTrans, M, N, K, one, tmp3.data(),M,Environment::Sz[myID](row,0).data(),N,zero,LSz.data(),N);

         //now contract x,y and z with R for local expectation values:
         dim = R[row].size();

         auxvec[myID][row*Lx][0] = blas::dot(dim,LSx.data(),1,R[row].data(),1);
         auxvec[myID][row*Lx][1] = blas::dot(dim,LSy.data(),1,R[row].data(),1);
         auxvec[myID][row*Lx][2] = blas::dot(dim,LSz.data(),1,R[row].data(),1);

         // 4) finally construct new unity on the left
         LU.resize(Environment::r[myID][0][row].shape(2),Environment::U[myID](row,0).shape(3));
         blas::gemm(CblasRowMajor, CblasTrans, CblasNoTrans, M, N, K, one, tmp3.data(),M,Environment::U[myID](row,0).data(),N,zero,LU.data(),N);

      }

      //last site of left column: close down the left x,y and z

      //1) Sx to close down LSx
      Contract(one,Environment::r[myID][0][Ly-1],shape(1),Environment::Sx[myID](Ly-1,0),shape(1),zero,tmp5);

      R[Ly-2] = tmp5.reshape_clear(shape(Environment::r[myID][0][Ly-1].shape(0),Environment::Sx[myID](Ly-1,0).shape(0)));

      energy += Dot(LSx,R[Ly-2]);

      dim = R[Ly-2].size();
      auxvec[myID][(Ly-1)*Lx][0] = blas::dot(dim,LU.data(),1,R[Ly-2].data(),1);

      //2) Sy to close down LSy
      Contract(one,Environment::r[myID][0][Ly-1],shape(1),Environment::Sy[myID](Ly-1,0),shape(1),zero,tmp5);

      R[Ly-2] = tmp5.reshape_clear(shape(Environment::r[myID][0][Ly-1].shape(0),Environment::Sy[myID](Ly-1,0).shape(0)));

      energy += Dot(LSy,R[Ly-2]);

      auxvec[myID][(Ly-1)*Lx][1] = blas::dot(dim,LU.data(),1,R[Ly-2].data(),1);

      //3) Sz to close down Lz
      Contract(one,Environment::r[myID][0][Ly-1],shape(1),Environment::Sz[myID](Ly-1,0),shape(1),zero,tmp5);

      R[Ly-2] = tmp5.reshape_clear(shape(Environment::r[myID][0][Ly-1].shape(0),Environment::Sz[myID](Ly-1,0).shape(0)));

      energy += Dot(LSz,R[Ly-2]);

      auxvec[myID][(Ly-1)*Lx][2] = blas::dot(dim,LU.data(),1,R[Ly-2].data(),1);

      // -- (2) -- now move from left to right calculating everything like an MPO/MPS expectation value

      //Right renormalized operators
      vector< TArray<complex<double>,3> > RO(Ly - 1);

      //4 left renormalized operators needed
      TArray<complex<double>,3> LOSx;
      TArray<complex<double>,3> LOSy;
      TArray<complex<double>,3> LOSz;
      TArray<complex<double>,3> LOU;

      for(int col = 1;col < Lx - 1;++col){

         //first create right renormalized operator

         //paste right environment on
         tmp5.clear();
         Contract(one,Environment::r[myID][col][Ly - 1],shape(1),Environment::U[myID](Ly-1,col),shape(1),zero,tmp5);

         //then left enviroment
         TArray<complex<double>,6> tmp6;
         Contract(one,tmp5,shape(3),Environment::l[myID][col-1][Ly-1],shape(1),zero,tmp6);

         //move to a ZArray<3> object
         RO[Ly - 2] = tmp6.reshape_clear(shape(Environment::r[myID][col][Ly - 1].shape(0),Environment::U[myID](Ly-1,col).shape(0),Environment::l[myID][col-1][Ly - 1].shape(0)));

         ZArray<4> I4;
         ZArray<4> I4bis;

         //now construct the middle operators
         for(int row = Ly-2;row > 0;--row){

            I4.clear();
            Contract(one,Environment::r[myID][col][row],shape(2),RO[row],shape(0),zero,I4);

            enum {i,j,k,o,m,n};

            I4bis.clear();
            Contract(one,I4,shape(i,j,k,o),Environment::U[myID](row,col),shape(m,j,n,k),zero,I4bis,shape(i,m,n,o));

            RO[row-1].clear();
            Contract(one,I4bis,shape(2,3),Environment::l[myID][col-1][row],shape(1,2),zero,RO[row-1]);

         }

         // --- now move from left to right to get the expecation value of the interactions ---
         // --- First construct the left going operators for the first site -----

         // 1) Sx

         //paste top environment on local Sx
         tmp5.clear();
         Contract(one,Environment::r[myID][col][0],shape(1),Environment::Sx[myID](0,col),shape(1),zero,tmp5);

         //then bottom enviroment on that
         tmp6.clear();
         Contract(one,tmp5,shape(3),Environment::l[myID][col-1][0],shape(1),zero,tmp6);

         //move to a ZArray<3> object: order (top-env,peps-row,bottom-env)
         LOSx = tmp6.reshape_clear(shape(Environment::r[myID][col][0].shape(2),Environment::Sx[myID](0,col).shape(3),Environment::l[myID][col-1][0].shape(2)));

         // 2) Sy

         //paste top environment on local Sy
         tmp5.clear();
         Contract(one,Environment::r[myID][col][0],shape(1),Environment::Sy[myID](0,col),shape(1),zero,tmp5);

         //then bottom enviroment on that
         tmp6.clear();
         Contract(one,tmp5,shape(3),Environment::l[myID][col-1][0],shape(1),zero,tmp6);

         //move to a ZArray<3> object: order (top-env,peps-row,bottom-env)
         LOSy = tmp6.reshape_clear(shape(Environment::r[myID][col][0].shape(2),Environment::Sy[myID](0,col).shape(3),Environment::l[myID][col-1][0].shape(2)));

         // 3) Sz

         //paste top environment on local Sz
         tmp5.clear();
         Contract(one,Environment::r[myID][col][0],shape(1),Environment::Sz[myID](0,col),shape(1),zero,tmp5);

         //then bottom enviroment on that
         tmp6.clear();
         Contract(one,tmp5,shape(3),Environment::l[myID][col-1][0],shape(1),zero,tmp6);

         //move to a ZArray<3> object: order (top-env,peps-row,bottom-env)
         LOSz = tmp6.reshape_clear(shape(Environment::r[myID][col][0].shape(2),Environment::Sz[myID](0,col).shape(3),Environment::r[myID][col-1][0].shape(2)));

         // 4) 1 -- finally construct left renormalized operator with unity

         //paste top environment on local Sz
         tmp5.clear();
         Contract(one,Environment::r[myID][col][0],shape(1),Environment::U[myID](0,col),shape(1),zero,tmp5);

         //then bottom enviroment on that
         tmp6.clear();
         Contract(one,tmp5,shape(3),Environment::l[myID][col-1][0],shape(1),zero,tmp6);

         //move to a ZArray<3> object: order (top-env,peps-row,bottom-env)
         LOU = tmp6.reshape_clear(shape(Environment::r[myID][col][0].shape(2),Environment::U[myID](0,col).shape(3),Environment::l[myID][col-1][0].shape(2)));

         //now contract x,y and z with R for local expectation values:
         dim = RO[0].size();

         auxvec[myID][col][0] = blas::dot(dim,LOSx.data(),1,RO[0].data(),1);
         auxvec[myID][col][1] = blas::dot(dim,LOSy.data(),1,RO[0].data(),1);
         auxvec[myID][col][2] = blas::dot(dim,LOSz.data(),1,RO[0].data(),1);

         // --- now for the middle sites, close down the operators on the left and construct new ones --- 
         for(int row = 1;row < Ly - 1;++row){

            //first add top to the right side, put it in I4
            I4.clear();
            Contract(one,Environment::r[myID][col][row],shape(2),RO[row],shape(0),zero,I4);

            enum {i,j,k,o,m,n};

            //1) close down LOSx with Sx
            I4bis.clear();
            Contract(one,I4,shape(i,j,k,o),Environment::Sx[myID](row,col),shape(m,j,n,k),zero,I4bis,shape(i,m,n,o));

            Contract(one,I4bis,shape(2,3),Environment::l[myID][col-1][row],shape(1,2),zero,RO[row-1]);

            //expectation value:
            energy += Dot(LOSx,RO[row-1]);

            //2) close down LOSy with Sy
            I4bis.clear();
            Contract(one,I4,shape(i,j,k,o),Environment::Sy[myID](row,col),shape(m,j,n,k),zero,I4bis,shape(i,m,n,o));

            Contract(one,I4bis,shape(2,3),Environment::l[myID][col-1][row],shape(1,2),zero,RO[row-1]);

            //expectation value:
            energy += Dot(LOSy,RO[row-1]);

            //3) finally close down LOSz with Sz
            I4bis.clear();
            Contract(one,I4,shape(i,j,k,o),Environment::Sz[myID](row,col),shape(m,j,n,k),zero,I4bis,shape(i,m,n,o));

            Contract(one,I4bis,shape(2,3),Environment::l[myID][col-1][row],shape(1,2),zero,RO[row-1]);

            //expectation value:
            energy += Dot(LOSz,RO[row-1]);

            // now construct the new left going renormalized operators

            //first attach top to left unity
            I4.clear();
            Contract(one,Environment::r[myID][col][row],shape(0),LOU,shape(0),zero,I4);

            // 1) construct new left Sx operator
            I4bis.clear();
            Contract(one,I4,shape(i,j,k,o),Environment::Sx[myID](row,col),shape(k,i,m,n),zero,I4bis,shape(j,n,o,m));

            LOSx.clear();
            Contract(one,I4bis,shape(2,3),Environment::l[myID][col-1][row],shape(0,1),zero,LOSx);

            // 2) construct new left Sy operator
            Contract(one,I4,shape(i,j,k,o),Environment::Sy[myID](row,col),shape(k,i,m,n),zero,I4bis,shape(j,n,o,m));

            LOSy.clear();
            Contract(one,I4bis,shape(2,3),Environment::l[myID][col-1][row],shape(0,1),zero,LOSy);

            // 3) construct new left Sz operator
            Contract(one,I4,shape(i,j,k,o),Environment::Sz[myID](row,col),shape(k,i,m,n),zero,I4bis,shape(j,n,o,m));

            LOSz.clear();
            Contract(one,I4bis,shape(2,3),Environment::l[myID][col-1][row],shape(0,1),zero,LOSz);

            // 4) finally construct new left unity
            Contract(one,I4,shape(i,j,k,o),Environment::U[myID](row,col),shape(k,i,m,n),zero,I4bis,shape(j,n,o,m));

            LOU.clear();
            Contract(one,I4bis,shape(2,3),Environment::l[myID][col-1][row],shape(0,1),zero,LOU);

            //now contract x,y and z with R for local expectation values:
            dim = RO[row].size();

            auxvec[myID][row*Lx + col][0] = blas::dot(dim,LOSx.data(),1,RO[row].data(),1);
            auxvec[myID][row*Lx + col][1] = blas::dot(dim,LOSy.data(),1,RO[row].data(),1);
            auxvec[myID][row*Lx + col][2] = blas::dot(dim,LOSz.data(),1,RO[row].data(),1);

         }

         //last site on the right: close down on the incomings

         //1) first LSx with Sx

         //paste top environment on
         tmp5.clear();
         Contract(one,Environment::r[myID][col][Ly - 1],shape(1),Environment::Sx[myID](Ly-1,col),shape(1),zero,tmp5);

         //then bottom enviroment
         Contract(one,tmp5,shape(3),Environment::l[myID][col-1][Ly-1],shape(1),zero,tmp6);

         //move to a ZArray<3> object
         RO[Ly - 2] = tmp6.reshape_clear(shape(Environment::r[myID][col][Ly - 1].shape(0),Environment::Sx[myID](Ly-1,col).shape(0),Environment::l[myID][col-1][Ly - 1].shape(0)));

         //add to energy
         energy += Dot(LOSx,RO[Ly - 2]);

         //local expectation value
         dim = RO[Ly-2].size();
         auxvec[myID][(Ly-1)*Lx + col][0] = blas::dot(dim,LOU.data(),1,RO[Ly-2].data(),1);

         //2) then LSy with Sy

         //paste top environment on
         tmp5.clear();
         Contract(one,Environment::r[myID][col][Ly - 1],shape(1),Environment::Sy[myID](Ly-1,col),shape(1),zero,tmp5);

         //then bottom enviroment
         Contract(one,tmp5,shape(3),Environment::l[myID][col-1][Ly-1],shape(1),zero,tmp6);

         //move to a DArray<3> object
         RO[Ly - 2] = tmp6.reshape_clear(shape(Environment::r[myID][col][Ly - 1].shape(0),Environment::Sy[myID](Ly-1,col).shape(0),Environment::l[myID][col-1][Ly - 1].shape(0)));

         //add to value
         energy += Dot(LOSy,RO[Ly - 2]);

         //local expectation value
         auxvec[myID][(Ly-1)*Lx + col][1] = blas::dot(dim,LOU.data(),1,RO[Ly-2].data(),1);

         //3) then LSz with Sz

         //paste top environment on
         tmp5.clear();
         Contract(one,Environment::r[myID][col][Ly - 1],shape(1),Environment::Sz[myID](Ly-1,col),shape(1),zero,tmp5);

         //then bottom enviroment
         Contract(one,tmp5,shape(3),Environment::l[myID][col-1][Ly-1],shape(1),zero,tmp6);

         //move to a DArray<3> object
         RO[Ly - 2] = tmp6.reshape_clear(shape(Environment::r[myID][col][Ly - 1].shape(0),Environment::Sz[myID](Ly-1,col).shape(0),Environment::l[myID][col-1][Ly - 1].shape(0)));

         //add to value
         energy += Dot(LOSz,RO[Ly - 2]);

         //local expectation value
         auxvec[myID][(Ly-1)*Lx + col][2] = blas::dot(dim,LOU.data(),1,RO[Ly-2].data(),1);

      }

      // -- (3) -- || top row = Ly-1: again similar to overlap calculation

      //first construct the right renormalized operators

      //tmp comes out index (r,l)
      tmp4.clear();
      Contract(one,Environment::r[myID][Lx-2][Ly - 1],shape(1),Environment::l[myID][Lx-2][Ly - 1],shape(1),zero,tmp4);

      //reshape tmp to a 2-index array
      R[Ly - 2] = tmp4.reshape_clear(shape(Environment::r[myID][Lx-2][Ly - 1].shape(0),Environment::l[myID][Lx-2][Ly - 1].shape(0)));

      //now construct the rest
      for(int row = Ly - 2;row > 0;--row){

         tmp3.clear();
         Contract(one,Environment::r[myID][Lx-2][row],shape(2),R[row],shape(0),zero,tmp3);

         R[row - 1].clear();
         Contract(one,tmp3,shape(1,2),Environment::l[myID][Lx-2][row],shape(1,2),zero,R[row-1]);

      }

      //construct the left going operators on the first top site

      //first Sx

      //tmp comes out index (r,l)
      tmp5.clear();
      Contract(one,Environment::Sx[myID](0,Lx-1),shape(2),Environment::l[myID][Lx-2][0],shape(1),zero,tmp5);

      LSx = tmp5.reshape_clear(shape(Environment::Sx[myID](0,Lx-1).shape(3),Environment::l[myID][Lx-2][0].shape(2)));

      //then Sy

      //tmp5 comes out index (r,l)
      Contract(one,Environment::Sy[myID](0,Lx-1),shape(2),Environment::l[myID][Lx-2][0],shape(1),zero,tmp5);

      LSy = tmp5.reshape_clear(shape(Environment::Sy[myID](0,Lx-1).shape(3),Environment::l[myID][Lx-2][0].shape(2)));

      //then Sz 

      //tmp comes out index (r,l)
      Contract(one,Environment::Sz[myID](0,Ly-1),shape(2),Environment::l[myID][Ly-2][0],shape(1),zero,tmp5);

      LSz = tmp5.reshape_clear(shape(Environment::Sz[myID](0,Lx-1).shape(3),Environment::l[myID][Lx-2][0].shape(2)));

      //and finally unity
      Contract(one,Environment::r[myID][Lx-2][0],shape(1),Environment::l[myID][Lx-2][0],shape(1),zero,tmp4);

      LU = tmp4.reshape_clear(shape(Environment::r[myID][Lx-2][0].shape(2),Environment::l[myID][Lx-2][0].shape(2)));

      dim = R[0].size();

      //now contract x,y and z with R for local expectation values:
      auxvec[myID][Lx-1][0] = blas::dot(dim,LSx.data(),1,R[0].data(),1);
      auxvec[myID][Lx-1][1] = blas::dot(dim,LSy.data(),1,R[0].data(),1);
      auxvec[myID][Lx-1][2] = blas::dot(dim,LSz.data(),1,R[0].data(),1);

      //middle of the chain:
      for(int row = 1;row < Ly-1;++row){

         //first close down the x,y and z terms from the previous site for the energy

         //construct the right intermediate contraction (paste bottom to right)
         tmp3.clear();
         Contract(one,Environment::l[myID][Lx-2][row],shape(2),R[row],shape(1),zero,tmp3);

         // 1) paste Sx to the right
         int M = Environment::Sx[myID](row,Lx-1).shape(0);
         int N = tmp3.shape(0);
         int K = tmp3.shape(1) * tmp3.shape(2);

         blas::gemm(CblasRowMajor, CblasNoTrans, CblasTrans, M, N, K, one, Environment::Sx[myID](row,Lx-1).data(),K,tmp3.data(),K,zero,R[row-1].data(),N);

         //contract with left Sx
         energy += Dot(LSx,R[row - 1]);

         // 2) paste Sy to the right
         M = Environment::Sy[myID](row,Lx-1).shape(0);
         N = tmp3.shape(0);
         K = tmp3.shape(1) * tmp3.shape(2);

         blas::gemm(CblasRowMajor, CblasNoTrans, CblasTrans, M, N, K, one, Environment::Sy[myID](row,Lx-1).data(),K,tmp3.data(),K,zero,R[row-1].data(),N);

         //contract with left Sy
         energy += Dot(LSy,R[row - 1]);

         // 3) paste Sz to the right
         M = Environment::Sz[myID](row,Lx-1).shape(0);
         N = tmp3.shape(0);
         K = tmp3.shape(1) * tmp3.shape(2);

         blas::gemm(CblasRowMajor, CblasNoTrans, CblasTrans, M, N, K, one, Environment::Sz[myID](row,Lx-1).data(),K,tmp3.data(),K,zero,R[row-1].data(),N);

         //contract with left Sy
         energy += Dot(LSz,R[row - 1]);

         //construct left renormalized operators for next site: first paste bottom to Left unity
         tmp3.clear();
         Contract(one,LU,shape(1),Environment::l[myID][Lx-2][row],shape(0),zero,tmp3);

         // 1) construct new Sx left operator
         LSx.resize(Environment::Sx[myID](row,Lx-1).shape(3),Environment::l[myID][Lx-2][row].shape(2));

         M = Environment::Sx[myID](row,Lx-1).shape(3);
         N = tmp3.shape(2);
         K = tmp3.shape(0) * tmp3.shape(1);

         blas::gemm(CblasRowMajor, CblasTrans, CblasNoTrans, M, N, K, one, Environment::Sx[myID](row,Lx-1).data(),M,tmp3.data(),N,zero,LSx.data(),N);

         // 2) construct new Sy left operator
         LSy.resize(Environment::Sy[myID](row,Lx-1).shape(3),Environment::l[myID][Lx-2][row].shape(2));
         blas::gemm(CblasRowMajor, CblasTrans, CblasNoTrans, M, N, K, one, Environment::Sy[myID](row,Lx-1).data(),M,tmp3.data(),N,zero,LSy.data(),N);

         // 3) construct new Sz left operator
         LSz.resize(Environment::Sz[myID](row,Lx-1).shape(3),Environment::l[myID][Lx-2][row].shape(2));
         blas::gemm(CblasRowMajor, CblasTrans, CblasNoTrans, M, N, K, one, Environment::Sz[myID](row,Lx-1).data(),M,tmp3.data(),N,zero,LSz.data(),N);

         // 4) finally construct new unity on the left
         LU.resize(Environment::U[myID](row,Lx-1).shape(3),Environment::l[myID][Lx-2][row].shape(2));
         blas::gemm(CblasRowMajor, CblasTrans, CblasNoTrans, M, N, K, one, Environment::r[myID][Lx-2][row].data(),M,tmp3.data(),N,zero,LU.data(),N);

         //now contract x,y and z with R for local expectation values:
         dim = R[row].size();

         auxvec[myID][row*Lx + Lx-1][0] = blas::dot(dim,LSx.data(),1,R[row].data(),1);
         auxvec[myID][row*Lx + Lx-1][1] = blas::dot(dim,LSy.data(),1,R[row].data(),1);
         auxvec[myID][row*Lx + Lx-1][2] = blas::dot(dim,LSz.data(),1,R[row].data(),1);

      }

      //finally close down on last top site

      //first calculate overlap
      overlap = Dot(LU,R[Ly-2]);

      //1) Sx to close down LSx

      //tmp comes out index (r,l)
      tmp5.clear();
      Contract(one,Environment::Sx[myID](Ly-1,Lx-1),shape(2),Environment::l[myID][Lx-2][Ly - 1],shape(1),zero,tmp5);

      //reshape tmp to a 2-index array
      R[Ly - 2] = tmp5.reshape_clear(shape(Environment::Sx[myID](Ly-1,Lx-1).shape(0),Environment::l[myID][Lx-2][Ly - 1].shape(0)));

      //energy
      energy += Dot(LSx,R[Ly-2]);

      //local expectation value
      auxvec[myID][(Ly-1)*Lx + Lx-1][0] = blas::dot(dim,LU.data(),1,R[Ly-2].data(),1);

      //2) Sy to close down Ly
      Contract(one,Environment::Sy[myID](Ly-1,Lx-1),shape(2),Environment::l[myID][Lx-2][Ly - 1],shape(1),zero,tmp5);

      //reshape tmp to a 2-index array
      R[Ly - 2] = tmp5.reshape_clear(shape(Environment::Sy[myID](Ly-1,Lx-1).shape(0),Environment::l[myID][Lx-2][Ly - 1].shape(0)));

      //energy
      energy += Dot(LSy,R[Ly-2]);

      //local expectation value
      auxvec[myID][(Ly-1)*Lx + Lx-1][1] = blas::dot(dim,LU.data(),1,R[Ly-2].data(),1);

      //3) Sz to close down Lz

      //tmp comes out index (t,b)
      Contract(one,Environment::Sz[myID](Ly-1,Lx-1),shape(2),Environment::l[myID][Lx-2][Ly - 1],shape(1),zero,tmp5);

      //reshape tmp to a 2-index array
      R[Ly - 2] = tmp5.reshape_clear(shape(Environment::Sz[myID](Ly-1,Lx-1).shape(0),Environment::l[myID][Lx-2][Ly - 1].shape(0)));

      //energy
      energy += Dot(LSz,R[Ly-2]);

      //local expectation value
      auxvec[myID][(Ly-1)*Lx + Lx-1][2] = blas::dot(dim,LU.data(),1,R[Ly-2].data(),1);

      // #################################################################
      // ### ----               SET THE PROPERTIES                ---- ### 
      // #################################################################

      EL = energy/overlap;

      //set VL: Sx, Sy and Sz on every site: ready to construct auxiliary expectation values quickly!
      for(int k = 0;k < Trotter::n_trot;++k)
         for(int r = 0;r < 3;++r){

            VL[r*Trotter::n_trot + k] = auxvec[myID][0][r] * Trotter::gV()(k,0);


            for(int row = 0;row < Ly;++row)
               for(int col = 0;col < Lx;++col)
                  VL[r*Trotter::n_trot + k] += auxvec[myID][row*Lx + col][r] * Trotter::gV()(k,row*Lx + col);

            VL[r*Trotter::n_trot + k] /= overlap;

         }

   }

}

/**
 * Apply the propagator (a D=1 MPO) to the current state of the walker. The walker is changed when calling this function
 * @param P the propagator
 */
void Walker::propagate(const Propagator &P){

   complex<double> one(1.0,0.0);
   complex<double> zero(0.0,0.0);

   TArray<complex<double>,1> tmp(d);

   for(int i = 0;i < this->size();++i){

      blas::gemv(CblasRowMajor,CblasNoTrans, d, d, one, P[i].data(), d, (*this)[i].data(), 1, zero, tmp.data(), 1);

      (*this)[i] = std::move(tmp);

   }

}


/** 
 * normalize the walker state, will not physically effect the afqmc walk, just keeps the numbers from blowing up.
 */
void Walker::normalize(){

   for(int i = 0;i < this->size();++i){

      double nrm = blas::nrm2(d, (*this)[i].data(), 1);

      blas::scal(d, static_cast< complex<double> >(1.0/nrm), (*this)[i].data(), 1);

   }

}
