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
 * empty constructor:
 */
Walker::Walker() : std::vector< TArray<complex<double>,1> >( Lx * Ly ){ }

/**
 * construct a Walker object: initialize on AF state
 * @param n_trot_in number of trotter terms
 */
Walker::Walker(int n_trot_in) : std::vector< TArray<complex<double>,1> >( Lx * Ly ){

   this->n_trot = n_trot_in;

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

   VL.resize(3*n_trot);

}

/**
 * copy constructor
 * @param walker input Walker object to be copied
 */
Walker::Walker(const Walker &walker) : std::vector< TArray<complex<double>,1> >(walker) {

   this->weight = walker.gWeight();
   this->n_trot = walker.gn_trot();
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
 * @return the number of trotter terms
 */
int Walker::gn_trot() const {

   return n_trot;

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
complex<double> Walker::calc_properties(char option,const PEPS< complex<double> >& peps){

   // ---- || evaluate the expectation values in an MPO/MPS manner, first from bottom to top, then left to right || ----

   complex<double> energy(0.0,0.0);

   complex<double> one(1.0,0.0);
   complex<double> zero(0.0,0.0);

   if(option == 'H'){

      //calculate the single layer contractions first:
      Environment::U.fill('H',peps,*this);

      Environment::Sx.fill('H',peps,Sx,*this);
      Environment::Sy.fill('H',peps,Sy,*this);
      Environment::Sz.fill('H',peps,Sz,*this);

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
      Contract(one,Environment::t[0][Lx - 1],shape(1),Environment::b[0][Lx - 1],shape(1),zero,tmp4);

      //reshape tmp to a 2-index array
      R[Lx - 2] = tmp4.reshape_clear(shape(Environment::t[0][Lx - 1].shape(0),Environment::b[0][Lx - 1].shape(0)));

      //now construct the rest
      for(int col = Lx - 2;col > 0;--col){

         tmp3.clear();
         Contract(one,Environment::t[0][col],shape(2),R[col],shape(0),zero,tmp3);

         Contract(one,tmp3,shape(1,2),Environment::b[0][col],shape(1,2),zero,R[col-1]);

      }

      //4 left going operators: Sx, Sy, Sz, and 1
      ZArray<2> LSx;
      ZArray<2> LSy;
      ZArray<2> LSz;
      ZArray<2> LU;

      TArray<complex<double>,5> tmp5;

      //tmp comes out index (t,b)
      Contract(one,Environment::t[0][0],shape(1),Environment::Sx(0,0),shape(1),zero,tmp5);

      LSx = tmp5.reshape_clear(shape(Environment::t[0][0].shape(2),Environment::Sx(0,0).shape(3)));

      //then Sy
      Contract(one,Environment::t[0][0],shape(1),Environment::Sy(0,0),shape(1),zero,tmp5);

      LSy = tmp5.reshape_clear(shape(Environment::t[0][0].shape(2),Environment::Sy(0,0).shape(3)));

      //then Sz
      Contract(one,Environment::t[0][0],shape(1),Environment::Sz(0,0),shape(1),zero,tmp5);

      LSz = tmp5.reshape_clear(shape(Environment::t[0][0].shape(2),Environment::Sz(0,0).shape(3)));

      //finally unity
      Contract(one,Environment::t[0][0],shape(1),Environment::U(0,0),shape(1),zero,tmp5);

      LU = tmp5.reshape_clear(shape(Environment::t[0][0].shape(2),Environment::Sz(0,0).shape(3)));

      int dim = R[0].size();

      //now contract x,y and z with R for local expectation values:
      auxvec[0][0] = blas::dot(dim,LSx.data(),1,R[0].data(),1);
      auxvec[0][1] = blas::dot(dim,LSy.data(),1,R[0].data(),1);
      auxvec[0][2] = blas::dot(dim,LSz.data(),1,R[0].data(),1);

      //now for the middle terms
      for(int col = 1;col < Lx - 1;++col){

         //first close down the x,y and z terms from the previous site for the energy

         //construct the right intermediate contraction (paste top to right)
         tmp3.clear();
         Contract(one,Environment::t[0][col],shape(2),R[col],shape(0),zero,tmp3);

         // 1) paste Sx to the right
         int M = tmp3.shape(0);
         int N = Environment::Sx(0,col).shape(0);
         int K = tmp3.shape(1) * tmp3.shape(2);

         blas::gemm(CblasRowMajor, CblasNoTrans, CblasTrans, M, N, K, one, tmp3.data(),K,Environment::Sx(0,col).data(),K,zero,R[col-1].data(),N);

         //contract with left Sx
         energy += Dot(LSx,R[col - 1]);

         // 2) then paste Sy to the right
         blas::gemm(CblasRowMajor, CblasNoTrans, CblasTrans, M, N, K, one, tmp3.data(),K,Environment::Sy(0,col).data(),K,zero,R[col-1].data(),N);

         //contract with left Sy
         energy += Dot(LSy,R[col - 1]);

         // 3) finally Sz
         blas::gemm(CblasRowMajor, CblasNoTrans, CblasTrans, M, N, K, one, tmp3.data(),K,Environment::Sz(0,col).data(),K,zero,R[col-1].data(),N);

         //contract with left Sz
         energy += Dot(LSz,R[col - 1]);

         //construct left renormalized operators for next site: first paste top to Left unity
         tmp3.clear();
         Contract(one,LU,shape(0),Environment::t[0][col],shape(0),zero,tmp3);

         // 1) construct new Sx left operator
         LSx.resize(Environment::t[0][col].shape(2),Environment::Sx(0,col).shape(3));

         M = tmp3.shape(2);
         N = Environment::Sx(0,col).shape(3);
         K = tmp3.shape(0) * tmp3.shape(1);

         blas::gemm(CblasRowMajor, CblasTrans, CblasNoTrans, M, N, K, one, tmp3.data(),M,Environment::Sx(0,col).data(),N,zero,LSx.data(),N);

         // 2) construct new Sy left operator
         LSy.resize(Environment::t[0][col].shape(2),Environment::Sy(0,col).shape(3));
         blas::gemm(CblasRowMajor, CblasTrans, CblasNoTrans, M, N, K, one, tmp3.data(),M,Environment::Sy(0,col).data(),N,zero,LSx.data(),N);

         // 3) construct new Sz left operator
         LSz.resize(Environment::t[0][col].shape(2),Environment::Sz(0,col).shape(3));
         blas::gemm(CblasRowMajor, CblasTrans, CblasNoTrans, M, N, K, one, tmp3.data(),M,Environment::Sz(0,col).data(),N,zero,LSx.data(),N);

         //now contract x,y and z with R for local expectation values:
         dim = R[col].size();

         auxvec[col][0] = blas::dot(dim,LSx.data(),1,R[col].data(),1);
         auxvec[col][1] = blas::dot(dim,LSy.data(),1,R[col].data(),1);
         auxvec[col][2] = blas::dot(dim,LSz.data(),1,R[col].data(),1);

         // 4) finally construct new unity on the left
         LU.resize(Environment::t[0][col].shape(2),Environment::U(0,col).shape(3));
         blas::gemm(CblasRowMajor, CblasTrans, CblasNoTrans, M, N, K, one, tmp3.data(),M,Environment::U(0,col).data(),N,zero,LSx.data(),N);

      }
      
      //last site of bottom row:close down the left x,y and z

      //1) Sx to close down LSx
      Contract(one,Environment::t[0][Lx-1],shape(1),Environment::Sx(0,Lx-1),shape(1),zero,tmp5);

      R[Lx-2] = tmp5.reshape_clear(shape(Environment::t[0][Lx-1].shape(0),Environment::Sx(0,Lx-1).shape(0)));

      energy += Dot(LSx,R[Lx-2]);

      dim = R[Lx-2].size();
      auxvec[Lx-1][0] = blas::dot(dim,LSx.data(),1,R[Lx-2].data(),1);

      //2) Sy to close down Ly
      Contract(one,Environment::t[0][Lx-1],shape(1),Environment::Sy(0,Lx-1),shape(1),zero,tmp5);

      R[Lx-2] = tmp5.reshape_clear(shape(Environment::t[0][Lx-1].shape(0),Environment::Sy(0,Lx-1).shape(0)));

      energy += Dot(LSy,R[Lx-2]);

      auxvec[Lx-1][0] = blas::dot(dim,LSy.data(),1,R[Lx-2].data(),1);

      //3) Sz to close down Lz
      Contract(one,Environment::t[0][Lx-1],shape(1),Environment::Sz(0,Lx-1),shape(1),zero,tmp5);

      R[Lx-2] = tmp5.reshape_clear(shape(Environment::t[0][Lx-1].shape(0),Environment::Sz(0,Lx-1).shape(0)));

      energy += Dot(LSz,R[Lx-2]);

      auxvec[Lx-1][0] = blas::dot(dim,LSz.data(),1,R[Lx-2].data(),1);
/*
      // -- (2) -- now move from bottom to top calculating everything like an MPO/MPS expectation value

      //Right renormalized operators
      vector< DArray<3> > RO(Lx - 2);

      //4 left renormalized operators needed
      DArray<3> LOp;
      DArray<3> LOm;
      DArray<3> LOz;
      DArray<3> LOu;

      //double layer objects
      DArray<4> dlop;
      DArray<4> dlom;
      DArray<4> dloz;
      DArray<4> dlou;

      for(int row = 1;row < Ly - 1;++row){

      //first create right renormalized operator

      //first site make double layer object from peps
      Environment::construct_double_layer('H',peps(row,Lx-1),dlou);

      //paste top environment on
      DArray<5> tmp5;
      Contract(1.0,Environment::t[row][Lx - 1],shape(1),dlou,shape(1),0.0,tmp5);

      //then bottom enviroment
      DArray<6> tmp6;
      Contract(1.0,tmp5,shape(3),Environment::b[row-1][Lx-1],shape(1),0.0,tmp6);

      //move to a DArray<3> object
      RO[Lx - 3] = tmp6.reshape_clear(shape(Environment::t[row][Lx - 1].shape(0),dlou.shape(0),Environment::b[row-1][Lx - 1].shape(0)));

      DArray<4> I4;
      DArray<4> I4bis;

      //now construct the middle operators
      for(int col = Lx-2;col > 1;--col){

         I4.clear();
         Contract(1.0,Environment::t[row][col],shape(2),RO[col-1],shape(0),0.0,I4);

         enum {i,j,k,o,m,n};

         Environment::construct_double_layer('H',peps(row,col),dlou);

         I4bis.clear();
         Contract(1.0,I4,shape(i,j,k,o),dlou,shape(m,j,n,k),0.0,I4bis,shape(i,m,n,o));

         RO[col-2].clear();
         Contract(1.0,I4bis,shape(2,3),Environment::b[row-1][col],shape(1,2),0.0,RO[col-2]);

      }

      // --- now move from left to right to get the expecation value of the interactions ---
      // --- First construct the left going operators for the first site -----

      // 1) S+ -- make double layer object from peps with Sp
      Environment::construct_double_layer('H',peps(row,0),Sp,dlop);

      //paste top environment on
      tmp5.clear();
      Contract(1.0,Environment::t[row][0],shape(1),dlop,shape(1),0.0,tmp5);

      //then bottom enviroment
      tmp6.clear();
      Contract(1.0,tmp5,shape(3),Environment::b[row-1][0],shape(1),0.0,tmp6);

      //move to a DArray<3> object: order (top-env,peps-row,bottom-env)
      LOp = tmp6.reshape_clear(shape(Environment::t[row][0].shape(2),dlop.shape(3),Environment::b[row-1][0].shape(2)));

      // 2) S- -- make double layer object from peps with Sm
      Environment::construct_double_layer('H',peps(row,0),Sm,dlom);

      //paste top environment on
      tmp5.clear();
      Contract(1.0,Environment::t[row][0],shape(1),dlom,shape(1),0.0,tmp5);

      //then bottom enviroment
      tmp6.clear();
      Contract(1.0,tmp5,shape(3),Environment::b[row-1][0],shape(1),0.0,tmp6);

      //move to a DArray<3> object: order (top-env,peps-row,bottom-env)
      LOm = tmp6.reshape_clear(shape(Environment::t[row][0].shape(2),dlom.shape(3),Environment::b[row-1][0].shape(2)));

      // 3) Sz -- make double layer object from peps with Sz
      Environment::construct_double_layer('H',peps(row,0),Sz,dloz);

      //paste top environment on
      tmp5.clear();
      Contract(1.0,Environment::t[row][0],shape(1),dloz,shape(1),0.0,tmp5);

      //then bottom enviroment
      tmp6.clear();
      Contract(1.0,tmp5,shape(3),Environment::b[row-1][0],shape(1),0.0,tmp6);

      //move to a DArray<3> object: order (top-env,peps-row,bottom-env)
      LOz = tmp6.reshape_clear(shape(Environment::t[row][0].shape(2),dlom.shape(3),Environment::b[row-1][0].shape(2)));

      // 4) 1 -- finally construct left renormalized operator with unity
      Environment::construct_double_layer('H',peps(row,0),dlou);

      //paste top environment on
      tmp5.clear();
      Contract(1.0,Environment::t[row][0],shape(1),dlou,shape(1),0.0,tmp5);

      //then bottom enviroment
      tmp6.clear();
      Contract(1.0,tmp5,shape(3),Environment::b[row-1][0],shape(1),0.0,tmp6);

      //move to a DArray<3> object: order (top-env,peps-row,bottom-env)
      LOu = tmp6.reshape_clear(shape(Environment::t[row][0].shape(2),dlou.shape(3),Environment::b[row-1][0].shape(2)));

      // --- now for the middle sites, close down the operators on the left and construct new ones --- 
      for(int col = 1;col < Lx - 1;++col){

         //first add top to the right side, put it in I4
         I4.clear();
         Contract(1.0,Environment::t[row][col],shape(2),RO[col-1],shape(0),0.0,I4);

         enum {i,j,k,o,m,n};

         //1) close down LOp with Sm
         Environment::construct_double_layer('H',peps(row,col),Sm,dlom);

         I4bis.clear();
         Contract(1.0,I4,shape(i,j,k,o),dlom,shape(m,j,n,k),0.0,I4bis,shape(i,m,n,o));

         RO[col-1].clear();
         Contract(1.0,I4bis,shape(2,3),Environment::b[row-1][col],shape(1,2),0.0,RO[col-1]);

         //expectation value:
         val += 0.5 * Dot(LOp,RO[col-1]);

         //2) close down LOm with Sp
         Environment::construct_double_layer('H',peps(row,col),Sp,dlop);

         I4bis.clear();
         Contract(1.0,I4,shape(i,j,k,o),dlop,shape(m,j,n,k),0.0,I4bis,shape(i,m,n,o));

         RO[col-1].clear();
         Contract(1.0,I4bis,shape(2,3),Environment::b[row-1][col],shape(1,2),0.0,RO[col-1]);

         //expectation value:
         val += 0.5 * Dot(LOm,RO[col-1]);

         //3) finally close down LOz with Sz
         Environment::construct_double_layer('H',peps(row,col),Sz,dloz);

         I4bis.clear();
         Contract(1.0,I4,shape(i,j,k,o),dloz,shape(m,j,n,k),0.0,I4bis,shape(i,m,n,o));

         RO[col-1].clear();
         Contract(1.0,I4bis,shape(2,3),Environment::b[row-1][col],shape(1,2),0.0,RO[col-1]);

         //expectation value:
         val += Dot(LOz,RO[col-1]);

         // now construct the new left going renormalized operators
         //first attach top to left unity
         I4.clear();
         Contract(1.0,Environment::t[row][col],shape(0),LOu,shape(0),0.0,I4);

         // 1) construct left Sp operator
         I4bis.clear();
         Contract(1.0,I4,shape(i,j,k,o),dlop,shape(k,i,m,n),0.0,I4bis,shape(j,n,o,m));

         LOp.clear();
         Contract(1.0,I4bis,shape(2,3),Environment::b[row-1][col],shape(0,1),0.0,LOp);

         // 2) construct left Sm operator
         I4bis.clear();
         Contract(1.0,I4,shape(i,j,k,o),dlom,shape(k,i,m,n),0.0,I4bis,shape(j,n,o,m));

         LOm.clear();
         Contract(1.0,I4bis,shape(2,3),Environment::b[row-1][col],shape(0,1),0.0,LOm);

         // 3) construct left Sz operator
         I4bis.clear();
         Contract(1.0,I4,shape(i,j,k,o),dloz,shape(k,i,m,n),0.0,I4bis,shape(j,n,o,m));

         LOz.clear();
         Contract(1.0,I4bis,shape(2,3),Environment::b[row-1][col],shape(0,1),0.0,LOz);

         // 4) finally construct new left unity
         Environment::construct_double_layer('H',peps(row,col),dlou);

         I4bis.clear();
         Contract(1.0,I4,shape(i,j,k,o),dlou,shape(k,i,m,n),0.0,I4bis,shape(j,n,o,m));

         LOu.clear();
         Contract(1.0,I4bis,shape(2,3),Environment::b[row-1][col],shape(0,1),0.0,LOu);

      }

      //last site on the right: close down on the incomings

      //1) first Lp with Sm
      Environment::construct_double_layer('H',peps(row,Lx-1),Sm,dlom);

      //paste top environment on
      tmp5.clear();
      Contract(1.0,Environment::t[row][Lx - 1],shape(1),dlom,shape(1),0.0,tmp5);

      //then bottom enviroment
      tmp6.clear();
      Contract(1.0,tmp5,shape(3),Environment::b[row-1][Lx-1],shape(1),0.0,tmp6);

      //move to a DArray<3> object
      RO[Lx - 3] = tmp6.reshape_clear(shape(Environment::t[row][Lx - 1].shape(0),dlom.shape(0),Environment::b[row-1][Lx - 1].shape(0)));

      //add to value
      val += 0.5 * Dot(LOp,RO[Lx - 3]);

      //2) then Lm with Sp
      Environment::construct_double_layer('H',peps(row,Lx-1),Sp,dlop);

      //paste top environment on
      tmp5.clear();
      Contract(1.0,Environment::t[row][Lx - 1],shape(1),dlop,shape(1),0.0,tmp5);

      //then bottom enviroment
      tmp6.clear();
      Contract(1.0,tmp5,shape(3),Environment::b[row-1][Lx-1],shape(1),0.0,tmp6);

      //move to a DArray<3> object
      RO[Lx - 3] = tmp6.reshape_clear(shape(Environment::t[row][Lx - 1].shape(0),dlop.shape(0),Environment::b[row-1][Lx - 1].shape(0)));

      //add to value
      val += 0.5 * Dot(LOm,RO[Lx - 3]);

      //3) then Lz with Sz
      Environment::construct_double_layer('H',peps(row,Lx-1),Sz,dloz);

      //paste top environment on
      tmp5.clear();
      Contract(1.0,Environment::t[row][Lx - 1],shape(1),dloz,shape(1),0.0,tmp5);

      //then bottom enviroment
      tmp6.clear();
      Contract(1.0,tmp5,shape(3),Environment::b[row-1][Lx-1],shape(1),0.0,tmp6);

      //move to a DArray<3> object
      RO[Lx - 3] = tmp6.reshape_clear(shape(Environment::t[row][Lx - 1].shape(0),dloz.shape(0),Environment::b[row-1][Lx - 1].shape(0)));

      //add to value
      val += Dot(LOz,RO[Lx - 3]);

}

// -- (3) -- || top row = Ly-1: again similar to overlap calculation

//first construct the right renormalized operators

//tmp comes out index (t,b)
tmp.clear();
Contract(1.0,Environment::t[Ly-2][Lx - 1],shape(1),Environment::b[Ly-2][Lx - 1],shape(1),0.0,tmp);

//reshape tmp to a 2-index array
R[Lx - 3] = tmp.reshape_clear(shape(Environment::t[Ly-2][Lx - 1].shape(0),Environment::b[Ly-2][Lx - 1].shape(0)));

//now construct the rest
for(int col = Lx - 2;col > 1;--col){

   I.clear();
   Contract(1.0,Environment::t[Ly-2][col],shape(2),R[col-1],shape(0),0.0,I);

   R[col-2].clear();
   Contract(1.0,I,shape(1,2),Environment::b[Ly-2][col],shape(1,2),0.0,R[col-2]);

}

//construct the left going operators on the first top site

//first S+
Environment::construct_double_layer('H',peps(Ly-1,0),Sp,dlsp);

//tmp comes out index (t,b)
Contract(1.0,dlsp,shape(1),Environment::b[Ly-2][0],shape(1),0.0,tmp);

Lp = tmp.reshape_clear(shape(dlsp.shape(2),Environment::b[Ly-2][0].shape(2)));

//then S-
Environment::construct_double_layer('H',peps(Ly-1,0),Sm,dlsm);

//tmp comes out index (t,b)
Contract(1.0,dlsm,shape(1),Environment::b[Ly-2][0],shape(1),0.0,tmp);

Lm = tmp.reshape_clear(shape(dlsm.shape(2),Environment::b[Ly-2][0].shape(2)));

//then Sz 
Environment::construct_double_layer('H',peps(Ly-1,0),Sz,dlsz);

//tmp comes out index (t,b)
Contract(1.0,dlsz,shape(1),Environment::b[Ly-2][0],shape(1),0.0,tmp);

Lz = tmp.reshape_clear(shape(dlsz.shape(2),Environment::b[Ly-2][0].shape(2)));

//and finally unity
Contract(1.0,Environment::t[Ly-2][0],shape(1),Environment::b[Ly-2][0],shape(1),0.0,tmp);

Lu = tmp.reshape_clear(shape(Environment::t[Ly-2][0].shape(2),Environment::b[Ly-2][0].shape(2)));

//middle of the chain:
for(int col = 1;col < Lx-1;++col){

   //first close down the +,- and z terms from the previous site

   //construct the right intermediate contraction (paste bottom to right)
   I.clear();
   Contract(1.0,Environment::b[Ly-2][col],shape(2),R[col - 1],shape(1),0.0,I);

   // 1) construct Sm double layer
   Environment::construct_double_layer('H',peps(Ly-1,col),Sm,dlsm);

   R[col-1].clear();
   Contract(1.0,dlsm,shape(1,2),I,shape(1,2),0.0,R[col - 1]);

   //contract with left S+
   val += 0.5 * Dot(Lp,R[col - 1]);

   // 2) construct Sp double layer
   Environment::construct_double_layer('H',peps(Ly-1,col),Sp,dlsp);

   R[col-1].clear();
   Contract(1.0,dlsp,shape(1,2),I,shape(1,2),0.0,R[col - 1]);

   //contract with left S-
   val += 0.5 * Dot(Lm,R[col - 1]);

   // 3) construct Sz double layer
   Environment::construct_double_layer('H',peps(Ly-1,col),Sz,dlsz);

   R[col-1].clear();
   Contract(1.0,dlsz,shape(1,2),I,shape(1,2),0.0,R[col - 1]);

   //contract with left Sz
   val += Dot(Lz,R[col - 1]);

   //construct left renormalized operators for next site: first paste bottom to Left unity
   I.clear();
   Contract(1.0,Lu,shape(1),Environment::b[Ly-2][col],shape(0),0.0,I);

   // 1) construct new Sm left operator
   Lm.clear();
   Contract(1.0,dlsm,shape(0,1),I,shape(0,1),0.0,Lm);

   // 2) construct new Sp left operator
   Lp.clear();
   Contract(1.0,dlsp,shape(0,1),I,shape(0,1),0.0,Lp);

   // 3) construct new Sz left operator
   Lz.clear();
   Contract(1.0,dlsz,shape(0,1),I,shape(0,1),0.0,Lz);

   // 4) finally construct new unity on the left
   Lu.clear();
   Contract(1.0,Environment::t[Ly-2][col],shape(0,1),I,shape(0,1),0.0,Lu);

}

//finally close down on last top site

//1) Sm to close down Lp
Environment::construct_double_layer('H',peps(Ly-1,Lx-1),Sm,dlsm);

//tmp comes out index (t,b)
tmp.clear();
Contract(1.0,dlsm,shape(1),Environment::b[Ly-2][Lx - 1],shape(1),0.0,tmp);

//reshape tmp to a 2-index array
R[Lx - 3] = tmp.reshape_clear(shape(dlsm.shape(0),Environment::b[Ly-2][Lx - 1].shape(0)));

val += 0.5 * Dot(Lp,R[Lx-3]);

//2) Sp to close down Lm
Environment::construct_double_layer('H',peps(Ly-1,Lx-1),Sp,dlsp);

//tmp comes out index (t,b)
tmp.clear();
Contract(1.0,dlsp,shape(1),Environment::b[Ly-2][Lx - 1],shape(1),0.0,tmp);

//reshape tmp to a 2-index array
R[Lx - 3] = tmp.reshape_clear(shape(dlsp.shape(0),Environment::b[Ly-2][Lx - 1].shape(0)));

val += 0.5 * Dot(Lm,R[Lx-3]);

//3) Sz to close down Lz
Environment::construct_double_layer('H',peps(Ly-1,Lx-1),Sz,dlsz);

//tmp comes out index (t,b)
tmp.clear();
Contract(1.0,dlsz,shape(1),Environment::b[Ly-2][Lx - 1],shape(1),0.0,tmp);

//reshape tmp to a 2-index array
R[Lx - 3] = tmp.reshape_clear(shape(dlsz.shape(0),Environment::b[Ly-2][Lx - 1].shape(0)));

val += Dot(Lz,R[Lx-3]);
*/
}
else{//VERTICAL: Left to right
   /*
   // #################################################################
   // ### ---- from left to right: contract in mps/mpo fashion ---- ### 
   // #################################################################

   // -- (1) -- || left column: similar to overlap calculation

   //first construct the right renormalized operators

   //first the rightmost operator

   //tmp comes out index (r,l)
   Contract(1.0,Environment::r[0][Ly - 1],shape(1),Environment::l[0][Ly - 1],shape(1),0.0,tmp);

   //reshape tmp to a 2-index array
   R[Ly - 3] = tmp.reshape_clear(shape(Environment::r[0][Ly - 1].shape(0),Environment::l[0][Ly - 1].shape(0)));

   //now construct the rest
   for(int row = Ly - 2;row > 1;--row){

   I.clear();
   Contract(1.0,Environment::r[0][row],shape(2),R[row-1],shape(0),0.0,I);

   R[row-2].clear();
   Contract(1.0,I,shape(1,2),Environment::l[0][row],shape(1,2),0.0,R[row-2]);

   }

   //4 left going operators: S+, S-, Sz, and 1

   //first S+
   Environment::construct_double_layer('V',peps(0,0),Sp,dlsp);

   //tmp comes out index (r,l)
   Contract(1.0,Environment::r[0][0],shape(1),dlsp,shape(1),0.0,tmp);

   Lp = tmp.reshape_clear(shape(Environment::r[0][0].shape(2),dlsp.shape(2)));

   //then S-
   Environment::construct_double_layer('V',peps(0,0),Sm,dlsm);

   //tmp comes out index (r,l)
   Contract(1.0,Environment::r[0][0],shape(1),dlsm,shape(1),0.0,tmp);

   Lm = tmp.reshape_clear(shape(Environment::r[0][0].shape(2),dlsm.shape(2)));

   //then Sz 
   Environment::construct_double_layer('V',peps(0,0),Sz,dlsz);

   //tmp comes out index (r,l)
   Contract(1.0,Environment::r[0][0],shape(1),dlsz,shape(1),0.0,tmp);

   Lz = tmp.reshape_clear(shape(Environment::r[0][0].shape(2),dlsz.shape(2)));

   //and finally unity
   Contract(1.0,Environment::r[0][0],shape(1),Environment::l[0][0],shape(1),0.0,tmp);

   Lu = tmp.reshape_clear(shape(Environment::r[0][0].shape(2),Environment::l[0][0].shape(2)));

   //now for the middle terms
   for(int row = 1;row < Ly - 1;++row){

   //first close down the +,- and z terms from the previous site

   //construct the right intermediate contraction (paste 'right' to R)
   I.clear();
   Contract(1.0,Environment::r[0][row],shape(2),R[row - 1],shape(0),0.0,I);

   // 1) construct Sm double layer
   Environment::construct_double_layer('V',peps(row,0),Sm,dlsm);

   R[row-1].clear();
   Contract(1.0,I,shape(1,2),dlsm,shape(1,2),0.0,R[row - 1]);

   //contract with left S+
   val += 0.5 * Dot(Lp,R[row - 1]);

   // 2) then construct Sp double layer
   Environment::construct_double_layer('V',peps(row,0),Sp,dlsp);

   R[row-1].clear();
   Contract(1.0,I,shape(1,2),dlsp,shape(1,2),0.0,R[row - 1]);

   //contract with left S-
   val += 0.5 * Dot(Lm,R[row - 1]);

   // 3) then construct Sz double layer
   Environment::construct_double_layer('V',peps(row,0),Sz,dlsz);

   R[row-1].clear();
   Contract(1.0,I,shape(1,2),dlsz,shape(1,2),0.0,R[row - 1]);

   //contract with left Sz
   val += Dot(Lz,R[row - 1]);

   //construct left renormalized operators for next site: first paste top to Left unity
   I.clear();
   Contract(1.0,Lu,shape(0),Environment::r[0][row],shape(0),0.0,I);

   // 1) construct new Sm left operator
   Lm.clear();
   Contract(1.0,I,shape(0,1),dlsm,shape(0,1),0.0,Lm);

   // 2) construct new Sp left operator
   Lp.clear();
   Contract(1.0,I,shape(0,1),dlsp,shape(0,1),0.0,Lp);

   // 3) construct new Sz left operator
   Lz.clear();
   Contract(1.0,I,shape(0,1),dlsz,shape(0,1),0.0,Lz);

   // 4) finally construct new unity on the left
   Lu.clear();
   Contract(1.0,I,shape(0,1),Environment::l[0][row],shape(0,1),0.0,Lu);

}

//last site of left column: close down the left +,- and z

//1) Sm to close down Lp
Environment::construct_double_layer('V',peps(Ly-1,0),Sm,dlsm);

//tmp comes out index (r,l)
tmp.clear();
Contract(1.0,Environment::r[0][Ly - 1],shape(1),dlsm,shape(1),0.0,tmp);

//reshape tmp to a 2-index array
R[Ly - 3] = tmp.reshape_clear(shape(Environment::r[0][Ly - 1].shape(0),dlsm.shape(0)));

val += 0.5 * Dot(Lp,R[Ly-3]);

//2) Sp to close down Lm
Environment::construct_double_layer('V',peps(Ly-1,0),Sp,dlsp);

//tmp comes out index (r,l)
tmp.clear();
Contract(1.0,Environment::r[0][Ly - 1],shape(1),dlsp,shape(1),0.0,tmp);

//reshape tmp to a 2-index array
R[Ly - 3] = tmp.reshape_clear(shape(Environment::r[0][Ly - 1].shape(0),dlsp.shape(0)));

val += 0.5 * Dot(Lm,R[Ly-3]);

//3) Sz to close down Lz
Environment::construct_double_layer('V',peps(Ly-1,0),Sz,dlsz);

//tmp comes out index (t,b)
tmp.clear();
Contract(1.0,Environment::r[0][Ly - 1],shape(1),dlsz,shape(1),0.0,tmp);

//reshape tmp to a 2-index array
R[Ly - 3] = tmp.reshape_clear(shape(Environment::r[0][Ly - 1].shape(0),dlsz.shape(0)));

val += Dot(Lz,R[Ly-3]);

// -- (2) -- now move from left to right calculating everything like an MPO/MPS expectation value
for(int col = 1;col < Lx - 1;++col){

   //first create right renormalized operator

   //first site make double layer object from peps
   Environment::construct_double_layer('V',peps(Ly-1,col),dlou);

   //paste right environment on
   DArray<5> tmp5;
   Contract(1.0,Environment::r[col][Ly - 1],shape(1),dlou,shape(1),0.0,tmp5);

   //then left enviroment
   DArray<6> tmp6;
   Contract(1.0,tmp5,shape(3),Environment::l[col-1][Ly-1],shape(1),0.0,tmp6);

   //move to a DArray<3> object
   RO[Lx - 3] = tmp6.reshape_clear(shape(Environment::r[col][Ly - 1].shape(0),dlou.shape(0),Environment::l[col-1][Ly - 1].shape(0)));

   DArray<4> I4;
   DArray<4> I4bis;

   //now construct the middle operators
   for(int row = Ly-2;row > 1;--row){

      I4.clear();
      Contract(1.0,Environment::r[col][row],shape(2),RO[row-1],shape(0),0.0,I4);

      enum {i,j,k,o,m,n};

      Environment::construct_double_layer('V',peps(row,col),dlou);

      I4bis.clear();
      Contract(1.0,I4,shape(i,j,k,o),dlou,shape(m,j,n,k),0.0,I4bis,shape(i,m,n,o));

      RO[row-2].clear();
      Contract(1.0,I4bis,shape(2,3),Environment::l[col-1][row],shape(1,2),0.0,RO[row-2]);

   }

   // --- now move from left to right to get the expecation value of the interactions ---
   // --- First construct the left going operators for the first site -----

   // 1) S+ -- make double layer object from peps with Sp
   Environment::construct_double_layer('V',peps(0,col),Sp,dlop);

   //paste right environment on
   tmp5.clear();
   Contract(1.0,Environment::r[col][0],shape(1),dlop,shape(1),0.0,tmp5);

   //then left enviroment
   tmp6.clear();
   Contract(1.0,tmp5,shape(3),Environment::l[col-1][0],shape(1),0.0,tmp6);

   //move to a DArray<3> object: order (right-env,peps-col,left-env)
   LOp = tmp6.reshape_clear(shape(Environment::r[col][0].shape(2),dlop.shape(3),Environment::l[col-1][0].shape(2)));

   // 2) S- -- make double layer object from peps with Sm
   Environment::construct_double_layer('V',peps(0,col),Sm,dlom);

   //paste right environment on
   tmp5.clear();
   Contract(1.0,Environment::r[col][0],shape(1),dlom,shape(1),0.0,tmp5);

   //then left enviroment
   tmp6.clear();
   Contract(1.0,tmp5,shape(3),Environment::l[col-1][0],shape(1),0.0,tmp6);

   //move to a DArray<3> object: 
   LOm = tmp6.reshape_clear(shape(Environment::r[col][0].shape(2),dlom.shape(3),Environment::l[col-1][0].shape(2)));

   // 3) Sz -- make double layer object from peps with Sz
   Environment::construct_double_layer('V',peps(0,col),Sz,dloz);

   //paste right environment on
   tmp5.clear();
   Contract(1.0,Environment::r[col][0],shape(1),dloz,shape(1),0.0,tmp5);

   //then bottom enviroment
   tmp6.clear();
   Contract(1.0,tmp5,shape(3),Environment::l[col-1][0],shape(1),0.0,tmp6);

   //move to a DArray<3> object: order (top-env,peps-row,bottom-env)
   LOz = tmp6.reshape_clear(shape(Environment::r[col][0].shape(2),dlom.shape(3),Environment::l[col-1][0].shape(2)));

   // 4) 1 -- finally construct left renormalized operator with unity
   Environment::construct_double_layer('V',peps(0,col),dlou);

   //paste right environment on
   tmp5.clear();
   Contract(1.0,Environment::r[col][0],shape(1),dlou,shape(1),0.0,tmp5);

   //then left enviroment
   tmp6.clear();
   Contract(1.0,tmp5,shape(3),Environment::l[col-1][0],shape(1),0.0,tmp6);

   //move to a DArray<3> object: 
   LOu = tmp6.reshape_clear(shape(Environment::r[col][0].shape(2),dlou.shape(3),Environment::l[col-1][0].shape(2)));

   // --- now for the middle sites, close down the operators on the left and construct new ones --- 
   for(int row = 1;row < Ly - 1;++row){

      //first add right to the right side, put it in I4
      I4.clear();
      Contract(1.0,Environment::r[col][row],shape(2),RO[row-1],shape(0),0.0,I4);

      enum {i,j,k,o,m,n};

      //1) close down LOp with Sm
      Environment::construct_double_layer('V',peps(row,col),Sm,dlom);

      I4bis.clear();
      Contract(1.0,I4,shape(i,j,k,o),dlom,shape(m,j,n,k),0.0,I4bis,shape(i,m,n,o));

      RO[row-1].clear();
      Contract(1.0,I4bis,shape(2,3),Environment::l[col-1][row],shape(1,2),0.0,RO[row-1]);

      //expectation value:
      val += 0.5 * Dot(LOp,RO[row-1]);

      //2) close down LOm with Sp
      Environment::construct_double_layer('V',peps(row,col),Sp,dlop);

      I4bis.clear();
      Contract(1.0,I4,shape(i,j,k,o),dlop,shape(m,j,n,k),0.0,I4bis,shape(i,m,n,o));

      RO[row-1].clear();
      Contract(1.0,I4bis,shape(2,3),Environment::l[col-1][row],shape(1,2),0.0,RO[row-1]);

      //expectation value:
      val += 0.5 * Dot(LOm,RO[row-1]);

      //3) finally close down LOz with Sz
      Environment::construct_double_layer('V',peps(row,col),Sz,dloz);

      I4bis.clear();
      Contract(1.0,I4,shape(i,j,k,o),dloz,shape(m,j,n,k),0.0,I4bis,shape(i,m,n,o));

      RO[row-1].clear();
      Contract(1.0,I4bis,shape(2,3),Environment::l[col-1][row],shape(1,2),0.0,RO[row-1]);

      //expectation value:
      val += Dot(LOz,RO[row-1]);

      // now construct the new left going renormalized operators
      //first attach top to left unity
      I4.clear();
      Contract(1.0,Environment::r[col][row],shape(0),LOu,shape(0),0.0,I4);

      // 1) construct left Sp operator
      I4bis.clear();
      Contract(1.0,I4,shape(i,j,k,o),dlop,shape(k,i,m,n),0.0,I4bis,shape(j,n,o,m));

      LOp.clear();
      Contract(1.0,I4bis,shape(2,3),Environment::l[col-1][row],shape(0,1),0.0,LOp);

      // 2) construct left Sm operator
      I4bis.clear();
      Contract(1.0,I4,shape(i,j,k,o),dlom,shape(k,i,m,n),0.0,I4bis,shape(j,n,o,m));

      LOm.clear();
      Contract(1.0,I4bis,shape(2,3),Environment::l[col-1][row],shape(0,1),0.0,LOm);

      // 3) construct left Sz operator
      I4bis.clear();
      Contract(1.0,I4,shape(i,j,k,o),dloz,shape(k,i,m,n),0.0,I4bis,shape(j,n,o,m));

      LOz.clear();
      Contract(1.0,I4bis,shape(2,3),Environment::l[col-1][row],shape(0,1),0.0,LOz);

      // 4) finally construct new left unity
      Environment::construct_double_layer('V',peps(row,col),dlou);

      I4bis.clear();
      Contract(1.0,I4,shape(i,j,k,o),dlou,shape(k,i,m,n),0.0,I4bis,shape(j,n,o,m));

      LOu.clear();
      Contract(1.0,I4bis,shape(2,3),Environment::l[col-1][row],shape(0,1),0.0,LOu);

   }

   //last site on the right: close down on the incomings

   //1) first Lp with Sm
   Environment::construct_double_layer('V',peps(Ly-1,col),Sm,dlom);

   //paste right environment on
   tmp5.clear();
   Contract(1.0,Environment::r[col][Ly - 1],shape(1),dlom,shape(1),0.0,tmp5);

   //then left enviroment
   tmp6.clear();
   Contract(1.0,tmp5,shape(3),Environment::l[col-1][Ly-1],shape(1),0.0,tmp6);

   //move to a DArray<3> object
   RO[Ly - 3] = tmp6.reshape_clear(shape(Environment::r[col][Ly - 1].shape(0),dlom.shape(0),Environment::l[col-1][Ly - 1].shape(0)));

   //add to value
   val += 0.5 * Dot(LOp,RO[Ly - 3]);

   //2) then Lm with Sp
   Environment::construct_double_layer('V',peps(Ly-1,col),Sp,dlop);

   //paste right environment on
   tmp5.clear();
   Contract(1.0,Environment::r[col][Ly - 1],shape(1),dlop,shape(1),0.0,tmp5);

   //then bottom enviroment
   tmp6.clear();
   Contract(1.0,tmp5,shape(3),Environment::l[col-1][Ly-1],shape(1),0.0,tmp6);

   //move to a DArray<3> object
   RO[Ly - 3] = tmp6.reshape_clear(shape(Environment::r[col][Ly - 1].shape(0),dlop.shape(0),Environment::l[col-1][Ly - 1].shape(0)));

   //add to value
   val += 0.5 * Dot(LOm,RO[Ly - 3]);

   //3) then Lz with Sz
   Environment::construct_double_layer('V',peps(Ly-1,col),Sz,dloz);

   //paste top environment on
   tmp5.clear();
   Contract(1.0,Environment::r[col][Ly - 1],shape(1),dloz,shape(1),0.0,tmp5);

   //then bottom enviroment
   tmp6.clear();
   Contract(1.0,tmp5,shape(3),Environment::l[col-1][Ly-1],shape(1),0.0,tmp6);

   //move to a DArray<3> object
   RO[Ly - 3] = tmp6.reshape_clear(shape(Environment::r[col][Ly - 1].shape(0),dloz.shape(0),Environment::l[col-1][Ly - 1].shape(0)));

   //add to value
   val += Dot(LOz,RO[Ly - 3]);

}

// -- (3) -- || right column = Lx-1: again similar to overlap calculation

//first construct the right renormalized operators

//tmp comes out index (r,l)
tmp.clear();
Contract(1.0,Environment::r[Lx-2][Ly - 1],shape(1),Environment::l[Lx-2][Ly - 1],shape(1),0.0,tmp);

//reshape tmp to a 2-index array
R[Lx - 3] = tmp.reshape_clear(shape(Environment::r[Lx-2][Ly - 1].shape(0),Environment::l[Lx-2][Ly - 1].shape(0)));

//now construct the rest
for(int row = Ly - 2;row > 1;--row){

   I.clear();
   Contract(1.0,Environment::r[Lx-2][row],shape(2),R[row-1],shape(0),0.0,I);

   R[row-2].clear();
   Contract(1.0,I,shape(1,2),Environment::l[Lx-2][row],shape(1,2),0.0,R[row-2]);

}

//construct the left going operators on the first top site

//first S+
Environment::construct_double_layer('V',peps(0,Lx-1),Sp,dlsp);

//tmp comes out index (r,l)
Contract(1.0,dlsp,shape(1),Environment::l[Lx-2][0],shape(1),0.0,tmp);

Lp = tmp.reshape_clear(shape(dlsp.shape(2),Environment::l[Lx-2][0].shape(2)));

//then S-
Environment::construct_double_layer('V',peps(0,Lx-1),Sm,dlsm);

//tmp comes out index (r,l)
Contract(1.0,dlsm,shape(1),Environment::l[Lx-2][0],shape(1),0.0,tmp);

Lm = tmp.reshape_clear(shape(dlsm.shape(2),Environment::l[Lx-2][0].shape(2)));

//then Sz 
Environment::construct_double_layer('V',peps(0,Lx-1),Sz,dlsz);

//tmp comes out index (r,l)
Contract(1.0,dlsz,shape(1),Environment::l[Lx-2][0],shape(1),0.0,tmp);

Lz = tmp.reshape_clear(shape(dlsz.shape(2),Environment::l[Lx-2][0].shape(2)));

//and finally unity
Contract(1.0,Environment::r[Lx-2][0],shape(1),Environment::l[Lx-2][0],shape(1),0.0,tmp);

Lu = tmp.reshape_clear(shape(Environment::r[Lx-2][0].shape(2),Environment::l[Lx-2][0].shape(2)));

//middle of the chain:
for(int row = 1;row < Ly-1;++row){

   //first close down the +,- and z terms from the previous site

   //construct the right intermediate contraction (paste left to 'right')
   I.clear();
   Contract(1.0,Environment::l[Lx-2][row],shape(2),R[row - 1],shape(1),0.0,I);

   // 1) construct Sm double layer
   Environment::construct_double_layer('V',peps(row,Lx-1),Sm,dlsm);

   R[row-1].clear();
   Contract(1.0,dlsm,shape(1,2),I,shape(1,2),0.0,R[row - 1]);

   //contract with left S+
   val += 0.5 * Dot(Lp,R[row - 1]);

   // 2) construct Sp double layer
   Environment::construct_double_layer('V',peps(row,Lx-1),Sp,dlsp);

   R[row-1].clear();
   Contract(1.0,dlsp,shape(1,2),I,shape(1,2),0.0,R[row - 1]);

   //contract with left S-
   val += 0.5 * Dot(Lm,R[row - 1]);

   // 3) construct Sz double layer
   Environment::construct_double_layer('V',peps(row,Lx-1),Sz,dlsz);

   R[row-1].clear();
   Contract(1.0,dlsz,shape(1,2),I,shape(1,2),0.0,R[row - 1]);

   //contract with left Sz
   val += Dot(Lz,R[row - 1]);

   //construct left renormalized operators for next site: first paste bottom to Left unity
   I.clear();
   Contract(1.0,Lu,shape(1),Environment::l[Lx-2][row],shape(0),0.0,I);

   // 1) construct new Sm left operator
   Lm.clear();
   Contract(1.0,dlsm,shape(0,1),I,shape(0,1),0.0,Lm);

   // 2) construct new Sp left operator
   Lp.clear();
   Contract(1.0,dlsp,shape(0,1),I,shape(0,1),0.0,Lp);

   // 3) construct new Sz left operator
   Lz.clear();
   Contract(1.0,dlsz,shape(0,1),I,shape(0,1),0.0,Lz);

   // 4) finally construct new unity on the left
   Lu.clear();
   Contract(1.0,Environment::r[Lx-2][row],shape(0,1),I,shape(0,1),0.0,Lu);

}

//finally close down on last 'right' site

//1) Sm to close down Lp
Environment::construct_double_layer('V',peps(Ly-1,Lx-1),Sm,dlsm);

//tmp comes out index (r,l)
tmp.clear();
Contract(1.0,dlsm,shape(1),Environment::l[Lx-2][Ly - 1],shape(1),0.0,tmp);

//reshape tmp to a 2-index array
R[Ly - 3] = tmp.reshape_clear(shape(dlsm.shape(0),Environment::l[Lx-2][Ly - 1].shape(0)));

val += 0.5 * Dot(Lp,R[Ly-3]);

//2) Sp to close down Lm
Environment::construct_double_layer('V',peps(Ly-1,Lx-1),Sp,dlsp);

//tmp comes out index (r,l)
tmp.clear();
Contract(1.0,dlsp,shape(1),Environment::l[Lx-2][Ly - 1],shape(1),0.0,tmp);

//reshape tmp to a 2-index array
R[Ly - 3] = tmp.reshape_clear(shape(dlsp.shape(0),Environment::l[Lx-2][Ly - 1].shape(0)));

val += 0.5 * Dot(Lm,R[Ly-3]);

//3) Sz to close down Lz
Environment::construct_double_layer('V',peps(Ly-1,Lx-1),Sz,dlsz);

//tmp comes out index (r,l)
tmp.clear();
Contract(1.0,dlsz,shape(1),Environment::l[Lx-2][Ly - 1],shape(1),0.0,tmp);

//reshape tmp to a 2-index array
R[Ly - 3] = tmp.reshape_clear(shape(dlsz.shape(0),Environment::l[Lx-2][Ly - 1].shape(0)));

val += Dot(Lz,R[Ly-3]);
*/
}

return energy;

}
