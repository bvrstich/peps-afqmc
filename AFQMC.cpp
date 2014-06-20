#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <complex>
#include <vector>

#include "include.h"

using std::cout;
using std::endl;
using std::complex;
using std::vector;
using std::ofstream;
using std::ifstream;
using std::ios;

using namespace global;

/**
 * constructor of the AFQMC object, takes input parameters that define the QMC walk.
 * @param peps_in input trialwavefunction
 * @param Nwalkers number of Walker states
 */
AFQMC::AFQMC(const PEPS< complex<double> > &peps_in,int Nw_in){

   this->peps = peps_in;
   this->Nw = Nw_in;

   SetupWalkers();

}

/**
 * unnecessary destructor...
 */
AFQMC::~AFQMC(){ }

/**
 * initialize the walkers
 */
void AFQMC::SetupWalkers(){

   walker.resize(Nw);

   for(int i = 0;i < walker.size();++i)
      walker[i].calc_properties('V',peps);

}

void AFQMC::walk(const int steps){

   complex<double> EP = gEP();

   char filename[200];
   sprintf(filename,"output/%dx%d/D=%d.txt",Lx,Ly,peps.gD());

   ofstream output(filename,ios::trunc);

   output << "#Step\t\tE_P\t\tE_T\t" << endl;
   output.close();

#ifdef _DEBUG
   cout << "Energy at start = " << 2.0 * EP << endl;
   cout << "---------------------------------------------------------" << endl;
#endif

   for(int step = 1;step <= steps;step++){

      //Propagate the walkers of each rank separately
      double wsum;
      
      if(step % 2 == 0)
         wsum = propagate('H');
      else
         wsum = propagate('V');

      //Form the total sum of the walker weights and calculate the scaling for population control
      double avgw = wsum / (double)walker.size();

      double scaling = Nw / wsum;

      double ET = log(scaling)/Trotter::dtau;

      EP = gEP();

#ifdef _DEBUG
      cout << "        Step = " << step << endl;
      cout << "   # walkers = " << walker.size() << endl;
      cout << " avg(weight) = " << avgw << endl;
      cout << "         E_P = " << EP << endl;
      cout << "         E_T = " << ET << endl;
      cout << "---------------------------------------------------------" << endl;
#endif

      write(step,std::real(EP), ET);
/*
      //Based on scaling, first control the population on each rank separately, and then balance the population over the ranks (uses MPI)
      PopulationControl(scaling);

      double min_en = 0.0;
      double min_ov = 1.0;

      for(int i = 0;i < walker.size();++i){

         if(min_en > std::real(walker[i]->gEL()))
            min_en = std::real(walker[i]->gEL());

         if(min_ov > std::abs(walker[i]->gOverlap()))
            min_ov = std::abs(walker[i]->gOverlap());

      }

#ifdef _DEBUG
      cout << "Minimal Energy:\t" << min_en << endl;
      cout << "Minimal Overlap:\t" << min_ov << endl;
#endif
 */
   }

}

/**
 * Here the trotter terms, propagator terms are applied to every walker individually.
 * @param option 'H'orizontal or 'V'ertical energy evaluation
 */
double AFQMC::propagate(char option){

   double sum = 0.0;

#pragma omp parallel for reduction(+: sum,num_rej)
   for(int i = 0;i < walker.size();i++){

#ifdef _OPENMP
      int myID = omp_get_thread_num();
#else
      int myID = 0;
#endif

      //now loop over the auxiliary fields:
      for(int k = 0;k < Trotter::n_trot;++k)
         for(int r = 0;r < 3;++r){

            double x = RN.normal();

            complex<double> shift(0.0,0.0);// = walker[i].gVL(k,r);

            //set the values
            P.set(x + shift,k,r);

            //and fill the propagator
            P.fill();

            //and apply it to the walker:
            walker[i].propagate(P);

         }

      //previous values of overlap and EL
      complex<double> prev_overlap = walker[i].gOverlap();
      complex<double> prev_EL = walker[i].gEL();

      //just for stability, doesn't do anything physical
      walker[i].normalize();

      //horizontal or vertical energy depending on iteration: (for speed reasons)
      walker[i].calc_properties(option,peps);

      complex<double> overlap = walker[i].gOverlap();
      complex<double> EL = walker[i].gEL();

      //energy term for weight importance sampling
      double scale = exp( - Trotter::dtau * std::real(walker[i].gEL() + prev_EL));

      //phase free projection
      scale *= std::max(0.0,cos(std::arg(overlap/prev_overlap)));

      walker[i].multWeight(scale);

      sum += walker[i].gWeight();

   }

   return sum;

}

/**
 * redistribute the weights to stabilize the walk, keep the population in check
 */
/*
   void AFQMC::PopulationControl(double scaling){

   double minw = 1.0;
   double maxw = 1.0;

   double sum = 0.0;

   for(int i = 0;i < walker.size();i++){

   walker[i]->multWeight(scaling);

   double weight = walker[i]->gWeight();

   if(weight < minw)
   minw = weight;

   if(weight > maxw)
   maxw = weight;

   if (weight < 0.25){ //Energy doesn't change statistically

   int nCopies = (int) ( weight + Global::RN());

   if(nCopies == 0){

#ifdef _DEBUG
cout << "Walker with weight " << weight << " will be deleted." << endl;
#endif

delete walker[i];

walker.erase(walker.begin() + i);

}
else
walker[i]->sWeight(1.0);

}

if(weight > 1.5){ //statically energy doesn't change

int nCopies =(int) ( weight + Global::RN());
double new_weight = weight / (double) nCopies;

walker[i]->sWeight(new_weight);

#ifdef _DEBUG
cout << "Walker with weight " << weight << " will be " << nCopies << " times copied with weight " << new_weight << "." << endl;
#endif

for(int n = 1;n < nCopies;++n){

Walker *nw = new Walker(*walker[i]);

walker.push_back(nw);

}

}

sum += weight;

}

#ifdef _DEBUG
cout << endl;
cout << "total weight:\t" << sum << endl;
cout << endl;

cout << "The min. encountered weight is " << minw << " ." << endl;
cout << "The max. encountered weight is " << maxw << " ." << endl;
#endif

}
*/
/**
 * @return the total projected energy of the walkers at a certain timestep
 */
complex<double> AFQMC::gEP(){

   complex<double> projE_num = 0.0;
   complex<double> projE_den = 0.0;

   for(int wi = 0;wi < walker.size();wi++){

      complex<double> w_loc_en = walker[wi].gEL(); // <Psi_T | H | walk > / <Psi_T | walk >

      //For the projected energy
      projE_num   += walker[wi].gWeight() * w_loc_en;
      projE_den += walker[wi].gWeight();

   }

   complex<double> EP = 0.0;

   EP = projE_num / projE_den;

   return EP;

}

/**
 * write output to file
 */
void AFQMC::write(const int step,const double EP, const double ET){

   char filename[200];
   sprintf(filename,"output/%dx%d/D=%d.txt",Lx,Ly,peps.gD());

   ofstream output(filename,ios::app);
   output.precision(16);
   output << step << "\t\t" << walker.size() << "\t" << EP*2 << "\t\t" << ET << endl;
   output.close();

}
