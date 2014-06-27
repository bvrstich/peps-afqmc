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
AFQMC::AFQMC(int Nw_in){

   this->Nw = Nw_in;

   P.resize(omp_num_threads);

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

#pragma omp parallel for
   for(int i = 0;i < walker.size();++i)
      walker[i].calc_properties('V',peps);

}

void AFQMC::walk(const int steps){

   complex<double> EP = gEP();

   char filename[200];
   sprintf(filename,"output/%dx%d/D=%d.txt",Lx,Ly,DT);

   ofstream output(filename,ios::trunc);

   output << "#Step\t\tE_P\t\tE_T\t" << endl;
   output.close();

#ifdef _DEBUG
   cout << "Energy at start = " << EP << endl;
   cout << "---------------------------------------------------------" << endl;
#endif

   for(int step = 0;step < steps;step++){

      //Propagate the walkers of each rank separately
      double wsum;

      complex<double> prev_EP = EP;
      
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

      write(step,0.5 * std::real(EP + prev_EP), ET);

      //Based on scaling, first control the population on each rank separately, and then balance the population over the ranks (uses MPI)
      PopulationControl(scaling);

      double max_ov = 0.0;
      double min_ov = 1.0;

      for(int i = 0;i < walker.size();++i){
         
         if(max_ov < std::abs(walker[i].gOverlap()))
            max_ov = std::abs(walker[i].gOverlap());

         if(min_ov > std::abs(walker[i].gOverlap()))
            min_ov = std::abs(walker[i].gOverlap());

      }

#ifdef _DEBUG
      cout << "Minimal Overlap:\t" << min_ov << endl;
      cout << "Maximal Overlap:\t" << max_ov << endl;
#endif

   }

}

/**
 * Here the trotter terms, propagator terms are applied to every walker individually.
 * @param option 'H'orizontal or 'V'ertical energy evaluation
 */
double AFQMC::propagate(char option){

   double sum = 0.0;

   double width = sqrt(2.0/Trotter::dtau);

   int num_rej = 0;

#pragma omp parallel for reduction(+: sum,num_rej)
   for(int i = 0;i < walker.size();i++){

#ifdef _OPENMP
      int myID = omp_get_thread_num();
#else
      int myID = 0;
#endif

      //backup the walker for stability
      backup_walker[myID].copy_essential(walker[i]);

      //now loop over the auxiliary fields:
      for(int k = 0;k < Trotter::n_trot;++k)
         for(int r = 0;r < 3;++r){

            double x = RN.normal();

            complex<double> shift = walker[i].gVL(k,r);

            //set the values
            P[myID].set(x + shift,k,r);

            //and fill the propagator
            P[myID].fill();

            //and apply it to the walker:
            walker[i].propagate(P[myID]);

         }

      //previous values of overlap and EL
      complex<double> prev_overlap = walker[i].gOverlap();
      complex<double> prev_EL = walker[i].gEL();

      //just for stability, doesn't do anything physical
      walker[i].normalize();

      //horizontal or vertical energy depending on iteration: (for speed reasons)
      walker[i].calc_properties(option,peps);

      complex<double> EL = walker[i].gEL();
      complex<double> overlap = walker[i].gOverlap();

      if( (std::real(EL) < std::real(prev_EL) - width) || (std::real(EL) > std::real(prev_EL) + width) ){//very rare event, will cause numerical unstability

         num_rej++;

         //copy the state back!
         walker[i].copy_essential(backup_walker[myID]);

      }
      else{

         //energy term for weight importance sampling
         double scale = exp( - Trotter::dtau * std::real(EL + prev_EL));

         //phase free projection
         scale *= std::max(0.0,cos(std::arg(overlap/prev_overlap)));

         walker[i].multWeight(scale);

      }

      sum += walker[i].gWeight();

   }

   return sum;

}

/**
 * redistribute the weights to stabilize the walk, keep the population in check
 */
void AFQMC::PopulationControl(double scaling){

   double minw = 1.0;
   double maxw = 1.0;

   double sum = 0.0;

   for(int i = 0;i < walker.size();i++){

      walker[i].multWeight(scaling);

      double weight = walker[i].gWeight();

      if(weight < minw)
         minw = weight;

      if(weight > maxw)
         maxw = weight;

      if (weight < 0.25){ //Energy doesn't change statistically

         int nCopies = (int) ( weight + rgen_pos<double>());

         if(nCopies == 0){

#ifdef _DEBUG
            cout << "Walker with weight " << weight << " will be deleted." << endl;
#endif

            walker.erase(walker.begin() + i);

         }
         else
            walker[i].sWeight(1.0);

      }

      if(weight > 1.5){ //statically energy doesn't change

         int nCopies =(int) ( weight + rgen_pos<double>());
         double new_weight = weight / (double) nCopies;

         walker[i].sWeight(new_weight);

#ifdef _DEBUG
         cout << "Walker with weight " << weight << " will be " << nCopies << " times copied with weight " << new_weight << "." << endl;
#endif

         for(int n = 1;n < nCopies;++n){

            Walker nw = walker[i];

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

   EP = 2.0 * projE_num / projE_den;

   return EP;

}

/**
 * write output to file
 */
void AFQMC::write(const int step,const double EP, const double ET){

   char filename[200];
   sprintf(filename,"output/%dx%d/D=%d.txt",Lx,Ly,DT);

   ofstream output(filename,ios::app);
   output.precision(16);
   output << step << "\t\t" << walker.size() << "\t" << EP << "\t\t" << ET << endl;
   output.close();

}
