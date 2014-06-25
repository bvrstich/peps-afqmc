#ifndef AFQMC_H
#define AFQMC_H

#include <iostream>
#include <iomanip>
#include <complex>
#include <vector>

using std::complex;
using std::vector;

#include "Propagator.h"

class AFQMC {

   public:
   
      //constructor with input trialwavefunction
      AFQMC(int);
      
      //Destructor
      virtual ~AFQMC();
      
      //Let the walkers propagate for steps steps
      void walk(int);

      //Propagate my population of walkers for 1 timestep. Return the sum of the coeff of my walkers.
      double propagate(char);
      
      //Control the population of walkers based on scaling * weight
      void PopulationControl(double);

      //Calculate the single walker projected energies, update the energy history, calculate the fluctuation metric, and the total projected energy
      complex<double> gEP();

      //Write the projected energy, target energy
      void write(int,double,double);

      //Setup the walkers
      void SetupWalkers();

   private:
      
      //!The total desired number of walkers
      int Nw;
      
      //propagator
      std::vector< Propagator > P;
      
      //!The walkers
      std::vector<Walker> walker;
      
};

#endif
