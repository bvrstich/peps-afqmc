#ifndef RANDOM_H
#define RANDOM_H

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>

/*  Written by Sebastian Wouters <sebastianwouters@gmail.com> on September 4, 2013 */

class Random {

   public:
   
      //Constructor
      Random();
      
      //Destructor
      virtual ~Random();
      
      //Draw OpenMP and MPI thread safe random numbers from the uniform distribution (0,1)
      double operator()();

      //Draw OpenMP and MPI thread safe random numbers from the normal distribution (mean = 0, sigma = 1)
      double normal();
      
      //Tester of the Random number generator
      void test();
      
   private:
   
      //Number of OpenMP threads
      int num_omp_threads;
      
      //My MPI rank
      int rank;
   
      //The Mersenne Twister RN generator
      boost::random::mt19937 *mersenne;
      
      //The uniform real distribution
      boost::random::uniform_real_distribution<double> **dists;

      //The normal distribution
      boost::random::normal_distribution<double> **gauss;
      
};

#endif
