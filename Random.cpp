#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <omp.h>
#include <time.h>

#include "include.h"

/*  Written by Sebastian Wouters <sebastianwouters@gmail.com> on September 4, 2013 */

using namespace std;

Random::Random(){

#ifdef _OPENMP
   num_omp_threads = omp_get_max_threads();
#else
   num_omp_threads = 1;
#endif

   mersenne = new boost::random::mt19937 [num_omp_threads];

   for (int cnt=0; cnt<num_omp_threads; cnt++)
      mersenne[cnt].seed( time(0) + cnt*cnt*23 );

   dists = new boost::random::uniform_real_distribution<double> * [num_omp_threads];

   for (int cnt=0; cnt<num_omp_threads; cnt++)
      dists[cnt] = new boost::random::uniform_real_distribution<double> (0,1); 

   gauss = new boost::random::normal_distribution<double> * [num_omp_threads];

   for (int cnt=0; cnt<num_omp_threads; cnt++)
      gauss[cnt] = new boost::random::normal_distribution<double> (0.0,1.0); 

}

Random::~Random(){

   for (int cnt=0; cnt<num_omp_threads; cnt++)
      delete dists[cnt];

   delete [] dists;

   for (int cnt=0; cnt<num_omp_threads; cnt++)
      delete gauss[cnt];

   delete [] gauss;

   delete [] mersenne;

}

double Random::operator()(){

#ifdef _OPENMP
   int tid = omp_get_thread_num();
#else
   int tid = 0;
#endif

   return (*dists[tid])(mersenne[tid]);

}

double Random::normal(){

#ifdef _OPENMP
   int tid = omp_get_thread_num();
#else
   int tid = 0;
#endif

   return (*gauss[tid])(mersenne[tid]);

}
