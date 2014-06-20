#ifndef PROPAGATOR_H
#define PROPAGATOR_H

#include <iostream>
#include <fstream>
#include <vector>
#include <complex>

using std::ostream;
using std::vector;
using std::complex;

/**
 * @author Brecht Verstichel
 * @date 20-06-2014\n\n
 * This class Propagator is a class written for the construction and the action of the imaginary time propagator 
 * on a D=1 PEPS walker.
 */
class Propagator : public vector< TArray< complex<double> , 2 > > {

   public:
      
      //constructor
      Propagator();

      //copy constructor
      Propagator(const Propagator &);

      //destructor
      virtual ~Propagator();

      using vector::operator=;

      using vector::operator[];

      complex<double> gx() const;

      int gk() const;

      int gr() const;

      void set(complex<double>,int,int);

      void fill();

   private:

      //!random auxiliary field variable - shift
      complex<double> x;

      //!trotter index
      int k;

      //!type of operator (0=x,1=y or 2=z)
      int r;

};

#endif

/* vim: set ts=3 sw=3 expandtab :*/
