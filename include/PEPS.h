#ifndef PEPS_H
#define PEPS_H

#include <iostream>
#include <fstream>
#include <vector>

using std::ostream;
using std::vector;

using namespace btas;

/**
 * @author Brecht Verstichel
 * @date 26-03-2014\n\n
 * This class PEPS is a class written for the construction of projected entangled pair states on a rectangular lattice
 */
template<typename T>
class PEPS : public vector< TArray<T,5> > {

   public:

      //empty
      PEPS();

      //constructor
      PEPS(int);

      //copy constructor
      PEPS(const PEPS &);

      //destructor
      virtual ~PEPS();

      int gD() const;

      void sD(int);

      const TArray<T,5> &operator()(int,int) const;

      TArray<T,5> &operator()(int,int);

      void load(const char *);

   private:

      //!cutoff virtual dimension
      int D;

};

#endif

/* vim: set ts=3 sw=3 expandtab :*/
