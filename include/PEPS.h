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

/**
 * output stream operator overloaded for PEPS<T> 
 */
template<typename T>
ostream &operator<<(ostream &output,const PEPS<T> &peps_p){

   for(int r = 0;r < global::Ly;++r)
      for(int c = 0;c < global::Lx;++c){

         output << std::endl;
         output << "Tensor on site (" << r << "," << c << ")\t" << std::endl;
         output << std::endl;
         output << peps_p(r,c) << std::endl;

      }

   return output;

}

#endif

/* vim: set ts=3 sw=3 expandtab :*/
