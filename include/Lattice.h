#ifndef LATTICE_H
#define LATTICE_H

#include <iostream>
#include <fstream>
#include <vector>

using std::ostream;
using std::vector;

/**
 * class written to define some parameters and store some lists of the lattice, as only one of these will be needed.
 */
class Lattice {

   public:
   
      Lattice();

      //Constructor
      Lattice(int,int,int);

      Lattice(const Lattice &);
      
      //Destructor
      virtual ~Lattice();

      void set(int,int,int);

      int gLx() const;

      int gLy() const;
      
      int gd() const;

      int grc2i(int,int) const;

      int gi2rc(int,int) const;
      
   private:

      //!x dimension of the lattice, nr of cols
      int Lx;

      //!y dimension of the lattice, nr of rows
      int Ly;

      //!physical dimension of sites
      int d;

      //!list relating row and col index to lattice index
      vector< vector<int> > rc2i;

      //!inverse list relatiing lattice index to row and column index
      vector< vector<int> > i2rc;
   
};

#endif
