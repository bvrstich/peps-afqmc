#ifndef GLOBAL_H
#define GLOBAL_H

#include <iostream>
#include <fstream>
#include <vector>

using std::ostream;
using std::vector;

using namespace btas;

/**
 * @author Brecht Verstichel
 * @date 29-04-2014\n\n
 * This class contains some global variables used throughout the program
 */
class Global {

   public:

      template<typename T>
         static T rgen();

      //!static Lattice object containing the info about the lattice
      static Lattice lat;

      static Random RN;

};

#endif

/* vim: set ts=3 sw=3 expandtab :*/
