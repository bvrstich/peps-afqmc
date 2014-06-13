#ifndef COMPRESS_H
#define COMPRESS_H

#include <iostream>
#include <iomanip>

using namespace btas;

namespace compress {

   void init_ro(const BTAS_SIDE &,std::vector< TArray<complex<double>,2> > &,const MPS &,const MPS &);

   void update_L(int,std::vector< TArray<complex<double>,2> > &,const MPS &,const MPS &);

   void update_R(int,std::vector< TArray<complex<double>,2> > &,const MPS &,const MPS &);

}

#endif
