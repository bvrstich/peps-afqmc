#ifndef __BTAS_LAPACK_SYEV_IMPL_H
#define __BTAS_LAPACK_SYEV_IMPL_H 1

#include <lapack/types.h>

namespace btas {
namespace lapack {

template<typename T>
void syev (
   const int& order,
   const char& jobz,
   const char& uplo,
   const size_t& N,
         T* A,
   const size_t& ldA,
         T* W)
{
   BTAS_LAPACK_ASSERT(false, "syev must be specialized.");
}

inline void syev (
   const int& order,
   const char& jobz,
   const char& uplo,
   const size_t& N,
         float* A,
   const size_t& ldA,
         float* W)
{
   LAPACKE_ssyev(order, jobz, uplo, N, A, ldA, W);
}

inline void syev (
   const int& order,
   const char& jobz,
   const char& uplo,
   const size_t& N,
         double* A,
   const size_t& ldA,
         double* W)
{
   LAPACKE_dsyev(order, jobz, uplo, N, A, ldA, W);
}

} // namespace lapack
} // namespace btas

#endif // __BTAS_LAPACK_SYEV_IMPL_H
