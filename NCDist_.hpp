#ifndef NCDIST_HPP
#define NCDIST_HPP

#include "G6.h"
#include "NCDist.h"

void MakeArray(const G6& g, double a[6])
{
   for (size_t i = 0; i < 6; ++i) a[i] = g[i];
}

template<typename T1, typename T2>
double NCDist_(const T1& t1, const T2& t2)
{
   const G6 g1 = t1;
   const G6 g2 = t2;
   double array1[6], array2[6];
   MakeArray(g1, array1);
   MakeArray(g2, array2);
   return NCDist(array1, array2);
}


#endif  // NCDIST_HPP