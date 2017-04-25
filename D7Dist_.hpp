#ifndef D7DIST_HPP
#define D7DIST_HPP

#include "D7.h"
#include "NCDist.h"

void MakeArray(const D7& g, double a[])
{
   for (size_t i = 0; i < 7; ++i) a[i] = g[i];
}

template<typename T1, typename T2>
double D7Dist_(const T1& t1, const T2& t2)
{
   const D7 d1 = t1;
   const D7 d2 = t2;
   double array1[7], array2[7];
   MakeArray(d1, array1);
   MakeArray(d2, array2);
   return D7CDist(array1, array2);
}


#endif  // D7DIST_HPP