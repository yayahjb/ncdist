#ifndef D7_H
#define D7_H

#include "VecN.h"


class D6;
class DelaunayTetrahedron;
class G6;
class Cell;

class D7: public VecN {
public:
   friend class D6;

   D7( void );
   explicit D7( const D7& v );
   explicit D7( const G6& v );
   D7( const VecN& v );
   explicit D7(const D6& ds);
   explicit D7(const DelaunayTetrahedron& dt);
   D7( const std::string& s );
   D7( const std::vector<double>& v );
   D7(const Cell& c);

   D7& operator= (const D7& v);
   D7& operator= (const std::string& s);
   static double DistanceBetween( const D7& v1, const D7& v2 );
   friend std::ostream& operator<< ( std::ostream&, const D7& );
   double operator[](const size_t n) const { return m_vec[n]; }
   double& operator[](const size_t n) { return m_vec[n]; }
   //static size_t size(void) { return 7; }

};

#endif

