#include "Cell.h"
#include "D6.h"
#include "G6.h"
#include "D7.h"
#include "Vec_N_Tools.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

D7::D7( void )
: VecN( 7 ) {
   m_dim = 7;
   m_valid = true;
}

D7::D7( const D7& v )
: VecN( 7 ) {
   m_dim = 7;
   m_vec = v.m_vec;
   m_valid = v.m_valid;
}

double D7::DistanceBetween( const D7& v1, const D7& v2 ) {
   double sum = 0.0;
   const size_t size1 = sizeof(v1);
   const size_t size2 = sizeof(v1[0]);
   for( size_t i=0; i<7; ++i )  sum += (v1[i]-v2[i])*(v1[i]-v2[i]);
   return sqrt( sum );
}
D7::D7( const G6& v6 )
: VecN( 7 ) {
   const double& g1 = v6[0];
   const double& g2 = v6[1];
   const double& g3 = v6[2];
   const double& g4 = v6[3];
   const double& g5 = v6[4];
   const double& g6 = v6[5];
   double& d1 = m_vec[0];
   double& d2 = m_vec[1];
   double& d3 = m_vec[2];
   double& d4 = m_vec[3];
   double& d5 = m_vec[4];
   double& d6 = m_vec[5];
   double& d7 = m_vec[6];
   d1 = g1;
   d2 = g2;
   d3 = g3;
   d4 = g1 + g2 + g3 + g6 + g5 + g4; // (a+b+c).(a+b+c)
   d5 = g2 + g3 + g4;  // (b+c).(b+c)
   d6 = g1 + g3 + g5;  // (a+c).(a+c)
   d7 = g1 + g2 + g6; // (a+b).(a+b)
   m_valid = true;
}

D7::D7(const D6& ds) {
   m_dim = 7;
   (*this) = G6(ds);
   m_valid = true;
}

D7::D7(const DelaunayTetrahedron& dt) {
   m_dim = 7;
   (*this) = G6(dt);
   m_valid = true;
}

D7::D7(const Cell& c) {
   (*this) = G6(c.Cell2V6());
}

D7::D7( const VecN& v )
: VecN( v ) {
   if (v.size( ) == 6) {
      G6 g6(v);
      const D7 v7(g6);
      m_vec.resize(7);
      for( size_t i=0; i<7; ++i ) m_vec[i] = v7[i];
   }
  m_valid = true;
   //   *this = v;
}

D7::D7( const std::string& s ) 
: VecN(7)
{
  m_vec = Vec_N_Tools::FromString( s );
  m_dim = m_vec.size();
  m_valid = true;
  if ( m_dim != 7 ) throw "bad dimension in D7 from a string";
}

D7::D7( const std::vector<double>& v ) 
: VecN( v.size() )
{
  m_dim = v.size();
  if ( m_dim != 7 ) throw "bad dimension in D7 from a string";
  m_vec = v;
  m_valid = true;
}

D7& D7::operator= ( const D7& v ) {
   const double d = v[0];
   m_vec = v.m_vec;
   m_dim = v.m_dim;
   m_valid = v.m_valid;
   return *this;
}

D7& D7::operator= (const std::string& s)
{
   (*this) = D7(s);
   return *this;
}

std::ostream& operator<< ( std::ostream& o, const D7& v ) {
   std::streamsize oldPrecision = o.precision();
   o << std::fixed <<std::setprecision(3);
   for( size_t i=0; i<v.size( ); ++i )
      o << std::setw(9) << v[i] << " ";
   o << std::setprecision(oldPrecision);
   o.unsetf(std::ios::floatfield);
   return o;
}
