#ifndef NC_VEC_H
#define NC_VEC_H

#include "VecN.h"
#include "D7.h"

class D6;
class D7;
class DelaunayTetrahedron;
class Cell;

class G6 : public VecN {
public:
   friend class D6;

   G6(void);
   G6(const G6& v);
   G6(const D7& v);
   G6(const Cell& c);
   G6(const VecN& v);
   G6(const D6& ds);
   G6(const DelaunayTetrahedron& dt);
   G6(const std::string& s);
   G6(const std::vector<double>& v);

   G6& operator= (const G6& v);
   G6& operator= (const D7& v);
   G6& operator= (const D6& v);
   G6& operator= (const DelaunayTetrahedron& v);
   G6& operator= (const Cell& v);
   G6& operator= (const std::string& s);
   //G6& operator= ( const VecN& v );
   G6 InverseG6Vector(void) const;
   //static size_t size(void) { return 6; }
};

#endif



/*
#include <string>

#include "Mat66.h"

class Mat66;
class Matrix_3x3;


class G6 {
public:
   friend G6 operator* ( const double d, const G6&v );
   friend std::ostream& operator<< ( std::ostream&, const G6& );
   friend double dot( const G6& v1, const G6& v2 );

   G6( );
   G6( const std::string& s );
   G6( const Matrix_3x3& m );
   G6( const G6& v );

   ~G6( void );

   G6 operator* ( const G6& v2 ) const;
   G6 operator+ ( const G6& v2 ) const;
   G6 operator- ( const G6& v2 ) const;
   G6 operator*=( const double d );
   G6& operator=( const G6& v );
   G6 operator* ( const double d ) const;
   G6 operator/ ( const double d ) const;
   double operator[]( const size_t i ) const;
   double at( const size_t i ) const { return ( *this )[i]; }
   double& operator[]( const size_t i );
   bool operator== ( const G6& v2 ) const { return ( ( *this )-v2 ).norm( ) == 0.0; };
   bool operator!= ( const G6& v2 ) const { return !( ( *this )==v2 ); };

   size_t size( void ) const { return m_vec.size( ); }
   void clear( void ) { m_vec.clear( ); }
   double norm( void ) const;
   double dot( const G6& v2 ) const;
   Mat66 VV_Transpose( void ) const;
   static G6 abs( const G6& v );

   G6 InverseG6Vector( void ) const;

protected:
   std::vector<double> m_vec;

};

G6 operator* ( const double d, const G6&v );

#endif
*/
