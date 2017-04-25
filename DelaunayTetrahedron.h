#ifndef DELAUNAYTETRAHEDRON_H
#define DELAUNAYTETRAHEDRON_H

#include "Cell.h"
#include "vector_3d.h"

#include <vector>

   class D6;
   class G6;
   class D7;
   class Cell;

class DelaunayTetrahedron
{
public:

   friend std::ostream& operator<< (std::ostream&, const DelaunayTetrahedron&);
   friend DelaunayTetrahedron operator* (const double d, const DelaunayTetrahedron& dt);

   DelaunayTetrahedron(void);
   DelaunayTetrahedron(const DelaunayTetrahedron& dt);
   explicit DelaunayTetrahedron(const G6& v6);
   explicit DelaunayTetrahedron(const D7& v7);
   explicit DelaunayTetrahedron(const Cell& c);
   explicit DelaunayTetrahedron(const D6& c);
   explicit DelaunayTetrahedron(const std::string& c);
   ~DelaunayTetrahedron(void) {}

   Vector_3  operator[] (const size_t n) const;
   Vector_3& operator[] (const size_t n);
   DelaunayTetrahedron& operator= (const std::string& s);
   bool operator== (const DelaunayTetrahedron& dt);
   bool operator!= (const DelaunayTetrahedron& dt);
   DelaunayTetrahedron operator* (const double d) const;
   DelaunayTetrahedron operator+ (const DelaunayTetrahedron& dt) const;
   DelaunayTetrahedron operator- (const DelaunayTetrahedron& dt) const;

   template<typename T>
   DelaunayTetrahedron operator+ (const DelaunayTetrahedron& dt) const {
      const T v1(*this);
      const T v2(Cell(dt));
      return DelaunayTetrahedron(v1 + v2);
   }

   template<typename T>
   DelaunayTetrahedron operator- (const DelaunayTetrahedron& dt) const {
      const T v1(*this);
      const T v2(dt);
      return DelaunayTetrahedron(v1 - v2);
   }


   static size_t size(void) { return 4; }

protected:
   std::vector<Vector_3> m_vec;
};


#endif