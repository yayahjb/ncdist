#ifndef DeloneTETRAHEDRON_H
#define DeloneTETRAHEDRON_H

#include "Cell.h"
#include "vector_3d.h"

#include <vector>

   class D6;
   class G6;
   class D7;
   class Cell;

class DeloneTetrahedron
{
public:

   friend std::ostream& operator<< (std::ostream&, const DeloneTetrahedron&);
   friend DeloneTetrahedron operator* (const double d, const DeloneTetrahedron& dt);

   DeloneTetrahedron(void);
   DeloneTetrahedron(const DeloneTetrahedron& dt);
   explicit DeloneTetrahedron(const G6& v6);
   explicit DeloneTetrahedron(const D7& v7);
   explicit DeloneTetrahedron(const Cell& c);
   explicit DeloneTetrahedron(const D6& c);
   explicit DeloneTetrahedron(const std::string& c);
   ~DeloneTetrahedron(void) {}

   Vector_3  operator[] (const size_t n) const;
   Vector_3& operator[] (const size_t n);
   DeloneTetrahedron& operator= (const std::string& s);
   bool operator== (const DeloneTetrahedron& dt);
   bool operator!= (const DeloneTetrahedron& dt);
   DeloneTetrahedron operator* (const double d) const;
   DeloneTetrahedron operator+ (const DeloneTetrahedron& dt) const;
   DeloneTetrahedron operator- (const DeloneTetrahedron& dt) const;

   template<typename T>
   DeloneTetrahedron operator+ (const DeloneTetrahedron& dt) const {
      const T v1(*this);
      const T v2(Cell(dt));
      return DeloneTetrahedron(v1 + v2);
   }

   template<typename T>
   DeloneTetrahedron operator- (const DeloneTetrahedron& dt) const {
      const T v1(*this);
      const T v2(dt);
      return DeloneTetrahedron(v1 - v2);
   }


   static size_t size(void) { return 4; }

protected:
   std::vector<Vector_3> m_vec;
};


#endif