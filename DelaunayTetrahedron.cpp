#include "Cell.h"
#include "D6.h"
#include "DelaunayTetrahedron.h"
#include "G6.h"
#include "D7.h"
#include "vector_3d.h"

std::ostream& operator<< (std::ostream& o, const DelaunayTetrahedron& dt) {
   for (size_t i = 0; i < 4; ++i)
      o << dt.m_vec[i];
   return o;
}

DelaunayTetrahedron::DelaunayTetrahedron(void) {
   m_vec.resize(4);
}

//                   DelaunayTetrahedron(const DelaunayTetrahedron& dy);
DelaunayTetrahedron::DelaunayTetrahedron(const DelaunayTetrahedron& dt) 
   : m_vec(dt.m_vec)
{
   m_vec.resize(4);
}

DelaunayTetrahedron::DelaunayTetrahedron(const G6& v6) {
   m_vec.resize(4);
   (*this) = DelaunayTetrahedron(Cell(v6));
}

DelaunayTetrahedron::DelaunayTetrahedron(const D7& v7) {
   m_vec.resize(4);
   (*this) = DelaunayTetrahedron( Cell(v7));
}

DelaunayTetrahedron::DelaunayTetrahedron(const std::string& v7) {
   m_vec.resize(4);
   throw;
}

const double rad = 180.0 / (4 * atan(1.0));
DelaunayTetrahedron::DelaunayTetrahedron(const Cell& c) {
   m_vec.resize(4);
   const Matrix_3x3 amat = Matrix_3x3::Cart(c[0], c[1], c[2], c[3] * rad, c[4] * rad, c[5] * rad);

   (*this).m_vec[0] = amat * Vector_3(1, 0, 0);
   (*this).m_vec[1] = amat * Vector_3(0, 1, 0);
   (*this).m_vec[2] = amat * Vector_3(0, 0, 1);
   (*this).m_vec[3] = -(m_vec[0]+m_vec[1]+m_vec[2]);
}

DelaunayTetrahedron::DelaunayTetrahedron(const D6& ds) {
   m_vec.resize(4);
   (*this) = DelaunayTetrahedron(Cell(ds));
}

Vector_3 DelaunayTetrahedron::operator[] (const size_t n) const {
   return m_vec[n];
}

Vector_3& DelaunayTetrahedron::operator[] (const size_t n) {
   return m_vec[n];
}

DelaunayTetrahedron& DelaunayTetrahedron::operator= (const std::string& s)
{
   throw;
   return *this;
}

DelaunayTetrahedron DelaunayTetrahedron::operator+ (const DelaunayTetrahedron& dt) const {
   const G6 v1(*this);
   const G6 v2(dt);
   return DelaunayTetrahedron(G6(v1 + v2));
}

DelaunayTetrahedron DelaunayTetrahedron::operator* (const double d) const
{
   DelaunayTetrahedron dt;
   for (size_t i = 0; i < dt.size(); ++i)
      dt[i] = ((*this)[i]) * d;
   return dt;
}

DelaunayTetrahedron operator* (const double d, const DelaunayTetrahedron& dt)
{
   return (dt*d);
}

DelaunayTetrahedron DelaunayTetrahedron::operator- (const DelaunayTetrahedron& dt) const {
   const G6 v1(*this);
   const G6 v2(dt);
   return DelaunayTetrahedron(G6(v1 - v2));
}

bool DelaunayTetrahedron::operator== (const DelaunayTetrahedron& dt)
{
   bool breturn = true;
   for (size_t i = 0; i < 4; ++i)
      breturn = breturn && (*this)[i] == dt[i];
   return breturn;
}

bool DelaunayTetrahedron::operator!= (const DelaunayTetrahedron& dt)
{
   return !((*this) == dt);
}
