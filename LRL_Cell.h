#ifndef LRL_CELL_H
#define LRL_CELL_H

#ifdef _MSC_VER
#define USE_LOCAL_HEADERS
#endif
#ifndef USE_LOCAL_HEADERS
//#include <rhrand.h>
#else
#include "rhrand.h"
#endif

#include "BasisBase.h"
#include "C3.h"
#include "S6.h"
#include "D7.h"
#include "B4.h"
#include "G6.h"
#include "MatG6.h"
#include "Mat66.h"
#include "VecN.h"

#include <string>
#include <vector>

class C3;
class S6;
class D7;
class B4;
class G6;
class MatG6;
class Mat66;

class LRL_Cell : private BasisBase<LRL_Cell>
{
public:

   static double randomLatticeNormalizationConstant;
   static double randomLatticeNormalizationConstantSquared;

   friend std::ostream& operator<< (std::ostream& o, const LRL_Cell& c);
   friend LRL_Cell operator* (const double d, const LRL_Cell& v);
   LRL_Cell(void);
   LRL_Cell(const LRL_Cell& c);
   LRL_Cell(const std::string& s); // cell with angles in degrees from text
   LRL_Cell(const S6& ds);
   LRL_Cell(const C3& c3);
   LRL_Cell(const B4& dt);
   LRL_Cell(const double a, const double b, const double c,
      const double alpha, const double beta, const double gamma);
   LRL_Cell(const G6& v6); // specialized for G6
   LRL_Cell(const D7& v7); // specialized for d7
   ~LRL_Cell(void);

   double Volume(void) const;

   bool operator== (const LRL_Cell& cl) const;
   bool operator!= (const LRL_Cell& cl) const;
   double operator[](const size_t n) const;
   double& operator[](const size_t n);
   double DistanceBetween(const LRL_Cell& v1, const LRL_Cell& v2);
   LRL_Cell Inverse( void ) const;
   G6 Cell2V6( void ) const;
   G6 Cell2V6( const LRL_Cell& c );
   Mat66 LatSymMat66( const std::string& latsym );
   Mat66 LatSymMat66( const std::string& latsym, const LRL_Cell& c );
   MatG6 LatSymMatG6( const std::string& latsym );
   MatG6 LatSymMatG6( const std::string& latsym, const LRL_Cell& c );
   LRL_Cell& operator= (const std::string& s);
   LRL_Cell& operator= (const D7& v);
   LRL_Cell& operator= (const C3& c3);
   LRL_Cell& operator= (const S6& v);
   LRL_Cell& operator= (const B4& v);
   LRL_Cell& operator= (const LRL_Cell& v);
   LRL_Cell& operator= (const G6& v);
   LRL_Cell& operator/= (const double d);
   LRL_Cell& operator*= (const double d);
   LRL_Cell operator+ (const LRL_Cell& c) const;
   LRL_Cell operator- (const LRL_Cell& c) const;
   LRL_Cell operator* (const double d) const;
   LRL_Cell operator/ (const double d) const;
   LRL_Cell operator- (void) const { return *this; } // unary
   LRL_Cell& operator+= (const LRL_Cell& cl);
   LRL_Cell& operator-= (const LRL_Cell& cl);

   size_t size(void) const { return 6; }
   double norm() const { return (*this).Volume();}
   double norm(const LRL_Cell& cell) const { return cell.Volume();}
   double Norm() const { return (*this).Volume();}
   double Norm(const LRL_Cell& cell) const { return cell.Volume();}

   std::vector<double> GetVector(void) const { return m_cell; }
   double* data() const { return const_cast<double*>(m_cell.data()); }

   void SetVector(const std::vector<double>& v) { m_cell = v; }
   bool GetValid(void) const { return m_valid; }
   bool IsValid(void) const { return m_valid; }
   void SetValid(const bool b) { m_valid = b; };
   double at(const size_t n) const { return (*this)[n]; }


   static LRL_Cell GetPrimitiveCell(const std::string& latsym, const LRL_Cell& c);
   static G6 GetPrimitiveV6Vector(const std::string& latsym, const LRL_Cell& c);

   bool IsRhomobhedralAsHex(void) const;
   bool IsRhomobhedralAsHex( const G6& ) const;

   static LRL_Cell rand();
   static LRL_Cell randDeloneReduced();
   static LRL_Cell randDeloneUnreduced();
   static LRL_Cell rand(const double d);
   static LRL_Cell randDeloneReduced(const double d);
   static LRL_Cell randDeloneUnreduced(const double d);

   static std::string GetName(void) { return "LRL_Cell, unit cell"; }

   static double GetNormalization(void) { return randomLatticeNormalizationConstant; }

protected:
   std::vector<double> m_cell;
   bool m_valid;
};

#endif // LRL_CELL_H
