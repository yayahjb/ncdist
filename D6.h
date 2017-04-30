#ifndef D6_H
#define D6_H

#include <vector>

class Cell;
class G6;
class D7;
class DeloneTetrahedron;

#include <ostream>
#include <string>

class D6 : public std::vector<double>
{
public:
   friend std::ostream& operator<< (std::ostream&, const D6&);
   friend D6 operator* (const double d, const D6& ds);

   D6(void);
   D6(const G6& v6);
    D6(const double v[6]);
   D6(const D7& v7 );
   D6(const Cell& c);
   D6(const std::string& s);
   D6(const DeloneTetrahedron& del);
      ~D6(void) {}

   D6 operator- (const D6& ds) const;
   D6 operator+ (const D6& ds) const;
   D6 operator* (const double& d) const;
   D6& operator=(const std::string& s);
   double norm(const D6& ds) const;
   double norm(void) const;
   double Norm(const D6& ds) const { return norm(ds); }
   double Norm(void) const { return (*this).norm(); }
   static size_t size(void) { return 6; }

};


#endif
