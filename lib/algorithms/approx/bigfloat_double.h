    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/algorithms/approx/bigfloat_double.h

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

    See the full license in the file "LICENSE" in the top level distribution directory
    *************************************************************************************/
    /*  END LEGAL */
#include <math.h>

typedef double mfloat; 
class bigfloat {
private:

  mfloat x;

public:

  bigfloat() { }
  bigfloat(const bigfloat& y) { x=y.x; }
  bigfloat(const unsigned long u) { x=u; }
  bigfloat(const long i) { x=i; }
  bigfloat(const int i) { x=i;}
  bigfloat(const float d) { x=d;}
  bigfloat(const double d) {  x=d;}
  bigfloat(const char *str) { x=std::stod(std::string(str));}
  ~bigfloat(void) { }
  operator double (void) const { return (double)x; }
  static void setDefaultPrecision(unsigned long dprec) {
  }

  void setPrecision(unsigned long dprec) {
  }
  
  unsigned long getPrecision(void) const { return 64; }
  unsigned long getDefaultPrecision(void) const { return 64; }

  bigfloat& operator=(const bigfloat& y)     { x=y.x;    return *this;  }
  bigfloat& operator=(const unsigned long y) { x=y; return *this; }
  bigfloat& operator=(const signed long y)   { x=y; return *this; }
  bigfloat& operator=(const float y)    { x=y; return *this; }
  bigfloat& operator=(const double y)   { x=y; return *this; }

  size_t write(void);
  size_t read(void);

  /* Arithmetic Functions */

  bigfloat& operator+=(const bigfloat& y) { return *this = *this + y; }
  bigfloat& operator-=(const bigfloat& y) { return *this = *this - y; }
  bigfloat& operator*=(const bigfloat& y) { return *this = *this * y; }
  bigfloat& operator/=(const bigfloat& y) { return *this = *this / y; }

  friend bigfloat operator+(const bigfloat& x, const bigfloat& y) { 
    bigfloat a;
    a.x=x.x+y.x;
    return a;
  }

  friend bigfloat operator+(const bigfloat& x, const unsigned long y) {
    bigfloat a;
    a.x=x.x+y;
    return a;
  }

  friend bigfloat operator-(const bigfloat& x, const bigfloat& y) {
    bigfloat a;
    a.x=x.x-y.x;
    return a;
  }
  
  friend bigfloat operator-(const unsigned long x, const bigfloat& y) {
    bigfloat bx(x);
    return bx-y;
  }
  
  friend bigfloat operator-(const bigfloat& x, const unsigned long y) {
    bigfloat by(y);
    return x-by;
  }

  friend bigfloat operator-(const bigfloat& x) {
    bigfloat a;
    a.x=-x.x;
    return a;
  }

  friend bigfloat operator*(const bigfloat& x, const bigfloat& y) {
    bigfloat a;
    a.x=x.x*y.x;
    return a;
  }

  friend bigfloat operator*(const bigfloat& x, const unsigned long y) {
    bigfloat a;
    a.x=x.x*y;
    return a;
  }

  friend bigfloat operator/(const bigfloat& x, const bigfloat& y){
    bigfloat a;
    a.x=x.x/y.x;
    return a;
  }

  friend bigfloat operator/(const unsigned long x, const bigfloat& y){
    bigfloat bx(x);
    return bx/y;
  }

  friend bigfloat operator/(const bigfloat& x, const unsigned long y){
    bigfloat by(y);
    return x/by;
  }

  friend bigfloat sqrt_bf(const bigfloat& x){
    bigfloat a;
    a.x= sqrt(x.x);
    return a;
  }

  friend bigfloat sqrt_bf(const unsigned long x){
    bigfloat a(x);
    return sqrt_bf(a);
  }

  friend bigfloat abs_bf(const bigfloat& x){
    bigfloat a;
    a.x=fabs(x.x);
    return a;
  }

  friend bigfloat pow_bf(const bigfloat& a, long power) {
    bigfloat b;
    b.x=pow(a.x,power);
    return b;
  }

  friend bigfloat pow_bf(const bigfloat& a, bigfloat &power) {
    bigfloat b;
    b.x=pow(a.x,power.x);
    return b;
  }

  friend bigfloat exp_bf(const bigfloat& a) {
    bigfloat b;
    b.x=exp(a.x);
    return b;
  }

  /* Comparison Functions */
  friend int operator>(const bigfloat& x, const bigfloat& y) {
    return x.x>y.x;
  }

  friend int operator<(const bigfloat& x, const bigfloat& y) {
    return x.x<y.x;
  }

  friend int sgn(const bigfloat& x) {
    if ( x.x>=0 )  return 1;   
    else return 0;
  }

  /* Miscellaneous Functions */

  //  friend bigfloat& random(void);
};


