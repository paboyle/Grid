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


