/*
  C.Kelly Jan 2020 based on implementation by M. Clark May 2005

  AlgRemezGeneral is an implementation of the Remez algorithm for approximating an arbitrary function by a rational polynomial 
  It includes optional restriction to odd/even polynomials for the numerator and/or denominator
*/

#ifndef INCLUDED_ALG_REMEZ_GENERAL_H
#define INCLUDED_ALG_REMEZ_GENERAL_H

#include <stddef.h>
#include <Grid/GridStd.h>

#ifdef HAVE_LIBGMP
#include "bigfloat.h"
#else
#include "bigfloat_double.h"
#endif


class AlgRemezGeneral{
 public:
  enum PolyType { Even, Odd, Full };

 private:

  // In GSL-style, pass the function as a function pointer. Any data required to evaluate the function is passed in as a void pointer
  bigfloat (*f)(bigfloat x, void *data);
  void *data;

  // The approximation parameters
  std::vector<bigfloat> param;
  bigfloat norm;

  // The number of non-zero terms in the numerator and denominator
  int n, d;
  // The numerator and denominator degree (i.e.  the largest power)
  int pow_n, pow_d;
  
  // Specify if the numerator and/or denominator are odd/even polynomials
  PolyType num_type;
  PolyType den_type;
  std::vector<int> num_pows; //contains the mapping, with -1 if not present
  std::vector<int> den_pows;

  // The bounds of the approximation
  bigfloat apstrt, apwidt, apend;

  // Variables used to calculate the approximation
  int nd1, iter;
  std::vector<bigfloat> xx;
  std::vector<bigfloat> mm;
  std::vector<bigfloat> step;

  bigfloat delta, spread;
  
  // Variables used in search
  std::vector<bigfloat> yy;

  // Variables used in solving linear equations
  std::vector<bigfloat> A;
  std::vector<bigfloat> B;
  std::vector<int> IPS;

  // The number of equations we must solve at each iteration (n+d+1)
  int neq;

  // The precision of the GNU MP library
  long prec;

  // Initialize member variables associated with the polynomial's properties
  void setupPolyProperties(int num_degree, int den_degree, PolyType num_type_in, PolyType den_type_in);

  // Initial values of maximal and minmal errors
  void initialGuess();

  // Initialise step sizes
  void stpini();

  // Initialize the algorithm
  void reinitializeAlgorithm();

  // Solve the equations
  void equations();

  // Search for error maxima and minima
  void search(); 

  // Calculate function required for the approximation
  inline bigfloat func(bigfloat x) const{
    return f(x, data);
  }

  // Compute size and sign of the approximation error at x
  bigfloat getErr(bigfloat x, int *sign) const;

  // Solve the system AX=B   where X = param
  int simq();

  // Evaluate the rational form P(x)/Q(x) using coefficients from the solution vector param
  bigfloat approx(bigfloat x) const;

 public:
  
  AlgRemezGeneral(double lower, double upper, long prec,
		  bigfloat (*f)(bigfloat x, void *data), void *data);

  inline int getDegree(void) const{ 
    assert(n==d);
    return n;
  }
  // Reset the bounds of the approximation
  inline void setBounds(double lower, double upper) {
    apstrt = lower;
    apend = upper;
    apwidt = apend - apstrt;
  }

  // Get the bounds of the approximation
  inline void getBounds(double &lower, double &upper) const{ 
    lower=(double)apstrt;
    upper=(double)apend;
  }

  // Run the algorithm to generate the rational approximation
  double generateApprox(int num_degree, int den_degree, 
			PolyType num_type, PolyType den_type,
			const double tolerance = 1e-15, const int report_freq = 1000);
  
  inline double generateApprox(int num_degree, int den_degree, 
			       const double tolerance = 1e-15, const int report_freq = 1000){
    return generateApprox(num_degree, den_degree, Full, Full, tolerance, report_freq);
  }
  
  // Evaluate the rational form P(x)/Q(x) using coefficients from the
  // solution vector param
  inline double evaluateApprox(double x) const{
    return (double)approx((bigfloat)x);
  }

  // Evaluate the rational form Q(x)/P(x) using coefficients from the solution vector param
  inline double evaluateInverseApprox(double x) const{
    return 1.0/(double)approx((bigfloat)x);
  }  

  // Calculate function required for the approximation
  inline double evaluateFunc(double x) const{
    return (double)func((bigfloat)x);
  }

  // Calculate inverse function required for the approximation
  inline double evaluateInverseFunc(double x) const{
    return 1.0/(double)func((bigfloat)x);
  }

  // Dump csv of function, approx and error
  void csv(std::ostream &os = std::cout) const;

  // Get the coefficient of the term x^i in the numerator
  inline double getCoeffNum(const int i) const{    
    return num_pows[i] == -1 ? 0. : double(param[num_pows[i]]);
  }
  // Get the coefficient of the term x^i in the denominator
  inline double getCoeffDen(const int i) const{ 
    if(i == pow_d) return 1.0;
    else return den_pows[i] == -1 ? 0. : double(param[den_pows[i]+n+1]); 
  }
};

#endif
