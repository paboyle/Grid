#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string>
#include<iostream>
#include<iomanip>
#include<cassert>

#include<Grid/algorithms/approx/RemezGeneral.h>


// Constructor
AlgRemezGeneral::AlgRemezGeneral(double lower, double upper, long precision,
				 bigfloat (*f)(bigfloat x, void *data), void *data): f(f), 
										     data(data), 
										     prec(precision),
										     apstrt(lower), apend(upper), apwidt(upper - lower),
										     n(0), d(0), pow_n(0), pow_d(0)
{
  bigfloat::setDefaultPrecision(prec);

  std::cout<<"Approximation bounds are ["<<apstrt<<","<<apend<<"]\n";
  std::cout<<"Precision of arithmetic is "<<precision<<std::endl;
}

//Determine the properties of the numerator and denominator polynomials
void AlgRemezGeneral::setupPolyProperties(int num_degree, int den_degree, PolyType num_type_in, PolyType den_type_in){
  pow_n = num_degree;
  pow_d = den_degree;

  if(pow_n % 2 == 0 && num_type_in == PolyType::Odd) assert(0);
  if(pow_n % 2 == 1 && num_type_in == PolyType::Even) assert(0);

  if(pow_d % 2 == 0 && den_type_in == PolyType::Odd) assert(0);
  if(pow_d % 2 == 1 && den_type_in == PolyType::Even) assert(0);

  num_type = num_type_in;
  den_type = den_type_in;

  num_pows.resize(pow_n+1);
  den_pows.resize(pow_d+1);

  int n_in = 0;
  bool odd = num_type == PolyType::Full || num_type == PolyType::Odd;
  bool even = num_type == PolyType::Full || num_type == PolyType::Even;
  for(int i=0;i<=pow_n;i++){
    num_pows[i] = -1;
    if(i % 2 == 0 && even) num_pows[i] = n_in++;
    if(i % 2 == 1 && odd) num_pows[i] = n_in++;
  }

  std::cout << n_in << " terms in numerator" << std::endl;
  --n_in; //power is 1 less than the number of terms, eg  pow=1   a x^1  + b x^0

  int d_in = 0;
  odd = den_type == PolyType::Full || den_type == PolyType::Odd;
  even = den_type == PolyType::Full || den_type == PolyType::Even;
  for(int i=0;i<=pow_d;i++){
    den_pows[i] = -1;
    if(i % 2 == 0 && even) den_pows[i] = d_in++;
    if(i % 2 == 1 && odd) den_pows[i] = d_in++;
  }

  std::cout << d_in << " terms in denominator" << std::endl;
  --d_in;

  n = n_in;
  d = d_in;
}

//Setup algorithm
void AlgRemezGeneral::reinitializeAlgorithm(){
  spread = 1.0e37;
  iter = 0;

  neq = n + d + 1; //not +2 because highest-power term in denominator is fixed to 1

  param.resize(neq);
  yy.resize(neq+1);

  //Initialize linear equation temporaries
  A.resize(neq*neq);
  B.resize(neq);
  IPS.resize(neq);

  //Initialize maximum and minimum errors
  xx.resize(neq+2);
  mm.resize(neq+1);
  initialGuess();

  //Initialize search steps
  step.resize(neq+1);
  stpini();
}

double AlgRemezGeneral::generateApprox(const int num_degree, const int den_degree, 
				       const PolyType num_type_in, const PolyType den_type_in, 
				       const double _tolerance, const int report_freq){
  //Setup the properties of the polynomial
  setupPolyProperties(num_degree, den_degree, num_type_in, den_type_in);

  //Setup the algorithm
  reinitializeAlgorithm();

  bigfloat tolerance = _tolerance;

  //Iterate until convergance
  while (spread > tolerance) { 
    if (iter++ % report_freq==0)
      std::cout<<"Iteration " <<iter-1<<" spread "<<(double)spread<<" delta "<<(double)delta << std::endl; 

    equations();
    if (delta < tolerance) {
      std::cout<<"Iteration " << iter-1 << " delta too small (" << delta << "<" << tolerance << "), try increasing precision\n";
      assert(0);
    };    
    assert( delta>= tolerance );

    search();
  }

  int sign;
  double error = (double)getErr(mm[0],&sign);
  std::cout<<"Converged at "<<iter<<" iterations; error = "<<error<<std::endl;

  // Return the maximum error in the approximation
  return error;
}


// Initial values of maximal and minimal errors
void AlgRemezGeneral::initialGuess(){
  // Supply initial guesses for solution points
  long ncheb = neq;			// Degree of Chebyshev error estimate

  // Find ncheb+1 extrema of Chebyshev polynomial
  bigfloat a = ncheb;
  bigfloat r;

  mm[0] = apstrt;
  for (long i = 1; i < ncheb; i++) {
    r = 0.5 * (1 - cos((M_PI * i)/(double) a));
    //r *= sqrt_bf(r);
    r = (exp((double)r)-1.0)/(exp(1.0)-1.0);
    mm[i] = apstrt + r * apwidt;
  }
  mm[ncheb] = apend;

  a = 2.0 * ncheb;
  for (long i = 0; i <= ncheb; i++) {
    r = 0.5 * (1 - cos(M_PI * (2*i+1)/(double) a));
    //r *= sqrt_bf(r); // Squeeze to low end of interval
    r = (exp((double)r)-1.0)/(exp(1.0)-1.0);
    xx[i] = apstrt + r * apwidt;
  }
}

// Initialise step sizes
void AlgRemezGeneral::stpini(){
  xx[neq+1] = apend;
  delta = 0.25;
  step[0] = xx[0] - apstrt;
  for (int i = 1; i < neq; i++) step[i] = xx[i] - xx[i-1];
  step[neq] = step[neq-1];
}

// Search for error maxima and minima
void AlgRemezGeneral::search(){
  bigfloat a, q, xm, ym, xn, yn, xx1;
  int emsign, ensign, steps;

  int meq = neq + 1;

  bigfloat eclose = 1.0e30;
  bigfloat farther = 0l;

  bigfloat xx0 = apstrt;

  for (int i = 0; i < meq; i++) {
    steps = 0;
    xx1 = xx[i]; // Next zero
    if (i == meq-1) xx1 = apend;
    xm = mm[i];
    ym = getErr(xm,&emsign);
    q = step[i];
    xn = xm + q;
    if (xn < xx0 || xn >= xx1) {	// Cannot skip over adjacent boundaries
      q = -q;
      xn = xm;
      yn = ym;
      ensign = emsign;
    } else {
      yn = getErr(xn,&ensign);
      if (yn < ym) {
	q = -q;
	xn = xm;
	yn = ym;
	ensign = emsign;
      }
    }
  
    while(yn >= ym) {		// March until error becomes smaller.
      if (++steps > 10)
      	break;
      
      ym = yn;
      xm = xn;
      emsign = ensign;
      a = xm + q;
      if (a == xm || a <= xx0 || a >= xx1)
	break;// Must not skip over the zeros either side.      

      xn = a;
      yn = getErr(xn,&ensign);
    }

    mm[i] = xm;			// Position of maximum
    yy[i] = ym;			// Value of maximum

    if (eclose > ym) eclose = ym;
    if (farther < ym) farther = ym;

    xx0 = xx1; // Walk to next zero.
  } // end of search loop

  q = (farther - eclose);	// Decrease step size if error spread increased

  if (eclose != 0.0) q /= eclose; // Relative error spread

  if (q >= spread)
    delta *= 0.5; // Spread is increasing; decrease step size
  
  spread = q;

  for (int i = 0; i < neq; i++) {
    q = yy[i+1];
    if (q != 0.0) q = yy[i] / q  - (bigfloat)1l;
    else q = 0.0625;
    if (q > (bigfloat)0.25) q = 0.25;
    q *= mm[i+1] - mm[i];
    step[i] = q * delta;
  }
  step[neq] = step[neq-1];
  
  for (int i = 0; i < neq; i++) {	// Insert new locations for the zeros.
    xm = xx[i] - step[i];

    if (xm <= apstrt)
      continue;

    if (xm >= apend)
      continue;

    if (xm <= mm[i])
      xm = (bigfloat)0.5 * (mm[i] + xx[i]);    

    if (xm >= mm[i+1])
      xm = (bigfloat)0.5 * (mm[i+1] + xx[i]);
    
    xx[i] = xm;
  }
}

// Solve the equations
void AlgRemezGeneral::equations(){
  bigfloat x, y, z;
  bigfloat *aa;
  
  for (int i = 0; i < neq; i++) {	// set up the equations for solution by simq()
    int ip = neq * i;		// offset to 1st element of this row of matrix
    x = xx[i];			// the guess for this row
    y = func(x);		// right-hand-side vector

    z = (bigfloat)1l;
    aa = A.data()+ip;
    int t = 0;
    for (int j = 0; j <= pow_n; j++) {
      if(num_pows[j] != -1){ *aa++ = z; t++; }
      z *= x;
    }
    assert(t == n+1);

    z = (bigfloat)1l;
    t = 0;
    for (int j = 0; j < pow_d; j++) {
      if(den_pows[j] != -1){ *aa++ = -y * z; t++; }
      z *= x;
    }
    assert(t == d);

    B[i] = y * z;		// Right hand side vector
  }

  // Solve the simultaneous linear equations.
  if (simq()){
    std::cout<<"simq failed\n";
    exit(0);
  }
}


// Evaluate the rational form P(x)/Q(x) using coefficients
// from the solution vector param
bigfloat AlgRemezGeneral::approx(const bigfloat x) const{
  // Work backwards toward the constant term.
  int c = n;
  bigfloat yn = param[c--];		// Highest order numerator coefficient
  for (int i = pow_n-1; i >= 0; i--) yn = x * yn  +  (num_pows[i] != -1 ? param[c--] : bigfloat(0l));  

  c = n+d;
  bigfloat yd = 1l; //Highest degree coefficient is 1.0
  for (int i = pow_d-1; i >= 0; i--) yd = x * yd  +  (den_pows[i] != -1 ? param[c--] : bigfloat(0l)); 

  return(yn/yd);
}

// Compute size and sign of the approximation error at x
bigfloat AlgRemezGeneral::getErr(bigfloat x, int *sign) const{
  bigfloat f = func(x);
  bigfloat e = approx(x) - f;
  if (f != 0) e /= f;
  if (e < (bigfloat)0.0) {
    *sign = -1;
    e = -e;
  }
  else *sign = 1;
  
  return(e);
}

// Solve the system AX=B
int AlgRemezGeneral::simq(){

  int ip, ipj, ipk, ipn;
  int idxpiv;
  int kp, kp1, kpk, kpn;
  int nip, nkp;
  bigfloat em, q, rownrm, big, size, pivot, sum;
  bigfloat *aa;
  bigfloat *X = param.data();

  int n = neq;
  int nm1 = n - 1;
  // Initialize IPS and X
  
  int ij = 0;
  for (int i = 0; i < n; i++) {
    IPS[i] = i;
    rownrm = 0.0;
    for(int j = 0; j < n; j++) {
      q = abs_bf(A[ij]);
      if(rownrm < q) rownrm = q;
      ++ij;
    }
    if (rownrm == (bigfloat)0l) {
      std::cout<<"simq rownrm=0\n";
      return(1);
    }
    X[i] = (bigfloat)1.0 / rownrm;
  }
  
  for (int k = 0; k < nm1; k++) {
    big = 0.0;
    idxpiv = 0;
    
    for (int i = k; i < n; i++) {
      ip = IPS[i];
      ipk = n*ip + k;
      size = abs_bf(A[ipk]) * X[ip];
      if (size > big) {
	big = size;
	idxpiv = i;
      }
    }
    
    if (big == (bigfloat)0l) {
      std::cout<<"simq big=0\n";
      return(2);
    }
    if (idxpiv != k) {
      int j = IPS[k];
      IPS[k] = IPS[idxpiv];
      IPS[idxpiv] = j;
    }
    kp = IPS[k];
    kpk = n*kp + k;
    pivot = A[kpk];
    kp1 = k+1;
    for (int i = kp1; i < n; i++) {
      ip = IPS[i];
      ipk = n*ip + k;
      em = -A[ipk] / pivot;
      A[ipk] = -em;
      nip = n*ip;
      nkp = n*kp;
      aa = A.data()+nkp+kp1;
      for (int j = kp1; j < n; j++) {
	ipj = nip + j;
	A[ipj] = A[ipj] + em * *aa++;
      }
    }
  }
  kpn = n * IPS[n-1] + n - 1;	// last element of IPS[n] th row
  if (A[kpn] == (bigfloat)0l) {
    std::cout<<"simq A[kpn]=0\n";
    return(3);
  }

  
  ip = IPS[0];
  X[0] = B[ip];
  for (int i = 1; i < n; i++) {
    ip = IPS[i];
    ipj = n * ip;
    sum = 0.0;
    for (int j = 0; j < i; j++) {
      sum += A[ipj] * X[j];
      ++ipj;
    }
    X[i] = B[ip] - sum;
  }
  
  ipn = n * IPS[n-1] + n - 1;
  X[n-1] = X[n-1] / A[ipn];
  
  for (int iback = 1; iback < n; iback++) {
    //i goes (n-1),...,1
    int i = nm1 - iback;
    ip = IPS[i];
    nip = n*ip;
    sum = 0.0;
    aa = A.data()+nip+i+1;
    for (int j= i + 1; j < n; j++) 
      sum += *aa++ * X[j];
    X[i] = (X[i] - sum) / A[nip+i];
  }
  
  return(0);
}

void AlgRemezGeneral::csv(std::ostream & os) const{
  os << "Numerator" << std::endl;
  for(int i=0;i<=pow_n;i++){
    os << getCoeffNum(i) << "*x^" << i;
    if(i!=pow_n) os << " + ";
  }
  os << std::endl;

  os << "Denominator" << std::endl;
  for(int i=0;i<=pow_d;i++){
    os << getCoeffDen(i) << "*x^" << i;
    if(i!=pow_d) os << " + ";
  }
  os << std::endl;

  //For a true minimax solution the errors should all be equal and the signs should oscillate +-+-+- etc
  int sign;
  os << "Errors at maxima: coordinate, error, (sign)" << std::endl;
  for(int i=0;i<neq+1;i++){ 
    os << mm[i] << " " << getErr(mm[i],&sign) << " (" << sign << ")" << std::endl;
  }

  os << "Scan over range:" << std::endl;
  int npt = 60;
  bigfloat dlt = (apend - apstrt)/bigfloat(npt-1);

  for (bigfloat x=apstrt; x<=apend; x = x + dlt) {
    double f = evaluateFunc(x);
    double r = evaluateApprox(x);
    os<< x<<","<<r<<","<<f<<","<<r-f<<std::endl;
  }
  return;
}
