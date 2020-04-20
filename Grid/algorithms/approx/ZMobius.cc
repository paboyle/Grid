/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/algorithms/approx/ZMobius.cc

    Copyright (C) 2015

Author: Christopher Kelly <ckelly@phys.columbia.edu>

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

#include <Grid/algorithms/approx/ZMobius.h>
#include <Grid/algorithms/approx/RemezGeneral.h>

NAMESPACE_BEGIN(Grid);
NAMESPACE_BEGIN(Approx);

//Compute the tanh approximation
inline double epsilonMobius(const double x, const std::vector<ComplexD> &w){
  int Ls = w.size();

  ComplexD fxp = 1., fmp = 1.;
  for(int i=0;i<Ls;i++){
    fxp = fxp * ( w[i] + x );
    fmp = fmp * ( w[i] - x );
  }
  return ((fxp - fmp)/(fxp + fmp)).real();
}
inline double epsilonMobius(const double x, const std::vector<RealD> &w){
  int Ls = w.size();

  double fxp = 1., fmp = 1.;
  for(int i=0;i<Ls;i++){
    fxp = fxp * ( w[i] + x );
    fmp = fmp * ( w[i] - x );
  }
  return (fxp - fmp)/(fxp + fmp);
}



//Compute the tanh approximation in a form suitable for the Remez
bigfloat epsilonMobius(bigfloat x, void* data){
  const std::vector<RealD> &omega = *( (std::vector<RealD> const*)data );
  bigfloat fxp(1.0);
  bigfloat fmp(1.0);

  for(int i=0;i<omega.size();i++){
    fxp = fxp * ( bigfloat(omega[i]) + x);
    fmp = fmp * ( bigfloat(omega[i]) - x);
  }
  return (fxp - fmp)/(fxp + fmp);
}

//Compute the Zmobius Omega parameters suitable for eigenvalue range   -lambda_bound <= lambda <= lambda_bound
//Note omega_i = 1/(b_i + c_i)   where b_i and c_i are the Mobius parameters
void computeZmobiusOmega(std::vector<ComplexD> &omega_out, const int Ls_out,
			 const std::vector<RealD> &omega_in, const int Ls_in,
			 const RealD lambda_bound){
  assert(omega_in.size() == Ls_in);
  omega_out.resize(Ls_out);

  //Use the Remez algorithm to generate the appropriate rational polynomial
  //For odd polynomial, to satisfy Haar condition must take either positive or negative half of range (cf https://arxiv.org/pdf/0803.0439.pdf page 6)  
  AlgRemezGeneral remez(0, lambda_bound, 64, &epsilonMobius, (void*)&omega_in); 
  remez.generateApprox(Ls_out-1, Ls_out,AlgRemezGeneral::Odd, AlgRemezGeneral::Even, 1e-15, 100);
  remez.csv(std::cout);

  //The rational approximation has the form  [ f(x) - f(-x) ] / [ f(x) + f(-x) ]  where  f(x) = \Prod_{i=0}^{L_s-1} ( \omega_i + x )
  //cf https://academiccommons.columbia.edu/doi/10.7916/D8T72HD7  pg 102
  //omega_i are therefore the negative of the complex roots of f(x)

  //We can find the roots by recognizing that the eigenvalues of a matrix A are the roots of the characteristic polynomial
  // \rho(\lambda) = det( A - \lambda I )    where I is the unit matrix
  //The matrix whose characteristic polynomial is an arbitrary monic polynomial a0 + a1 x + a2 x^2 + ... x^n   is the companion matrix 
  // A = | 0    1   0    0 0 .... 0 |
  //     | 0    0   1    0 0 .... 0 |
  //     | :    :   :    : :      : |
  //     | 0    0   0    0 0      1
  //     | -a0 -a1 -a2  ...  ... -an|


  //Note the Remez defines the largest power to have unit coefficient
  std::vector<RealD> coeffs(Ls_out+1);
  for(int i=0;i<Ls_out+1;i+=2) coeffs[i] = coeffs[i] = remez.getCoeffDen(i); //even powers
  for(int i=1;i<Ls_out+1;i+=2) coeffs[i] = coeffs[i] = remez.getCoeffNum(i); //odd powers

  std::vector<std::complex<RealD> > roots(Ls_out);

  //Form the companion matrix
  Eigen::MatrixXd compn(Ls_out,Ls_out);
  for(int i=0;i<Ls_out-1;i++) compn(i,0) = 0.;
  compn(Ls_out - 1, 0) = -coeffs[0];
  
  for(int j=1;j<Ls_out;j++){
    for(int i=0;i<Ls_out-1;i++) compn(i,j) = i == j-1 ? 1. : 0.;
    compn(Ls_out - 1, j) = -coeffs[j];
  }

  //Eigensolve
  Eigen::EigenSolver<Eigen::MatrixXd> slv(compn, false);

  const auto & ev = slv.eigenvalues();
  for(int i=0;i<Ls_out;i++)
    omega_out[i] = -ev(i);

  //Sort ascending (smallest at start of vector!)
  std::sort(omega_out.begin(), omega_out.end(), 
	    [&](const ComplexD &a, const ComplexD &b){ return a.real() < b.real() || (a.real() == b.real() && a.imag() < b.imag()); });

  //McGlynn thesis pg 122 suggest improved iteration counts if magnitude of omega diminishes towards the center of the 5th dimension
  std::vector<ComplexD> omega_tmp = omega_out;
  int s_low=0, s_high=Ls_out-1, ss=0;
  for(int s_from = Ls_out-1; s_from >= 0; s_from--){ //loop from largest omega
    int s_to;
    if(ss % 2 == 0){
      s_to = s_low++;
    }else{
      s_to = s_high--;
    }
    omega_out[s_to] = omega_tmp[s_from];
    ++ss;
  }
  
  std::cout << "Resulting omega_i:" << std::endl;  
  for(int i=0;i<Ls_out;i++)
    std::cout << omega_out[i] << std::endl;

  std::cout << "Test result matches the approximate polynomial found by the Remez" << std::endl;
  std::cout << "<x> <remez approx> <poly approx> <diff poly approx remez approx> <exact> <diff poly approx exact>\n";
  
  int npt = 60;
  double dlt = lambda_bound/double(npt-1);

  for (int i =0; i<npt; i++){
    double x = i*dlt;
    double r = remez.evaluateApprox(x);
    double p = epsilonMobius(x, omega_out);
    double e = epsilonMobius(x, omega_in);

    std::cout << x<< " " << r << " " << p <<" " <<r-p << " " << e << " " << e-p << std::endl;
  }

}
  
//mobius_param = b+c   with b-c=1
void computeZmobiusOmega(std::vector<ComplexD> &omega_out, const int Ls_out, const RealD mobius_param, const int Ls_in, const RealD lambda_bound){
  std::vector<RealD> omega_in(Ls_in, 1./mobius_param);
  computeZmobiusOmega(omega_out, Ls_out, omega_in, Ls_in, lambda_bound);
}

//ZMobius class takes  gamma_i = (b+c) omega_i as its input, where b, c are factored out
void computeZmobiusGamma(std::vector<ComplexD> &gamma_out, 
			 const RealD mobius_param_out, const int Ls_out, 
			 const RealD mobius_param_in, const int Ls_in,
			 const RealD lambda_bound){
  computeZmobiusOmega(gamma_out, Ls_out, mobius_param_in, Ls_in, lambda_bound);
  for(int i=0;i<Ls_out;i++) gamma_out[i] = gamma_out[i] * mobius_param_out;
}
//Assumes mobius_param_out == mobius_param_in
void computeZmobiusGamma(std::vector<ComplexD> &gamma_out, const int Ls_out, const RealD mobius_param, const int Ls_in, const RealD lambda_bound){
  computeZmobiusGamma(gamma_out, mobius_param, Ls_out, mobius_param, Ls_in, lambda_bound);
}

NAMESPACE_END(Approx);
NAMESPACE_END(Grid);
