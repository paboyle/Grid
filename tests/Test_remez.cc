#include <Grid.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

class MultiShiftFunction {
public:
  std::vector<double> poles;
  std::vector<double> residues;
  double norm;
  double lo,hi;
  MultiShiftFunction(int n,double _lo,double _hi): poles(n), residues(n), lo(_lo), hi(_hi) {;};
  double approx(double x);
  void csv(std::ostream &out);
  void gnuplot(std::ostream &out);
};
double MultiShiftFunction::approx(double x)
{
  double a = norm;
  for(int n=0;n<poles.size();n++){
    a = a + residues[n]/(x+poles[n]);
  }
  return a;
}
void MultiShiftFunction::gnuplot(std::ostream &out)
{
  out<<"f(x) = "<<norm<<"";
  for(int n=0;n<poles.size();n++){
    out<<"+("<<residues[n]<<"/(x+"<<poles[n]<<"))";
  }
  out<<";"<<std::endl;
}
void MultiShiftFunction::csv(std::ostream &out)
{
  for (double x=lo; x<hi; x*=1.05) {
    double f = approx(x);
    double r = std::sqrt(x);
    out<< x<<","<<r<<","<<f<<","<<r-f<<std::endl;
  }
  return;
}

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  std::cout << "Testing Remez"<<std::endl;

  double     lo=0.01;
  double     hi=1.0;
  int precision=64;
  int    degree=10;
  AlgRemez remez(0.001,1.0,precision);

  ////////////////////////////////////////
  // sqrt and inverse sqrt
  ////////////////////////////////////////
  MultiShiftFunction Sqrt(degree,lo,hi);
  MultiShiftFunction InvSqrt(degree,lo,hi);

  MultiShiftFunction SqrtSqrt(degree,lo,hi);
  MultiShiftFunction InvSqrtSqrt(degree,lo,hi);


  std::cout << "Generating degree "<<degree<<" for x^(1/2)"<<std::endl;
  remez.generateApprox(degree,1,2);
  remez.getPFE (&   Sqrt.residues[0],&   Sqrt.poles[0],&   Sqrt.norm);
  remez.getIPFE(&InvSqrt.residues[0],&InvSqrt.poles[0],&InvSqrt.norm);

  std::cout << "Generating degree "<<degree<<" for x^(1/4)"<<std::endl;
  remez.generateApprox(degree,1,4);
  remez.getPFE (&SqrtSqrt.residues[0],&SqrtSqrt.poles[0],&SqrtSqrt.norm);
  remez.getIPFE(&InvSqrtSqrt.residues[0],&InvSqrtSqrt.poles[0],&InvSqrtSqrt.norm);
  
  ofstream gnuplot(std::string("Sqrt.gnu"),std::ios::out|std::ios::trunc);
  Sqrt.gnuplot(gnuplot);

  ofstream gnuplot_inv(std::string("InvSqrt.gnu"),std::ios::out|std::ios::trunc);
  InvSqrt.gnuplot(gnuplot);

  double x=0.6789;
  double sx=std::sqrt(x);
  double ssx=std::sqrt(sx);
  double isx=1.0/sx;
  double issx=1.0/ssx;

  double asx  =Sqrt.approx(x);
  double assx =SqrtSqrt.approx(x);
  double aisx =InvSqrt.approx(x);
  double aissx=InvSqrtSqrt.approx(x);

  std::cout << "x^(1/2) : "<<sx<<" "<<asx<<std::endl;
  std::cout << "x^(1/4) : "<<ssx<<" "<<assx<<std::endl;
  std::cout << "x^(-1/2): "<<isx<<" "<<aisx<<std::endl;
  std::cout << "x^(-1/4): "<<issx<<" "<<aissx<<std::endl;
  assert(fabs(sx-asx)<1.0e-6);
  assert(fabs(ssx-assx)<1.0e-6);
  assert(fabs(isx-aisx)<1.0e-6);
  assert(fabs(issx-aissx)<1.0e-6);

  Grid_finalize();
}
