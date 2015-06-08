#include <Grid.h>

namespace Grid {
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
    double r = sqrt(x);
    out<< x<<","<<r<<","<<f<<","<<r-f<<std::endl;
  }
  return;
}
}
