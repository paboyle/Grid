#ifndef MULTI_SHIFT_FUNCTION
#define MULTI_SHIFT_FUNCTION
namespace Grid {
class MultiShiftFunction {
public:
  int order;
  std::vector<RealD> poles;
  std::vector<RealD> residues;
  std::vector<RealD> tolerances;
  RealD norm;
  RealD lo,hi;
  MultiShiftFunction(int n,RealD _lo,RealD _hi): poles(n), residues(n), lo(_lo), hi(_hi) {;};
  RealD approx(RealD x);
  void csv(std::ostream &out);
  void gnuplot(std::ostream &out);
  MultiShiftFunction(AlgRemez & remez,double tol,bool inverse) :
      order(remez.getDegree()),
      tolerances(remez.getDegree(),tol),
      poles(remez.getDegree()),
      residues(remez.getDegree())
  {
    remez.getBounds(lo,hi);
    if ( inverse ) remez.getIPFE (&residues[0],&poles[0],&norm);
    else remez.getPFE (&residues[0],&poles[0],&norm);
  }
};
}
#endif
