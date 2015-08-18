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

  void Init(AlgRemez & remez,double tol,bool inverse) 
  {
    order=remez.getDegree();
    tolerances.resize(remez.getDegree(),tol);
    poles.resize(remez.getDegree());
    residues.resize(remez.getDegree());
    remez.getBounds(lo,hi);
    if ( inverse ) remez.getIPFE (&residues[0],&poles[0],&norm);
    else           remez.getPFE (&residues[0],&poles[0],&norm);
  }
  MultiShiftFunction(AlgRemez & remez,double tol,bool inverse)
  {
    Init(remez,tol,inverse);
  }

};
}
#endif
