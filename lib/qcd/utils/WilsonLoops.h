#ifndef QCD_UTILS_WILSON_LOOPS_H
#define QCD_UTILS_WILSON_LOOPS_H
namespace Grid {
namespace QCD {

// Common wilson loop observables
template<class GaugeMat,class GaugeLorentz>
class WilsonLoops {
public:

  //////////////////////////////////////////////////
  // directed plaquette oriented in mu,nu plane
  //////////////////////////////////////////////////
  static void dirPlaquette(GaugeMat &plaq,const std::vector<GaugeMat> &U, const int mu, const int nu)
  {
    plaq=CovShiftForward(U[mu],mu,U[nu])*adj(CovShiftForward(U[nu],nu,U[mu]));
  }
  //////////////////////////////////////////////////
  // trace of directed plaquette oriented in mu,nu plane
  //////////////////////////////////////////////////
  static void traceDirPlaquette(LatticeComplex &plaq, const std::vector<GaugeMat> &U, const int mu, const int nu)
  {
    GaugeMat sp(U[0]._grid);
    dirPlaquette(sp,U,mu,nu);
    plaq=trace(sp);
  }
  //////////////////////////////////////////////////
  // sum over all planes of plaquette
  //////////////////////////////////////////////////
  static void sitePlaquette(LatticeComplex &Plaq,const std::vector<GaugeMat> &U)
  {
    LatticeComplex sitePlaq(U[0]._grid);
    for(int mu=1;mu<Nd;mu++){
      for(int nu=0;nu<mu;nu++){
	traceDirPlaquette(sitePlaq,U,mu,nu);
	Plaq = Plaq + sitePlaq;
      }
    }
  }
  //////////////////////////////////////////////////
  // sum over all x,y,z,t and over all planes of plaquette
  //////////////////////////////////////////////////
  static RealD sumPlaquette(const GaugeLorentz &Umu){
    std::vector<GaugeMat> U(4,Umu._grid);
    for(int mu=0;mu<Nd;mu++){
      U[mu] = peekIndex<LorentzIndex>(Umu,mu);
    }

    LatticeComplex Plaq(Umu._grid);
    
    sitePlaquette(Plaq,U);
    
    TComplex Tp = sum(Plaq);
    Complex p  = TensorRemove(Tp);
    return p.real();
  }
  //////////////////////////////////////////////////
  // average over all x,y,z,t and over all planes of plaquette
  //////////////////////////////////////////////////
  static RealD avgPlaquette(const GaugeLorentz &Umu){

    RealD sumplaq = sumPlaquette(Umu);
    
    double vol = Umu._grid->gSites();
    
    double faces = (1.0*Nd*(Nd-1))/2.0;
    
    return sumplaq/vol/faces/Nc; // Nd , Nc dependent... FIXME
  }
  //////////////////////////////////////////////////
  // the sum over all staples on each site
  //////////////////////////////////////////////////
  static void Staple(GaugeMat &staple,GaugeLorentz &Umu,int mu){

    std::vector<GaugeMat> U(4,Umu._grid);
    for(int d=0;d<Nd;d++){
      U[d] = peekIndex<LorentzIndex>(Umu,d);
    }

    staple = zero;
    GaugeMat tmp(U[0]._grid);

    for(int nu=0;nu<Nd;nu++){

      if(nu != mu) {

      // mu
      // ^
      // |__  nu

      //    __                                         __
      //      |         |                          
      //    __|  =    __|                         *
      //

	staple   += CovShiftForward(U[nu],nu,U[mu])*Cshift(adj(U[nu]),mu,+1);

      //
      //  __        __                         
      // |         |                          
      // |__     = |                              *       __
      //
      //
	tmp    = CovShiftForward (U[mu],mu,U[nu]);
	staple+= CovShiftBackward(U[nu],nu,tmp);
      
      }
    }
  }

  //////////////////////////////////////////////////////
  // Similar to above for rectangle is required
  //////////////////////////////////////////////////////
  /*
void siteRectangle(GaugeMat &plaq,const std::vector<GaugeMat> &U, const int mu, const int nu){
RealD avgRectangle(const std::vector<GaugeMat> &U){}
RealD avgRectangle(const std::vector<GaugeMat> &U, const int mu, const int nu){}
void traceRectangle(LatticeComplex &plaq,const std::vector<GaugeMat> &U, const int mu, const int nu){}
void siteRectangle(GaugeMat &plaq,const std::vector<GaugeMat> &U, const int mu, const int nu){}
  */

};


 typedef WilsonLoops<LatticeColourMatrix,LatticeGaugeField> ColourWilsonLoops;

}}






#endif
