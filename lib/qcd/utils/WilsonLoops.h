#ifndef QCD_UTILS_WILSON_LOOPS_H
#define QCD_UTILS_WILSON_LOOPS_H
namespace Grid {
namespace QCD {

// Common wilson loop observables
template<class GaugeLorentz>
class WilsonLoops {
public:

  typedef LorentzScalar<GaugeLorentz> GaugeMat;

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
    Plaq=zero;
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
      U[mu] = PeekIndex<LorentzIndex>(Umu,mu);
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
  static RealD linkTrace(const GaugeLorentz &Umu){
    std::vector<GaugeMat> U(4,Umu._grid);

    LatticeComplex Tr(Umu._grid); Tr=zero;
    for(int mu=0;mu<Nd;mu++){
      U[mu] = PeekIndex<LorentzIndex>(Umu,mu);
      Tr = Tr+trace(U[mu]);
    }
    
    TComplex Tp = sum(Tr);
    Complex p  = TensorRemove(Tp);

    double vol = Umu._grid->gSites();

    return p.real()/vol/4.0/3.0;
  };
  //////////////////////////////////////////////////
  // the sum over all staples on each site
  //////////////////////////////////////////////////
  static void Staple(GaugeMat &staple,const GaugeLorentz &Umu,int mu){

    GridBase *grid = Umu._grid;

    std::vector<GaugeMat> U(4,grid);
    for(int d=0;d<Nd;d++){
      U[d] = PeekIndex<LorentzIndex>(Umu,d);
    }

    staple = zero;
    GaugeMat tmp(grid);

    for(int nu=0;nu<Nd;nu++){

      if(nu != mu) {

      // mu
      // ^
      // |__>  nu

      //    __                                         __
      //      |         |                          
      //    __|  =    __|                         *
      //

	//	staple   += CovShiftForward(U[nu],nu,U[mu])*Cshift(adj(U[nu]),mu,+1);
	staple   += 
	  Cshift(CovShiftForward (U[nu],nu, 
		 CovShiftBackward(U[mu],mu,Cshift(adj(U[nu]),nu,-1))),mu,1);


	// Unu(x) Umu(x+nu) Unu^dag(x+mu) ; left mult by Umu^dag(x) to close ring

      //
      //  __        __                         
      // |         |                          
      // |__     = |                              *       __
      //
      //
	staple   += 
          Cshift(
	  CovShiftBackward(U[nu],nu,		  		  
	  CovShiftBackward(U[mu],mu,U[nu])),mu,1);

      //	tmp    = CovShiftForward (U[mu],mu,U[nu]);
      //	staple+= CovShiftBackward(U[nu],nu,tmp);

	// Unu^dag(x-nu) Umu(x-nu) Unu(x+mu-nu) ; left mult by Umu^dag(x) to close ring.
      
      }
    }
  }

  //////////////////////////////////////////////////////
  // Similar to above for rectangle is required
  //////////////////////////////////////////////////////
  static void dirRectangle(GaugeMat &rect,const std::vector<GaugeMat> &U, const int mu, const int nu)
  {
    rect = CovShiftForward(U[mu],mu,CovShiftForward(U[mu],mu,U[nu]))* // ->->|
	adj(CovShiftForward(U[nu],nu,CovShiftForward(U[mu],mu,U[mu]))) ;
    rect = rect + 
      CovShiftForward(U[mu],mu,CovShiftForward(U[nu],nu,U[nu]))* // ->||
      adj(CovShiftForward(U[nu],nu,CovShiftForward(U[nu],nu,U[mu]))) ;
  }
  static void traceDirRectangle(LatticeComplex &rect, const std::vector<GaugeMat> &U, const int mu, const int nu)
  {
    GaugeMat sp(U[0]._grid);
    dirRectangle(sp,U,mu,nu);
    rect=trace(sp);
  }
  static void siteRectangle(LatticeComplex &Rect,const std::vector<GaugeMat> &U)
  {
    LatticeComplex siteRect(U[0]._grid);
    Rect=zero;
    for(int mu=1;mu<Nd;mu++){
      for(int nu=0;nu<mu;nu++){
	traceDirRectangle(siteRect,U,mu,nu);
	Rect = Rect + siteRect;
      }
    }
  }
  //////////////////////////////////////////////////
  // sum over all x,y,z,t and over all planes of plaquette
  //////////////////////////////////////////////////
  static RealD sumRectangle(const GaugeLorentz &Umu){
    std::vector<GaugeMat> U(4,Umu._grid);

    for(int mu=0;mu<Nd;mu++){
      U[mu] = PeekIndex<LorentzIndex>(Umu,mu);
    }

    LatticeComplex Rect(Umu._grid);
    
    siteRectangle(Rect,U);
    
    TComplex Tp = sum(Rect);
    Complex p  = TensorRemove(Tp);
    return p.real();
  }
  //////////////////////////////////////////////////
  // average over all x,y,z,t and over all planes of plaquette
  //////////////////////////////////////////////////
  static RealD avgRectangle(const GaugeLorentz &Umu){

    RealD sumrect = sumRectangle(Umu);
    
    double vol = Umu._grid->gSites();
    
    double faces = (1.0*Nd*(Nd-1)); // 2 distinct orientations summed
    
    return sumrect/vol/faces/Nc; // Nd , Nc dependent... FIXME
  }

  //////////////////////////////////////////////////
  // the sum over all staples on each site
  //////////////////////////////////////////////////
  static void RectStaple(GaugeMat &Stap,const GaugeLorentz &Umu,int mu){

    GridBase *grid = Umu._grid;

    std::vector<GaugeMat> U(4,grid);
    for(int d=0;d<Nd;d++){
      U[d] = PeekIndex<LorentzIndex>(Umu,d);
    }

    Stap=zero;

    for(int nu=0;nu<Nd;nu++){
      if ( nu!=mu) {
    //           __ ___ 
    //          |    __ |
    //
    Stap+= Cshift(CovShiftForward (U[mu],mu,
		  CovShiftForward (U[nu],nu,
		  CovShiftBackward(U[mu],mu,
                  CovShiftBackward(U[mu],mu,
	          Cshift(adj(U[nu]),nu,-1))))) , mu, 1);

    //              __ 
    //          |__ __ |

    Stap+= Cshift(CovShiftForward (U[mu],mu,
		  CovShiftBackward(U[nu],nu,
		  CovShiftBackward(U[mu],mu,
                  CovShiftBackward(U[mu],mu, U[nu])))) , mu, 1);

    //           __ 
    //          |__ __ |

    Stap+= Cshift(CovShiftBackward(U[nu],nu,
		  CovShiftBackward(U[mu],mu,
		  CovShiftBackward(U[mu],mu,
		  CovShiftForward(U[nu],nu,U[mu])))) , mu, 1);

    //           __ ___ 
    //          |__    |

     Stap+= Cshift(CovShiftForward (U[nu],nu,
	           CovShiftBackward(U[mu],mu,
                   CovShiftBackward(U[mu],mu,
                   CovShiftBackward(U[nu],nu,U[mu])))) , mu, 1);

     //       --    
     //      |  |              
     //          
     //      |  | 

     Stap+= Cshift(CovShiftForward(U[nu],nu,
		   CovShiftForward(U[nu],nu,
                   CovShiftBackward(U[mu],mu,
                   CovShiftBackward(U[nu],nu,
		   Cshift(adj(U[nu]),nu,-1))))) , mu, 1);


     //      |  |              
     //          
     //      |  | 
     //       -- 
     
     Stap+= Cshift(CovShiftBackward(U[nu],nu,
		   CovShiftBackward(U[nu],nu,
                   CovShiftBackward(U[mu],mu,
                   CovShiftForward (U[nu],nu,U[nu])))) , mu, 1);
    }}
  }


};


 typedef WilsonLoops<LatticeGaugeField> ColourWilsonLoops;
 typedef WilsonLoops<LatticeGaugeField> U1WilsonLoops;
 typedef WilsonLoops<LatticeGaugeField> SU2WilsonLoops;
 typedef WilsonLoops<LatticeGaugeField> SU3WilsonLoops;

}}

#endif
