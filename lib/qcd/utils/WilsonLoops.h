    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/utils/WilsonLoops.h

    Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: neo <cossu@post.kek.jp>
Author: paboyle <paboyle@ph.ed.ac.uk>

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
#ifndef QCD_UTILS_WILSON_LOOPS_H
#define QCD_UTILS_WILSON_LOOPS_H
namespace Grid {
namespace QCD {

// Common wilson loop observables
template<class Gimpl>
class WilsonLoops : public Gimpl {
public:
  
  INHERIT_GIMPL_TYPES(Gimpl);

  typedef typename Gimpl::GaugeLinkField GaugeMat;
  typedef typename Gimpl::GaugeField     GaugeLorentz;

  //////////////////////////////////////////////////
  // directed plaquette oriented in mu,nu plane
  //////////////////////////////////////////////////
  static void dirPlaquette(GaugeMat &plaq,const std::vector<GaugeMat> &U, const int mu, const int nu)
  {
    // Annoyingly, must use either scope resolution to find dependent base class, 
    // or this-> ; there is no "this" in a static method. This forces explicit Gimpl scope
    // resolution throughout the usage in this file, and rather defeats the purpose of deriving
    // from Gimpl.
    plaq= Gimpl::CovShiftBackward(U[mu],mu,
	  Gimpl::CovShiftBackward(U[nu],nu,
          Gimpl::CovShiftForward (U[mu],mu,U[nu])));
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

      //    __ 
      //      |
      //    __|
      //

	staple+=Gimpl::ShiftStaple(
	        Gimpl::CovShiftForward (U[nu],nu, 
		Gimpl::CovShiftBackward(U[mu],mu,
		Gimpl::CovShiftIdentityBackward(U[nu],nu))),mu);

      //  __ 
      // |   
      // |__ 
      //
      //
	staple+=Gimpl::ShiftStaple(  
                Gimpl::CovShiftBackward(U[nu],nu,		  		  
		Gimpl::CovShiftBackward(U[mu],mu,U[nu])),mu);
      }
    }
  }

  //////////////////////////////////////////////////////
  // Similar to above for rectangle is required
  //////////////////////////////////////////////////////
  static void dirRectangle(GaugeMat &rect,const std::vector<GaugeMat> &U, const int mu, const int nu)
  {
    rect =  Gimpl::CovShiftForward(U[mu],mu,Gimpl::CovShiftForward(U[mu],mu,U[nu]))* // ->->|
	adj(Gimpl::CovShiftForward(U[nu],nu,Gimpl::CovShiftForward(U[mu],mu,U[mu]))) ;
    rect = rect + 
          Gimpl::CovShiftForward(U[mu],mu,Gimpl::CovShiftForward(U[nu],nu,U[nu]))* // ->||
      adj(Gimpl::CovShiftForward(U[nu],nu,Gimpl::CovShiftForward(U[nu],nu,U[mu]))) ;
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
  static void RectStapleDouble(GaugeMat &U2,const GaugeMat & U,int mu){
    U2 = U * Cshift(U,mu,1);
  }

  ////////////////////////////////////////////////////////////////////////////
  // Hop by two optimisation strategy does not work nicely with Gparity. (could do,
  // but need to track two deep where cross boundary and apply a conjugation).
  // Must differentiate this in Gimpl, and use Gimpl::isPeriodicGaugeField to do so .
  ////////////////////////////////////////////////////////////////////////////
  static void RectStapleOptimised(GaugeMat &Stap,std::vector<GaugeMat> &U2,std::vector<GaugeMat> &U,int mu){

    Stap = zero;

    GridBase *grid = U[0]._grid;

    GaugeMat Staple2x1 (grid);
    GaugeMat tmp (grid);

    for(int nu=0;nu<Nd;nu++){
      if ( nu!=mu) {

	// Up staple    ___ ___ 
	//             |       |
	tmp = Cshift(adj(U[nu]),nu,-1); 
	tmp = adj(U2[mu])*tmp;
	tmp = Cshift(tmp,mu,-2);

	Staple2x1 = Gimpl::CovShiftForward (U[nu],nu,tmp);


	// Down staple
	//             |___ ___|
	//
	tmp = adj(U2[mu])*U[nu];
	Staple2x1+= Gimpl::CovShiftBackward(U[nu],nu,Cshift(tmp,mu,-2));


	//              ___ ___
	//             |    ___|
	//             |___ ___|
	//

	Stap+= Cshift(Gimpl::CovShiftForward (U[mu],mu,Staple2x1),mu,1);

	//              ___ ___
	//             |___    |
	//             |___ ___|
	//

	//	tmp= Staple2x1* Cshift(U[mu],mu,-2);
	//	Stap+= Cshift(tmp,mu,1) ;
	Stap+= Cshift(Staple2x1,mu,1)*Cshift(U[mu],mu,-1); ;

	//       --    
	//      |  |              
	//          
	//      |  | 
	
	tmp = Cshift(adj(U2[nu]),nu,-2);
	tmp = Gimpl::CovShiftBackward(U[mu],mu,tmp);
	tmp = U2[nu]*Cshift(tmp,nu,2);
	Stap+= Cshift(tmp, mu, 1);

	//      |  |              
	//          
	//      |  | 
	//       -- 
	
	tmp = Gimpl::CovShiftBackward(U[mu],mu,U2[nu]);
	tmp = adj(U2[nu])*tmp;
	tmp = Cshift(tmp,nu,-2);
	Stap+=Cshift(tmp, mu, 1);
    }}


  }

  static void RectStaple(GaugeMat &Stap,const GaugeLorentz & Umu,int mu)
  {
    RectStapleUnoptimised(Stap,Umu,mu);
  }
  static void RectStaple(const GaugeLorentz & Umu,GaugeMat &Stap,
			 std::vector<GaugeMat> &U2,
			 std::vector<GaugeMat> &U, int mu)
  {
    if ( Gimpl::isPeriodicGaugeField() ){ 
      RectStapleOptimised(Stap,U2,U,mu);
    } else {
      RectStapleUnoptimised(Stap,Umu,mu);
    }
  }

  static void RectStapleUnoptimised(GaugeMat &Stap,const GaugeLorentz &Umu,int mu){
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
    Stap+= Gimpl::ShiftStaple(
		  Gimpl::CovShiftForward (U[mu],mu,
		  Gimpl::CovShiftForward (U[nu],nu,
		  Gimpl::CovShiftBackward(U[mu],mu,
                  Gimpl::CovShiftBackward(U[mu],mu,
		  Gimpl::CovShiftIdentityBackward(U[nu],nu))))) , mu);

    //              __ 
    //          |__ __ |

    Stap+= Gimpl::ShiftStaple(
                  Gimpl::CovShiftForward (U[mu],mu,
		  Gimpl::CovShiftBackward(U[nu],nu,
		  Gimpl::CovShiftBackward(U[mu],mu,
                  Gimpl::CovShiftBackward(U[mu],mu, U[nu])))) , mu);

    //           __ 
    //          |__ __ |

    Stap+= Gimpl::ShiftStaple(
		  Gimpl::CovShiftBackward(U[nu],nu,
		  Gimpl::CovShiftBackward(U[mu],mu,
		  Gimpl::CovShiftBackward(U[mu],mu,
		  Gimpl::CovShiftForward(U[nu],nu,U[mu])))) , mu);

    //           __ ___ 
    //          |__    |

    Stap+= Gimpl::ShiftStaple(
		   Gimpl::CovShiftForward (U[nu],nu,
	           Gimpl::CovShiftBackward(U[mu],mu,
                   Gimpl::CovShiftBackward(U[mu],mu,
                   Gimpl::CovShiftBackward(U[nu],nu,U[mu])))) , mu);

     //       --    
     //      |  |              
     //          
     //      |  | 
     
    Stap+= Gimpl::ShiftStaple(
		   Gimpl::CovShiftForward(U[nu],nu,
		   Gimpl::CovShiftForward(U[nu],nu,
                   Gimpl::CovShiftBackward(U[mu],mu,
                   Gimpl::CovShiftBackward(U[nu],nu,
		   Gimpl::CovShiftIdentityBackward(U[nu],nu))))) , mu);


     //      |  |              
     //          
     //      |  | 
     //       -- 
     
    Stap+= Gimpl::ShiftStaple(
		   Gimpl::CovShiftBackward(U[nu],nu,
		   Gimpl::CovShiftBackward(U[nu],nu,
                   Gimpl::CovShiftBackward(U[mu],mu,
                   Gimpl::CovShiftForward (U[nu],nu,U[nu])))) , mu);
    }}
  }


};


 typedef WilsonLoops<PeriodicGimplR> ColourWilsonLoops;
 typedef WilsonLoops<PeriodicGimplR> U1WilsonLoops;
 typedef WilsonLoops<PeriodicGimplR> SU2WilsonLoops;
 typedef WilsonLoops<PeriodicGimplR> SU3WilsonLoops;

}}

#endif
