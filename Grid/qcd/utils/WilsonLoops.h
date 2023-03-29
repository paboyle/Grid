/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./lib/qcd/utils/WilsonLoops.h

    Copyright (C) 2015

    Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
    Author: Peter Boyle <paboyle@ph.ed.ac.uk>
    Author: neo <cossu@post.kek.jp>
    Author: paboyle <paboyle@ph.ed.ac.uk>
    Author: James Harrison <J.Harrison@soton.ac.uk>
    Author: Antonin Portelli <antonin.portelli@me.com>

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

    See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
/*  END LEGAL */
#ifndef QCD_UTILS_WILSON_LOOPS_H
#define QCD_UTILS_WILSON_LOOPS_H

NAMESPACE_BEGIN(Grid);

// Common wilson loop observables
template <class Gimpl> class WilsonLoops : public Gimpl {
public:
  INHERIT_GIMPL_TYPES(Gimpl);

  typedef typename Gimpl::GaugeLinkField GaugeMat;
  typedef typename Gimpl::GaugeField GaugeLorentz;

  //////////////////////////////////////////////////
  // directed plaquette oriented in mu,nu plane
  //////////////////////////////////////////////////
  static void dirPlaquette(GaugeMat &plaq, const std::vector<GaugeMat> &U,
                           const int mu, const int nu) {
    // Annoyingly, must use either scope resolution to find dependent base
    // class,
    // or this-> ; there is no "this" in a static method. This forces explicit
    // Gimpl scope
    // resolution throughout the usage in this file, and rather defeats the
    // purpose of deriving
    // from Gimpl.
    /*
      plaq = Gimpl::CovShiftBackward(
      U[mu], mu, Gimpl::CovShiftBackward(
      U[nu], nu, Gimpl::CovShiftForward(U[mu], mu, U[nu])));
    */
    // _
    //|< _|
    plaq = Gimpl::CovShiftForward(U[mu],mu,
				  Gimpl::CovShiftForward(U[nu],nu,
							 Gimpl::CovShiftBackward(U[mu],mu,
										 Gimpl::CovShiftIdentityBackward(U[nu], nu))));




  }
  //////////////////////////////////////////////////
  // trace of directed plaquette oriented in mu,nu plane
  //////////////////////////////////////////////////
  static void traceDirPlaquette(ComplexField &plaq,
                                const std::vector<GaugeMat> &U, const int mu,
                                const int nu) {
    GaugeMat sp(U[0].Grid());
    dirPlaquette(sp, U, mu, nu);
    plaq = trace(sp);
  }
  //////////////////////////////////////////////////
  // sum over all planes of plaquette
  //////////////////////////////////////////////////
  static void sitePlaquette(ComplexField &Plaq,
                            const std::vector<GaugeMat> &U) {
    ComplexField sitePlaq(U[0].Grid());
    Plaq = Zero();
    for (int mu = 1; mu < Nd; mu++) {
      for (int nu = 0; nu < mu; nu++) {
        traceDirPlaquette(sitePlaq, U, mu, nu);
        Plaq = Plaq + sitePlaq;
      }
    }
  }
  //////////////////////////////////////////////////
  // sum over all x,y,z,t and over all planes of plaquette
  //////////////////////////////////////////////////
  static RealD sumPlaquette(const GaugeLorentz &Umu) {
    std::vector<GaugeMat> U(Nd, Umu.Grid());
    // inefficient here
    for (int mu = 0; mu < Nd; mu++) {
      U[mu] = PeekIndex<LorentzIndex>(Umu, mu);
    }

    ComplexField Plaq(Umu.Grid());

    sitePlaquette(Plaq, U);
    auto Tp = sum(Plaq);
    auto p = TensorRemove(Tp);
    return p.real();
  }


  //////////////////////////////////////////////////
  // average over all x,y,z,t and over all planes of plaquette
  //////////////////////////////////////////////////
  static RealD avgPlaquette(const GaugeLorentz &Umu) {
    RealD sumplaq = sumPlaquette(Umu);
    double vol = Umu.Grid()->gSites();
    double faces = (1.0 * Nd * (Nd - 1)) / 2.0;
    return sumplaq / vol / faces / Nc; // Nd , Nc dependent... FIXME
  }

  //////////////////////////////////////////////////
  // sum over all spatial planes of plaquette
  //////////////////////////////////////////////////
  static void siteSpatialPlaquette(ComplexField &Plaq,
                            const std::vector<GaugeMat> &U) {
    ComplexField sitePlaq(U[0].Grid());
    Plaq = Zero();
    for (int mu = 1; mu < Nd-1; mu++) {
      for (int nu = 0; nu < mu; nu++) {
        traceDirPlaquette(sitePlaq, U, mu, nu);
        Plaq = Plaq + sitePlaq;
      }
    }
  }

  ////////////////////////////////////
  // sum over all x,y,z and over all spatial planes of plaquette
  //////////////////////////////////////////////////
  static std::vector<RealD> timesliceSumSpatialPlaquette(const GaugeLorentz &Umu) {
    std::vector<GaugeMat> U(Nd, Umu.Grid());
    // inefficient here
    for (int mu = 0; mu < Nd; mu++) {
      U[mu] = PeekIndex<LorentzIndex>(Umu, mu);
    }

    ComplexField Plaq(Umu.Grid());

    siteSpatialPlaquette(Plaq, U);
    typedef typename ComplexField::scalar_object sobj;
    std::vector<sobj> Tq;
    sliceSum(Plaq, Tq, Nd-1);

    std::vector<Real> out(Tq.size());
    for(int t=0;t<Tq.size();t++) out[t] = TensorRemove(Tq[t]).real();
    return out;
  }
  
  //////////////////////////////////////////////////
  // average over all x,y,z and over all spatial planes of plaquette
  //////////////////////////////////////////////////
  static std::vector<RealD> timesliceAvgSpatialPlaquette(const GaugeLorentz &Umu) {
    std::vector<RealD> sumplaq = timesliceSumSpatialPlaquette(Umu);
    int Lt = Umu.Grid()->FullDimensions()[Nd-1];
    assert(sumplaq.size() == Lt);
    double vol = Umu.Grid()->gSites() / Lt;
    double faces = (1.0 * (Nd - 1)* (Nd - 2)) / 2.0;
    for(int t=0;t<Lt;t++)
      sumplaq[t] = sumplaq[t] / vol / faces / Nc; // Nd , Nc dependent... FIXME
    return sumplaq;
  }

  //////////////////////////////////////////////////
  // average over all x,y,z the temporal loop
  //////////////////////////////////////////////////
  static ComplexD avgPolyakovLoop(const GaugeField &Umu) {  //assume Nd=4
    GaugeMat Ut(Umu.Grid()), P(Umu.Grid());
    ComplexD out;
    int T = Umu.Grid()->GlobalDimensions()[3];
    int X = Umu.Grid()->GlobalDimensions()[0];
    int Y = Umu.Grid()->GlobalDimensions()[1];
    int Z = Umu.Grid()->GlobalDimensions()[2];

    Ut = peekLorentz(Umu,3); //Select temporal direction
    P = Ut;
    for (int t=1;t<T;t++){ 
      P = Gimpl::CovShiftForward(Ut,3,P);
    }
   RealD norm = 1.0/(Nc*X*Y*Z*T);
   out = sum(trace(P))*norm;
   return out;   
}

  //////////////////////////////////////////////////
  // average over traced single links
  //////////////////////////////////////////////////
  static RealD linkTrace(const GaugeLorentz &Umu) {
    std::vector<GaugeMat> U(Nd, Umu.Grid());

    ComplexField Tr(Umu.Grid());
    Tr = Zero();
    for (int mu = 0; mu < Nd; mu++) {
      U[mu] = PeekIndex<LorentzIndex>(Umu, mu);
      Tr = Tr + trace(U[mu]);
    }

    auto Tp = sum(Tr);
    auto p = TensorRemove(Tp);

    double vol = Umu.Grid()->gSites();

    return p.real() / vol / (4.0 * Nc ) ;
  };

  //////////////////////////////////////////////////
  // the sum over all staples on each site in direction mu,nu
  //////////////////////////////////////////////////
  static void Staple(GaugeMat &staple, const GaugeLorentz &Umu, int mu,
                     int nu) {

    GridBase *grid = Umu.Grid();

    std::vector<GaugeMat> U(Nd, grid);
    for (int d = 0; d < Nd; d++) {
      U[d] = PeekIndex<LorentzIndex>(Umu, d);
    }
    staple = Zero();

    if (nu != mu) {

      // mu
      // ^
      // |__>  nu

      //    __
      //      |
      //    __|
      //

      staple += Gimpl::ShiftStaple(
				   Gimpl::CovShiftForward(
							  U[nu], nu,
							  Gimpl::CovShiftBackward(
										  U[mu], mu, Gimpl::CovShiftIdentityBackward(U[nu], nu))),
				   mu);

      //  __
      // |
      // |__
      //
      //
      staple += Gimpl::ShiftStaple(
				   Gimpl::CovShiftBackward(U[nu], nu,
							   Gimpl::CovShiftBackward(U[mu], mu, U[nu])),
				   mu);
    }
  }


  // For the force term
/*
  static void StapleMult(GaugeMat &staple, const GaugeLorentz &Umu, int mu) {
    GridBase *grid = Umu.Grid();
    std::vector<GaugeMat> U(Nd, grid);
    for (int d = 0; d < Nd; d++) {
      // this operation is taking too much time
      U[d] = PeekIndex<LorentzIndex>(Umu, d);
    }
    staple = Zero();
    GaugeMat tmp1(grid);
    GaugeMat tmp2(grid);

    for (int nu = 0; nu < Nd; nu++) {
      if (nu != mu) {
        // this is ~10% faster than the Staple  -- PAB: so what it gives the WRONG answers for other BC's!
        tmp1 = Cshift(U[nu], mu, 1);
        tmp2 = Cshift(U[mu], nu, 1);
        staple += tmp1* adj(U[nu]*tmp2);
        tmp2 = adj(U[mu]*tmp1)*U[nu];
        staple += Cshift(tmp2, nu, -1);
      }
    }
    staple = U[mu]*staple;
  }
*/
  //////////////////////////////////////////////////
  // the sum over all staples on each site
  //////////////////////////////////////////////////
  static void Staple(GaugeMat &staple, const GaugeLorentz &Umu, int mu) {

    GridBase *grid = Umu.Grid();

    std::vector<GaugeMat> U(Nd, grid);
    for (int d = 0; d < Nd; d++) {
      U[d] = PeekIndex<LorentzIndex>(Umu, d);
    }
    staple = Zero();

    for (int nu = 0; nu < Nd; nu++) {

      if (nu != mu) {

        // mu
        // ^
        // |__>  nu

        //    __
        //      |
        //    __|
        //
     
        staple += Gimpl::ShiftStaple(
				     Gimpl::CovShiftForward(
							    U[nu], nu,
							    Gimpl::CovShiftBackward(
										    U[mu], mu, Gimpl::CovShiftIdentityBackward(U[nu], nu))),
				     mu);

        //  __
        // |
        // |__
        //
        //

        staple += Gimpl::ShiftStaple(
				     Gimpl::CovShiftBackward(U[nu], nu,
							     Gimpl::CovShiftBackward(U[mu], mu, U[nu])), mu);
      }
    }
  }

  //////////////////////////////////////////////////
  // the sum over all staples on each site in direction mu,nu, upper part
  //////////////////////////////////////////////////
  static void StapleUpper(GaugeMat &staple, const GaugeLorentz &Umu, int mu,
                          int nu) {
    if (nu != mu) {
      GridBase *grid = Umu.Grid();

      std::vector<GaugeMat> U(Nd, grid);
      for (int d = 0; d < Nd; d++) {
        U[d] = PeekIndex<LorentzIndex>(Umu, d);// some redundant copies
      }

      // mu
      // ^
      // |__>  nu

      //    __
      //      |
      //    __|
      //

      staple = Gimpl::ShiftStaple(
				  Gimpl::CovShiftForward(
							 U[nu], nu,
							 Gimpl::CovShiftBackward(
										 U[mu], mu, Gimpl::CovShiftIdentityBackward(U[nu], nu))),
				  mu);
    }
  }

  ////////////////////////////////////////////////////////////////////////
  // the sum over all staples on each site in direction mu,nu, lower part
  ////////////////////////////////////////////////////////////////////////
  static void StapleLower(GaugeMat &staple, const GaugeLorentz &Umu, int mu,
                          int nu) {
    if (nu != mu) {
      GridBase *grid = Umu.Grid();

      std::vector<GaugeMat> U(Nd, grid);
      for (int d = 0; d < Nd; d++) {
        U[d] = PeekIndex<LorentzIndex>(Umu, d);// some redundant copies
      }

      // mu
      // ^
      // |__>  nu

      //  __
      // |
      // |__
      //
      //
      staple = Gimpl::ShiftStaple(
				  Gimpl::CovShiftBackward(U[nu], nu,
                                  Gimpl::CovShiftBackward(U[mu], mu, U[nu])),
          mu);

    }
  }

  //////////////////////////////////////////////////////
  //  Field Strength
  //////////////////////////////////////////////////////
  static void FieldStrength(GaugeMat &FS, const GaugeLorentz &Umu, int mu, int nu){
    // Fmn +--<--+  Ut +--<--+
    //     |     |     |     |
      //  (x)+-->--+     +-->--+(x)  - h.c.
    //     |     |     |     |
    //     +--<--+     +--<--+

    GaugeMat Vup(Umu.Grid()), Vdn(Umu.Grid());
    StapleUpper(Vup, Umu, mu, nu);
    StapleLower(Vdn, Umu, mu, nu);
    GaugeMat v = Vup - Vdn;
    GaugeMat u = PeekIndex<LorentzIndex>(Umu, mu);  // some redundant copies
    GaugeMat vu = v*u;
      //FS = 0.25*Ta(u*v + Cshift(vu, mu, -1));
      FS = (u*v + Gimpl::CshiftLink(vu, mu, -1));
      FS = 0.125*(FS - adj(FS));
  }

  static Real TopologicalCharge(const GaugeLorentz &U){
    // 4d topological charge
    assert(Nd==4);
    // Bx = -iF(y,z), By = -iF(z,y), Bz = -iF(x,y)
    GaugeMat Bx(U.Grid()), By(U.Grid()), Bz(U.Grid());
    FieldStrength(Bx, U, Ydir, Zdir);
    FieldStrength(By, U, Zdir, Xdir);
    FieldStrength(Bz, U, Xdir, Ydir);

    // Ex = -iF(t,x), Ey = -iF(t,y), Ez = -iF(t,z)
    GaugeMat Ex(U.Grid()), Ey(U.Grid()), Ez(U.Grid());
    FieldStrength(Ex, U, Tdir, Xdir);
    FieldStrength(Ey, U, Tdir, Ydir);
    FieldStrength(Ez, U, Tdir, Zdir);

    double coeff = 8.0/(32.0*M_PI*M_PI);

    ComplexField qfield = coeff*trace(Bx*Ex + By*Ey + Bz*Ez);
    auto Tq = sum(qfield);
    return TensorRemove(Tq).real();
  }


  //Clover-leaf Wilson loop combination for arbitrary mu-extent M and nu extent N,  mu >= nu
  //cf  https://arxiv.org/pdf/hep-lat/9701012.pdf Eq 7  for 1x2 Wilson loop    
  //Clockwise ordering
  static void CloverleafMxN(GaugeMat &FS, const GaugeMat &Umu, const GaugeMat &Unu, int mu, int nu, int M, int N){  
#define Fmu(A) Gimpl::CovShiftForward(Umu, mu, A)
#define Bmu(A) Gimpl::CovShiftBackward(Umu, mu, A)
#define Fnu(A) Gimpl::CovShiftForward(Unu, nu, A)
#define Bnu(A) Gimpl::CovShiftBackward(Unu, nu, A)
#define FmuI Gimpl::CovShiftIdentityForward(Umu, mu)
#define BmuI Gimpl::CovShiftIdentityBackward(Umu, mu)
#define FnuI Gimpl::CovShiftIdentityForward(Unu, nu)
#define BnuI Gimpl::CovShiftIdentityBackward(Unu, nu)

    //Upper right loop
    GaugeMat tmp = BmuI;
    for(int i=1;i<M;i++)
      tmp = Bmu(tmp);
    for(int j=0;j<N;j++)
      tmp = Bnu(tmp);
    for(int i=0;i<M;i++)
      tmp = Fmu(tmp);
    for(int j=0;j<N;j++)
      tmp = Fnu(tmp);
      
    FS = tmp;

    //Upper left loop
    tmp = BnuI;
    for(int j=1;j<N;j++)
      tmp = Bnu(tmp);
    for(int i=0;i<M;i++)
      tmp = Fmu(tmp);
    for(int j=0;j<N;j++)
      tmp = Fnu(tmp);
    for(int i=0;i<M;i++)
      tmp = Bmu(tmp);
      
    FS = FS + tmp;

    //Lower right loop
    tmp = FnuI;
    for(int j=1;j<N;j++)
      tmp = Fnu(tmp);
    for(int i=0;i<M;i++)
      tmp = Bmu(tmp);
    for(int j=0;j<N;j++)
      tmp = Bnu(tmp);
    for(int i=0;i<M;i++)
      tmp = Fmu(tmp);
      
    FS = FS + tmp;

    //Lower left loop
    tmp = FmuI;
    for(int i=1;i<M;i++)
      tmp = Fmu(tmp);
    for(int j=0;j<N;j++)
      tmp = Fnu(tmp);
    for(int i=0;i<M;i++)
      tmp = Bmu(tmp);
    for(int j=0;j<N;j++)
      tmp = Bnu(tmp);

    FS = FS + tmp;

#undef Fmu
#undef Bmu
#undef Fnu
#undef Bnu
#undef FmuI
#undef BmuI
#undef FnuI
#undef BnuI
  }

  //Field strength from MxN Wilson loop
  //Note F_numu = - F_munu
  static void FieldStrengthMxN(GaugeMat &FS, const GaugeLorentz &U, int mu, int nu, int M, int N){  
    GaugeMat Umu = PeekIndex<LorentzIndex>(U, mu);
    GaugeMat Unu = PeekIndex<LorentzIndex>(U, nu);
    if(M == N){
      GaugeMat F(Umu.Grid());
      CloverleafMxN(F, Umu, Unu, mu, nu, M, N);
      FS = 0.125 * ( F - adj(F) );
    }else{
      //Average over both orientations
      GaugeMat horizontal(Umu.Grid()), vertical(Umu.Grid());
      CloverleafMxN(horizontal, Umu, Unu, mu, nu, M, N);
      CloverleafMxN(vertical, Umu, Unu, mu, nu, N, M);
      FS = 0.0625 * ( horizontal - adj(horizontal) + vertical - adj(vertical) );
    }
  }

  //Topological charge contribution from MxN Wilson loops
  //cf  https://arxiv.org/pdf/hep-lat/9701012.pdf  Eq 6
  //output is the charge by timeslice: sum over timeslices to obtain the total
  static std::vector<Real> TimesliceTopologicalChargeMxN(const GaugeLorentz &U, int M, int N){
    assert(Nd == 4);
    std::vector<std::vector<GaugeMat*> > F(Nd,std::vector<GaugeMat*>(Nd,nullptr));
    //Note F_numu = - F_munu
    //hence we only need to loop over mu,nu,rho,sigma that aren't related by permuting mu,nu  or rho,sigma
    //Use nu > mu
    for(int mu=0;mu<Nd-1;mu++){
      for(int nu=mu+1; nu<Nd; nu++){
	F[mu][nu] = new GaugeMat(U.Grid());
	FieldStrengthMxN(*F[mu][nu], U, mu, nu, M, N);
      }
    }
    Real coeff = -1./(32 * M_PI*M_PI * M*M * N*N); //overall sign to match CPS and Grid conventions, possibly related to time direction = 3 vs 0

    static const int combs[3][4] = { {0,1,2,3}, {0,2,1,3}, {0,3,1,2} };
    static const int signs[3] = { 1, -1, 1 }; //epsilon_{mu nu rho sigma}

    ComplexField fsum(U.Grid());
    fsum = Zero();
    for(int c=0;c<3;c++){
      int mu = combs[c][0], nu = combs[c][1], rho = combs[c][2], sigma = combs[c][3];
      int eps = signs[c];
      fsum = fsum + (8. * coeff * eps) * trace( (*F[mu][nu]) * (*F[rho][sigma]) ); 
    }

    for(int mu=0;mu<Nd-1;mu++)
      for(int nu=mu+1; nu<Nd; nu++)
	delete F[mu][nu];
    
    typedef typename ComplexField::scalar_object sobj;
    std::vector<sobj> Tq;
    sliceSum(fsum, Tq, Nd-1);

    std::vector<Real> out(Tq.size());
    for(int t=0;t<Tq.size();t++) out[t] = TensorRemove(Tq[t]).real();
    return out;
  }
  static Real TopologicalChargeMxN(const GaugeLorentz &U, int M, int N){
    std::vector<Real> Tq = TimesliceTopologicalChargeMxN(U,M,N);
    Real out(0);
    for(int t=0;t<Tq.size();t++) out += Tq[t];
    return out;
  }

  //Generate the contributions to the 5Li topological charge from Wilson loops of the following sizes
  //Use coefficients from hep-lat/9701012
  //1x1 : c1=(19.-55.*c5)/9.
  //2x2 : c2=(1-64.*c5)/9.
  //1x2 : c3=(-64.+640.*c5)/45.
  //1x3 : c4=1./5.-2.*c5
  //3x3 : c5=1./20.
  //Output array outer index contains the loops in the above order
  //Inner index is the time coordinate
  static std::vector<std::vector<Real> > TimesliceTopologicalCharge5LiContributions(const GaugeLorentz &U){
    static const int exts[5][2] = { {1,1}, {2,2}, {1,2}, {1,3}, {3,3} };       
    std::vector<std::vector<Real> > out(5);
    for(int i=0;i<5;i++){	
      out[i] = TimesliceTopologicalChargeMxN(U,exts[i][0],exts[i][1]);
    }
    return out;
  }   

  static std::vector<Real> TopologicalCharge5LiContributions(const GaugeLorentz &U){   
    static const int exts[5][2] = { {1,1}, {2,2}, {1,2}, {1,3}, {3,3} };
    std::vector<Real> out(5);
    std::cout << GridLogMessage << "Computing topological charge" << std::endl;
    for(int i=0;i<5;i++){
      out[i] = TopologicalChargeMxN(U,exts[i][0],exts[i][1]);
      std::cout << GridLogMessage << exts[i][0] << "x" << exts[i][1] << " Wilson loop contribution " << out[i] << std::endl;
    }
    return out;
  }

  //Compute the 5Li topological charge
  static std::vector<Real> TimesliceTopologicalCharge5Li(const GaugeLorentz &U){
    std::vector<std::vector<Real> > loops = TimesliceTopologicalCharge5LiContributions(U);

    double c5=1./20.;
    double c4=1./5.-2.*c5;
    double c3=(-64.+640.*c5)/45.;
    double c2=(1-64.*c5)/9.;
    double c1=(19.-55.*c5)/9.;

    int Lt = loops[0].size();
    std::vector<Real> out(Lt,0.);
    for(int t=0;t<Lt;t++)
      out[t] += c1*loops[0][t] + c2*loops[1][t] + c3*loops[2][t] + c4*loops[3][t] + c5*loops[4][t];
    return out;
  }

  static Real TopologicalCharge5Li(const GaugeLorentz &U){
    std::vector<Real> Qt = TimesliceTopologicalCharge5Li(U);
    Real Q = 0.;
    for(int t=0;t<Qt.size();t++) Q += Qt[t];
    std::cout << GridLogMessage << "5Li Topological charge: " << Q << std::endl;
    return Q;
  }




  //////////////////////////////////////////////////////
  // Similar to above for rectangle is required
  //////////////////////////////////////////////////////
  static void dirRectangle(GaugeMat &rect, const std::vector<GaugeMat> &U,
                           const int mu, const int nu) {
    rect = Gimpl::CovShiftForward(
				  U[mu], mu, Gimpl::CovShiftForward(U[mu], mu, U[nu])) * // ->->|
      adj(Gimpl::CovShiftForward(
				 U[nu], nu, Gimpl::CovShiftForward(U[mu], mu, U[mu])));
    rect = rect +
      Gimpl::CovShiftForward(
			     U[mu], mu, Gimpl::CovShiftForward(U[nu], nu, U[nu])) * // ->||
      adj(Gimpl::CovShiftForward(
				 U[nu], nu, Gimpl::CovShiftForward(U[nu], nu, U[mu])));
  }
  static void traceDirRectangle(ComplexField &rect,
                                const std::vector<GaugeMat> &U, const int mu,
                                const int nu) {
    GaugeMat sp(U[0].Grid());
    dirRectangle(sp, U, mu, nu);
    rect = trace(sp);
  }
  static void siteRectangle(ComplexField &Rect,
                            const std::vector<GaugeMat> &U) {
    ComplexField siteRect(U[0].Grid());
    Rect = Zero();
    for (int mu = 1; mu < Nd; mu++) {
      for (int nu = 0; nu < mu; nu++) {
        traceDirRectangle(siteRect, U, mu, nu);
        Rect = Rect + siteRect;
      }
    }
  }

  //////////////////////////////////////////////////
  // sum over all x,y,z,t and over all planes of plaquette
  //////////////////////////////////////////////////
  static RealD sumRectangle(const GaugeLorentz &Umu) {
    std::vector<GaugeMat> U(Nd, Umu.Grid());

    for (int mu = 0; mu < Nd; mu++) {
      U[mu] = PeekIndex<LorentzIndex>(Umu, mu);
    }

    ComplexField Rect(Umu.Grid());

    siteRectangle(Rect, U);

    auto Tp = sum(Rect);
    auto p = TensorRemove(Tp);
    return p.real();
  }
  //////////////////////////////////////////////////
  // average over all x,y,z,t and over all planes of plaquette
  //////////////////////////////////////////////////
  static RealD avgRectangle(const GaugeLorentz &Umu) {

    RealD sumrect = sumRectangle(Umu);

    double vol = Umu.Grid()->gSites();

    double faces = (1.0 * Nd * (Nd - 1)); // 2 distinct orientations summed

    return sumrect / vol / faces / Nc; // Nd , Nc dependent... FIXME
  }

  //////////////////////////////////////////////////
  // the sum over all staples on each site
  //////////////////////////////////////////////////
  static void RectStapleDouble(GaugeMat &U2, const GaugeMat &U, int mu) {
    U2 = U * Cshift(U, mu, 1);
  }

  ////////////////////////////////////////////////////////////////////////////
  // Hop by two optimisation strategy does not work nicely with Gparity. (could
  // do,
  // but need to track two deep where cross boundary and apply a conjugation).
  // Must differentiate this in Gimpl, and use Gimpl::isPeriodicGaugeField to do
  // so .
  ////////////////////////////////////////////////////////////////////////////
  static void RectStapleOptimised(GaugeMat &Stap, std::vector<GaugeMat> &U2,
                                  std::vector<GaugeMat> &U, int mu) {

    Stap = Zero();

    GridBase *grid = U[0].Grid();

    GaugeMat Staple2x1(grid);
    GaugeMat tmp(grid);

    for (int nu = 0; nu < Nd; nu++) {
      if (nu != mu) {

        // Up staple    ___ ___
        //             |       |
        tmp = Cshift(adj(U[nu]), nu, -1);
        tmp = adj(U2[mu]) * tmp;
        tmp = Cshift(tmp, mu, -2);

        Staple2x1 = Gimpl::CovShiftForward(U[nu], nu, tmp);

        // Down staple
        //             |___ ___|
        //
        tmp = adj(U2[mu]) * U[nu];
        Staple2x1 += Gimpl::CovShiftBackward(U[nu], nu, Cshift(tmp, mu, -2));

        //              ___ ___
        //             |    ___|
        //             |___ ___|
        //

        Stap += Cshift(Gimpl::CovShiftForward(U[mu], mu, Staple2x1), mu, 1);

        //              ___ ___
        //             |___    |
        //             |___ ___|
        //

        //  tmp= Staple2x1* Cshift(U[mu],mu,-2);
        //  Stap+= Cshift(tmp,mu,1) ;
        Stap += Cshift(Staple2x1, mu, 1) * Cshift(U[mu], mu, -1);
        ;

        //       --
        //      |  |
        //
        //      |  |

        tmp = Cshift(adj(U2[nu]), nu, -2);
        tmp = Gimpl::CovShiftBackward(U[mu], mu, tmp);
        tmp = U2[nu] * Cshift(tmp, nu, 2);
        Stap += Cshift(tmp, mu, 1);

        //      |  |
        //
        //      |  |
        //       --

        tmp = Gimpl::CovShiftBackward(U[mu], mu, U2[nu]);
        tmp = adj(U2[nu]) * tmp;
        tmp = Cshift(tmp, nu, -2);
        Stap += Cshift(tmp, mu, 1);
      }
    }
  }

  static void RectStaple(GaugeMat &Stap, const GaugeLorentz &Umu, int mu) {
    RectStapleUnoptimised(Stap, Umu, mu);
  }
  static void RectStaple(const GaugeLorentz &Umu, GaugeMat &Stap,
                         std::vector<GaugeMat> &U2, std::vector<GaugeMat> &U,
                         int mu) {
    if (Gimpl::isPeriodicGaugeField()) {
      RectStapleOptimised(Stap, U2, U, mu);
    } else {
      RectStapleUnoptimised(Stap, Umu, mu);
    }
  }

  static void RectStapleUnoptimised(GaugeMat &Stap, const GaugeLorentz &Umu,
                                    int mu) {
    GridBase *grid = Umu.Grid();

    std::vector<GaugeMat> U(Nd, grid);
    for (int d = 0; d < Nd; d++) {
      U[d] = PeekIndex<LorentzIndex>(Umu, d);
    }

    Stap = Zero();

    for (int nu = 0; nu < Nd; nu++) {
      if (nu != mu) {
        //           __ ___
        //          |    __ |
        //
        Stap += Gimpl::ShiftStaple(
				   Gimpl::CovShiftForward(
							  U[mu], mu,
							  Gimpl::CovShiftForward(
										 U[nu], nu,
										 Gimpl::CovShiftBackward(
													 U[mu], mu,
													 Gimpl::CovShiftBackward(
																 U[mu], mu,
																 Gimpl::CovShiftIdentityBackward(U[nu], nu))))),
				   mu);

        //              __
        //          |__ __ |

        Stap += Gimpl::ShiftStaple(
				   Gimpl::CovShiftForward(
							  U[mu], mu,
							  Gimpl::CovShiftBackward(
										  U[nu], nu,
										  Gimpl::CovShiftBackward(
													  U[mu], mu, Gimpl::CovShiftBackward(U[mu], mu, U[nu])))),
				   mu);

        //           __
        //          |__ __ |

        Stap += Gimpl::ShiftStaple(
				   Gimpl::CovShiftBackward(
							   U[nu], nu,
							   Gimpl::CovShiftBackward(
										   U[mu], mu,
										   Gimpl::CovShiftBackward(
													   U[mu], mu, Gimpl::CovShiftForward(U[nu], nu, U[mu])))),
				   mu);

        //           __ ___
        //          |__    |

        Stap += Gimpl::ShiftStaple(
				   Gimpl::CovShiftForward(
							  U[nu], nu,
							  Gimpl::CovShiftBackward(
										  U[mu], mu,
										  Gimpl::CovShiftBackward(
													  U[mu], mu, Gimpl::CovShiftBackward(U[nu], nu, U[mu])))),
				   mu);

        //       --
        //      |  |
        //
        //      |  |

        Stap += Gimpl::ShiftStaple(
				   Gimpl::CovShiftForward(
							  U[nu], nu,
							  Gimpl::CovShiftForward(
										 U[nu], nu,
										 Gimpl::CovShiftBackward(
													 U[mu], mu,
													 Gimpl::CovShiftBackward(
																 U[nu], nu,
																 Gimpl::CovShiftIdentityBackward(U[nu], nu))))),
				   mu);

        //      |  |
        //
        //      |  |
        //       --

        Stap += Gimpl::ShiftStaple(
				   Gimpl::CovShiftBackward(
							   U[nu], nu,
							   Gimpl::CovShiftBackward(
										   U[nu], nu,
										   Gimpl::CovShiftBackward(
													   U[mu], mu, Gimpl::CovShiftForward(U[nu], nu, U[nu])))),
				   mu);
      }
    }
  }

  //////////////////////////////////////////////////
  // Wilson loop of size (R1, R2), oriented in mu,nu plane
  //////////////////////////////////////////////////
  static void wilsonLoop(GaugeMat &wl, const std::vector<GaugeMat> &U,
                           const int Rmu, const int Rnu,
                           const int mu, const int nu) {
    wl = U[nu];

    for(int i = 0; i < Rnu-1; i++){
      wl = Gimpl::CovShiftForward(U[nu], nu, wl);
    }

    for(int i = 0; i < Rmu; i++){
      wl = Gimpl::CovShiftForward(U[mu], mu, wl);
    }

    for(int i = 0; i < Rnu; i++){
      wl = Gimpl::CovShiftBackward(U[nu], nu, wl);
    }

    for(int i = 0; i < Rmu; i++){
      wl = Gimpl::CovShiftBackward(U[mu], mu, wl);
    }
  }
  //////////////////////////////////////////////////
  // trace of Wilson Loop oriented in mu,nu plane
  //////////////////////////////////////////////////
  static void traceWilsonLoop(LatticeComplex &wl,
                                const std::vector<GaugeMat> &U,
                                const int Rmu, const int Rnu,
                                const int mu, const int nu) {
    GaugeMat sp(U[0].Grid());
    wilsonLoop(sp, U, Rmu, Rnu, mu, nu);
    wl = trace(sp);
  }
  //////////////////////////////////////////////////
  // sum over all planes of Wilson loop
  //////////////////////////////////////////////////
  static void siteWilsonLoop(LatticeComplex &Wl,
                            const std::vector<GaugeMat> &U,
                            const int R1, const int R2) {
    LatticeComplex siteWl(U[0].Grid());
    Wl = Zero();
    for (int mu = 1; mu < U[0].Grid()->_ndimension; mu++) {
      for (int nu = 0; nu < mu; nu++) {
        traceWilsonLoop(siteWl, U, R1, R2, mu, nu);
        Wl = Wl + siteWl;
        traceWilsonLoop(siteWl, U, R2, R1, mu, nu);
        Wl = Wl + siteWl;
      }
    }
  }
  //////////////////////////////////////////////////
  // sum over planes of Wilson loop with length R1
  // in the time direction
  //////////////////////////////////////////////////
  static void siteTimelikeWilsonLoop(LatticeComplex &Wl,
                            const std::vector<GaugeMat> &U,
                            const int R1, const int R2) {
    LatticeComplex siteWl(U[0].Grid());

    int ndim = U[0].Grid()->_ndimension;

    Wl = Zero();
    for (int nu = 0; nu < ndim - 1; nu++) {
      traceWilsonLoop(siteWl, U, R1, R2, ndim-1, nu);
      Wl = Wl + siteWl;
    }
  }
  //////////////////////////////////////////////////
  // sum Wilson loop over all planes orthogonal to the time direction
  //////////////////////////////////////////////////
  static void siteSpatialWilsonLoop(LatticeComplex &Wl,
                            const std::vector<GaugeMat> &U,
                            const int R1, const int R2) {
    LatticeComplex siteWl(U[0].Grid());

    Wl = Zero();
    for (int mu = 1; mu < U[0].Grid()->_ndimension - 1; mu++) {
      for (int nu = 0; nu < mu; nu++) {
        traceWilsonLoop(siteWl, U, R1, R2, mu, nu);
        Wl = Wl + siteWl;
        traceWilsonLoop(siteWl, U, R2, R1, mu, nu);
        Wl = Wl + siteWl;
      }
    }
  }
  //////////////////////////////////////////////////
  // sum over all x,y,z,t and over all planes of Wilson loop
  //////////////////////////////////////////////////
  static Real sumWilsonLoop(const GaugeLorentz &Umu,
                            const int R1, const int R2) {
    std::vector<GaugeMat> U(4, Umu.Grid());

    for (int mu = 0; mu < Umu.Grid()->_ndimension; mu++) {
      U[mu] = PeekIndex<LorentzIndex>(Umu, mu);
    }

    LatticeComplex Wl(Umu.Grid());

    siteWilsonLoop(Wl, U, R1, R2);

    TComplex Tp = sum(Wl);
    Complex p = TensorRemove(Tp);
    return p.real();
  }
  //////////////////////////////////////////////////
  // sum over all x,y,z,t and over all planes of timelike Wilson loop
  //////////////////////////////////////////////////
  static Real sumTimelikeWilsonLoop(const GaugeLorentz &Umu,
                            const int R1, const int R2) {
    std::vector<GaugeMat> U(4, Umu.Grid());

    for (int mu = 0; mu < Umu.Grid()->_ndimension; mu++) {
      U[mu] = PeekIndex<LorentzIndex>(Umu, mu);
    }

    LatticeComplex Wl(Umu.Grid());

    siteTimelikeWilsonLoop(Wl, U, R1, R2);

    TComplex Tp = sum(Wl);
    Complex p = TensorRemove(Tp);
    return p.real();
  }
  //////////////////////////////////////////////////
  // sum over all x,y,z,t and over all planes of spatial Wilson loop
  //////////////////////////////////////////////////
  static Real sumSpatialWilsonLoop(const GaugeLorentz &Umu,
                            const int R1, const int R2) {
    std::vector<GaugeMat> U(4, Umu.Grid());

    for (int mu = 0; mu < Umu.Grid()->_ndimension; mu++) {
      U[mu] = PeekIndex<LorentzIndex>(Umu, mu);
    }

    LatticeComplex Wl(Umu.Grid());

    siteSpatialWilsonLoop(Wl, U, R1, R2);

    TComplex Tp = sum(Wl);
    Complex p = TensorRemove(Tp);
    return p.real();
  }
  //////////////////////////////////////////////////
  // average over all x,y,z,t and over all planes of Wilson loop
  //////////////////////////////////////////////////
  static Real avgWilsonLoop(const GaugeLorentz &Umu,
                            const int R1, const int R2) {
    int ndim = Umu.Grid()->_ndimension;
    Real sumWl = sumWilsonLoop(Umu, R1, R2);
    Real vol = Umu.Grid()->gSites();
    Real faces = 1.0 * ndim * (ndim - 1);
    return sumWl / vol / faces / Nc; // Nc dependent... FIXME
  }
  //////////////////////////////////////////////////
  // average over all x,y,z,t and over all planes of timelike Wilson loop
  //////////////////////////////////////////////////
  static Real avgTimelikeWilsonLoop(const GaugeLorentz &Umu,
                            const int R1, const int R2) {
    int ndim = Umu.Grid()->_ndimension;
    Real sumWl = sumTimelikeWilsonLoop(Umu, R1, R2);
    Real vol = Umu.Grid()->gSites();
    Real faces = 1.0 * (ndim - 1);
    return sumWl / vol / faces / Nc; // Nc dependent... FIXME
  }
  //////////////////////////////////////////////////
  // average over all x,y,z,t and over all planes of spatial Wilson loop
  //////////////////////////////////////////////////
  static Real avgSpatialWilsonLoop(const GaugeLorentz &Umu,
                            const int R1, const int R2) {
    int ndim = Umu.Grid()->_ndimension;
    Real sumWl = sumSpatialWilsonLoop(Umu, R1, R2);
    Real vol = Umu.Grid()->gSites();
    Real faces = 1.0 * (ndim - 1) * (ndim - 2);
    return sumWl / vol / faces / Nc; // Nc dependent... FIXME
  }
};

typedef WilsonLoops<PeriodicGimplR> ColourWilsonLoops;
typedef WilsonLoops<PeriodicGimplR> U1WilsonLoops;
typedef WilsonLoops<PeriodicGimplR> SU2WilsonLoops;
typedef WilsonLoops<PeriodicGimplR> SU3WilsonLoops;

NAMESPACE_END(Grid);

#endif
