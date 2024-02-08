    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/core/Test_gpdwf_Xconj.cc

    Copyright (C) 2015

Author: Christopher Kelly <ckelly@bnl.gov>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>

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

/**
 * Tests for the implementation of X-conjugate BCs
 *
 */

#include <Grid/Grid.h>

using namespace std;
using namespace Grid;

//A test implementation of the X-conjugate action as a wrapper around the regular 2f G-parity action
template<typename GPaction>
class XconjugateDWF{
public:
  typedef typename GPaction::FermionField GPfermion;
  typedef LatticeFermionD FermionField;
  GPaction *action;
  bool Xbar;
  XconjugateDWF(GPaction *action, bool Xbar = false): action(action), Xbar(Xbar){}

  template<typename Op>
  LatticeFermionD op11(const LatticeFermionD &in, const Op &op){
    GPfermion tmp1(in.Grid()), tmp2(in.Grid());
    tmp1.Checkerboard() = in.Checkerboard();
    tmp1 = Zero();
    PokeIndex<GparityFlavourIndex>(tmp1, in, 0);
    op(tmp2, tmp1);
    return PeekIndex<GparityFlavourIndex>(tmp2,0);
  }
  template<typename Op>
  LatticeFermionD op12(const LatticeFermionD &in, const Op &op){
    GPfermion tmp1(in.Grid()), tmp2(in.Grid());
    tmp1.Checkerboard() = in.Checkerboard();
    tmp1 = Zero();
    PokeIndex<GparityFlavourIndex>(tmp1, in, 1);
    op(tmp2, tmp1);
    return PeekIndex<GparityFlavourIndex>(tmp2,0);
  }
  
  template<typename Op>
  LatticeFermionD opFull(const LatticeFermionD &in, const Op &op){
    static Gamma C = Gamma(Gamma::Algebra::MinusGammaY) * Gamma(Gamma::Algebra::GammaT);
    static Gamma g5 = Gamma(Gamma::Algebra::Gamma5);
    static Gamma X = C*g5;
    LatticeFermionD tmp1(in.Grid());
    tmp1.Checkerboard() = in.Checkerboard();
    tmp1 = -(X*conjugate(in));
    if(Xbar) tmp1 = -tmp1;
    
    LatticeFermionD out11 = op11(in, op);
    LatticeFermionD out12 = op12(tmp1, op);
    LatticeFermionD out = out11 + out12;
    return out;
  }

#define DEFOP(OP)  void OP(const FermionField &in, FermionField &out){ out=opFull(in, [&](GPfermion &out, const GPfermion &in){ return action->OP(in,out); }); }
  DEFOP(M);
  DEFOP(Mdag);
  DEFOP(Meooe);
  DEFOP(MeooeDag);
  DEFOP(Mooee);
  DEFOP(MooeeDag);
  DEFOP(MooeeInv);
  DEFOP(MooeeInvDag);
#undef DEFOP;
};

//A 2-flavor representation of the X-conjugate action acting on X-conjugate fermion fields
template<typename XconjAction, typename GPAction>
class Xconjugate2fWrapper{
public:
  typedef typename XconjAction::FermionField OneFlavorFermionField;
  typedef Lattice<iGparitySpinColourVector<typename OneFlavorFermionField::vector_type> > FermionField;

  XconjAction *action;
  GPAction *gaction;

  Xconjugate2fWrapper(XconjAction *action, GPAction *gaction): action(action), gaction(gaction){}
  
  template<typename Op>
  void opFull(FermionField &out, const FermionField &in, const Op &op){
    static Gamma C = Gamma(Gamma::Algebra::MinusGammaY) * Gamma(Gamma::Algebra::GammaT);
    static Gamma g5 = Gamma(Gamma::Algebra::Gamma5);
    static Gamma X = C*g5;

    OneFlavorFermionField tmp(in.Grid());
    op(tmp, PeekIndex<GparityFlavourIndex>(in,0));
    
    out.Checkerboard() = tmp.Checkerboard();
    PokeIndex<GparityFlavourIndex>(out, tmp, 0);
    tmp = -(X*conjugate(tmp));
    tmp = tmp * action->Params.boundary_phase;
    PokeIndex<GparityFlavourIndex>(out, tmp, 1);
  }

#define DEFOP(OP)  void OP(const FermionField &in, FermionField &out){ opFull(out, in, \
									      [&](OneFlavorFermionField &out1f, \
										  const OneFlavorFermionField &in1f){ \
										return action->OP(in1f,out1f); \
									      }	\
									      ); \
                                                                      } 
  DEFOP(M);
  DEFOP(Mdag);
  DEFOP(Meooe);
  DEFOP(MeooeDag);
  DEFOP(Mooee);
  DEFOP(MooeeDag);
  DEFOP(MooeeInv);
  DEFOP(MooeeInvDag);
#undef DEFOP;
};

const Gamma & Xmatrix(){
  static Gamma C = Gamma(Gamma::Algebra::MinusGammaY) * Gamma(Gamma::Algebra::GammaT);
  static Gamma g5 = Gamma(Gamma::Algebra::Gamma5);
  static Gamma X = C*g5;
  return X;
}

typedef typename GparityMobiusFermionD::FermionField FermionField2f;
typedef typename XconjugateMobiusFermionD::FermionField FermionField1f;


void boostXconjToGparity(FermionField2f &out, const FermionField1f &in){
  PokeIndex<GparityFlavourIndex>(out, in, 0);
  FermionField1f tmp = -(Xmatrix()*conjugate(in));
  PokeIndex<GparityFlavourIndex>(out, tmp, 1);
}
void boostXbarConjToGparity(FermionField2f &out, const FermionField1f &in){
  PokeIndex<GparityFlavourIndex>(out, in, 0);
  FermionField1f tmp = Xmatrix()*conjugate(in);
  PokeIndex<GparityFlavourIndex>(out, tmp, 1);
}

template<typename Field>
inline RealD norm2Diff(const Field &a, const Field &b){
  return norm2(Field(a-b));
}

void norm2DiffXconj(double &diff_f0, double &diff_f1, const FermionField2f &gp, const FermionField1f &xconj){
  FermionField1f tmp = PeekIndex<GparityFlavourIndex>(gp,0);
  diff_f0 = norm2(FermionField1f(tmp-xconj));
  tmp = PeekIndex<GparityFlavourIndex>(gp,1);
  FermionField1f xconj_f1 = -(Xmatrix()*conjugate(xconj));
  diff_f1 = norm2(FermionField1f(tmp-xconj_f1));
}
void norm2DiffXbarConj(double &diff_f0, double &diff_f1, const FermionField2f &gp, const FermionField1f &xconj){
  FermionField1f tmp = PeekIndex<GparityFlavourIndex>(gp,0);
  diff_f0 = norm2(FermionField1f(tmp-xconj));
  tmp = PeekIndex<GparityFlavourIndex>(gp,1);
  FermionField1f xconj_f1 = Xmatrix()*conjugate(xconj);
  diff_f1 = norm2(FermionField1f(tmp-xconj_f1));
}

//1/sqrt(2) * 
//|  X   -1 |
//| -i   iX |
void applyUdag(FermionField2f &out, const FermionField2f &in){
  FermionField1f u = PeekIndex<GparityFlavourIndex>(in,0);
  FermionField1f l = PeekIndex<GparityFlavourIndex>(in,1);
  FermionField1f tmp = ComplexD(1./sqrt(2),0) * (Xmatrix()*u - l);
  PokeIndex<GparityFlavourIndex>(out, tmp, 0);
  tmp = ComplexD(0,-1./sqrt(2))*(u - Xmatrix()*l);
  PokeIndex<GparityFlavourIndex>(out, tmp, 1);
}
//1/sqrt(2) * 
//|  -X   i |
//| -1   iX |
void applyU(FermionField2f &out, const FermionField2f &in){
  FermionField1f u = PeekIndex<GparityFlavourIndex>(in,0);
  FermionField1f l = PeekIndex<GparityFlavourIndex>(in,1);
  FermionField1f tmp = ComplexD(-1./sqrt(2)) * (Xmatrix()*u + ComplexD(0,-1)*l);
  PokeIndex<GparityFlavourIndex>(out, tmp, 0);
  tmp = ComplexD(-1./sqrt(2))*(u + ComplexD(0,-1)*(Xmatrix()*l));
  PokeIndex<GparityFlavourIndex>(out, tmp, 1);
}
template<typename GPAction>
void applyR(FermionField2f &out, const FermionField2f &in, GPAction &action){
  FermionField2f tmp(in.Grid()), tmp2(in.Grid());
  applyU(tmp, in);
  action.M(tmp, tmp2);
  applyUdag(out,tmp2); //Rv
}
template<typename GPAction>
void applyRij(FermionField1f &out, const FermionField1f &in, const int i, const int j, GPAction &action){
  FermionField2f tmp2f(in.Grid());
  tmp2f = Zero();
  PokeIndex<GparityFlavourIndex>(tmp2f, in, j);
  FermionField2f tmp2f_2(in.Grid());
  applyR(tmp2f_2, tmp2f, action);
  out = PeekIndex<GparityFlavourIndex>(tmp2f_2,i);
}
template<typename GPAction>
void applyMij(FermionField1f &out, const FermionField1f &in, const int i, const int j, GPAction &action){
  FermionField2f tmp2f(in.Grid());
  tmp2f = Zero();
  PokeIndex<GparityFlavourIndex>(tmp2f, in, j);
  FermionField2f tmp2f_2(in.Grid());
  action.M(tmp2f,tmp2f_2);
  out = PeekIndex<GparityFlavourIndex>(tmp2f_2,i);
}
template<typename GPAction>
void applyMijstar(FermionField1f &out, const FermionField1f &in, const int i, const int j, GPAction &action){
  FermionField1f tmp(in.Grid()), tmp2(in.Grid()); //Mij* v = (Mij v*)*
  tmp = conjugate(in);
  applyMij(tmp2,tmp,i,j,action);
  out = conjugate(tmp2);
}


//A reference implementation of the real, rotated GP Dirac operator
template<typename GPAction>
class RotatedGPrefAction{
public:
  typedef typename GPAction::FermionField TwoFlavorFermionField;
  typedef Lattice<iSpinColourVector<typename TwoFlavorFermionField::vector_type> > OneFlavorFermionField;

  GPAction *gaction;

  RotatedGPrefAction(GPAction *gaction): gaction(gaction){}

  //out = D_ij in
  template<typename Op>
  void DD(OneFlavorFermionField &out, const OneFlavorFermionField &in, const int i, const int j, const Op &op){
    TwoFlavorFermionField tmp2f(in.Grid());
    tmp2f = Zero();
    PokeIndex<GparityFlavourIndex>(tmp2f, in, j);
    TwoFlavorFermionField tmp2f_2(in.Grid());
    op(tmp2f_2, tmp2f);
    out = PeekIndex<GparityFlavourIndex>(tmp2f_2,i);
  }
  //out = D*_ij in
  template<typename Op>
  void DDstar(OneFlavorFermionField &out, const OneFlavorFermionField &in, const int i, const int j, const Op &op){
    OneFlavorFermionField tmp(in.Grid()), tmp2(in.Grid());
    tmp = conjugate(in);
    DD(tmp2,tmp,i,j,op);
    out = conjugate(tmp2);
  }

  //Use Re(U) V = 0.5 ( U + U^*) V = 0.5 ( UV + [U V*]* )
  template<typename U>
  void do_real(OneFlavorFermionField &out, const OneFlavorFermionField &in, const U &u){
    OneFlavorFermionField tmp(in.Grid()), tmp2(in.Grid());
    u(out,in);
    tmp = conjugate(in);
    u(tmp2,tmp);
    tmp = conjugate(tmp2);
    out = out + tmp;
    out = out * 0.5;
  }
  //Use Im(U) V = -0.5i ( U - U^*) V = -0.5i ( UV - [U V*]* )
  template<typename U>
  void do_imag(OneFlavorFermionField &out, const OneFlavorFermionField &in, const U &u){
    OneFlavorFermionField tmp(in.Grid()), tmp2(in.Grid());
    u(out,in);

    tmp = conjugate(in);
    u(tmp2,tmp);
    tmp = conjugate(tmp2);
    out = out - tmp;
    out = out * ComplexD(0,-0.5);
  }

  //-(X D_11 X + X D_12) in
  template<typename Op>
  void opAparen(OneFlavorFermionField &out, const OneFlavorFermionField &in, const Op &op){
    OneFlavorFermionField tmp(in.Grid()), tmp2(in.Grid());
    tmp = Xmatrix()*in;
    DD(tmp2,tmp,0,0,op);
    out = Xmatrix()*tmp2;
    
    DD(tmp,in,0,1,op);
    tmp2 = Xmatrix() * tmp;
    out = out + tmp2;
    out = -out;
  }
  template<typename Op>
  void opA(OneFlavorFermionField &out, const OneFlavorFermionField &in, const Op &op){
    do_real(out,in, [&](OneFlavorFermionField &out, const OneFlavorFermionField &in){ opAparen(out,in,op); } );
  }


  //(-X D11 -X D12 X) in
  template<typename Op>
  void opBparen(OneFlavorFermionField &out, const OneFlavorFermionField &in, const Op &op){
    OneFlavorFermionField tmp(in.Grid()), tmp2(in.Grid());
    DD(tmp,in,0,0,op);
    out = Xmatrix()*tmp;
    
    tmp = Xmatrix()*in;
    DD(tmp2,tmp,0,1,op);
    tmp = Xmatrix() * tmp2;
    out = out + tmp;

    out = -out;
  }
  template<typename Op>
  void opB(OneFlavorFermionField &out, const OneFlavorFermionField &in, const Op &op){
    do_imag(out,in, [&](OneFlavorFermionField &out, const OneFlavorFermionField &in){ opBparen(out,in,op); } );
  }

  //-(D11 X + D12 ) in
  template<typename Op>
  void opCparen(OneFlavorFermionField &out, const OneFlavorFermionField &in, const Op &op){
    OneFlavorFermionField tmp(in.Grid()), tmp2(in.Grid());
    tmp = Xmatrix()*in;
    DD(out,tmp,0,0,op);
    
    DD(tmp,in,0,1,op);
    out = out + tmp;

    out = -out;
  }
  template<typename Op>
  void opC(OneFlavorFermionField &out, const OneFlavorFermionField &in, const Op &op){
    do_imag(out,in, [&](OneFlavorFermionField &out, const OneFlavorFermionField &in){ opCparen(out,in,op); } );
  }

  //(D11 + D12 X ) in
  template<typename Op>
  void opDparen(OneFlavorFermionField &out, const OneFlavorFermionField &in, const Op &op){
    OneFlavorFermionField tmp(in.Grid()), tmp2(in.Grid());
    DD(out,in,0,0,op);
    
    tmp = Xmatrix()*in;
    DD(tmp2,tmp,0,1,op);
    out = out + tmp2;
  }
  template<typename Op>
  void opD(OneFlavorFermionField &out, const OneFlavorFermionField &in, const Op &op){
    do_real(out,in, [&](OneFlavorFermionField &out, const OneFlavorFermionField &in){ opDparen(out,in,op); } );
  }
  
  template<typename Op>
  void opFull(TwoFlavorFermionField &out, const TwoFlavorFermionField &in, const Op &op){
    OneFlavorFermionField u(in.Grid()), l(in.Grid()), tmp(in.Grid()), tmp2(in.Grid());
    u = PeekIndex<GparityFlavourIndex>(in,0);
    l = PeekIndex<GparityFlavourIndex>(in,1);
    opA(tmp,u,op);
    opB(tmp2,l,op);
    tmp = tmp + tmp2;
    PokeIndex<GparityFlavourIndex>(out, tmp, 0);

    opC(tmp,u,op);
    opD(tmp2,l,op);
    tmp = tmp + tmp2;
    PokeIndex<GparityFlavourIndex>(out, tmp, 1);
  }

#define DEFOP(OP)  void OP(const TwoFlavorFermionField &in, TwoFlavorFermionField &out){ opFull(out, in, \
												[&](TwoFlavorFermionField &out, \
												    const TwoFlavorFermionField &in){ \
												  return gaction->OP(in,out); \
												} \
												); \
                                                                                        } 
  DEFOP(M);
  DEFOP(Mdag);
  DEFOP(Meooe);
  DEFOP(MeooeDag);
  DEFOP(Mooee);
  DEFOP(MooeeDag);
  DEFOP(MooeeInv);
  DEFOP(MooeeInvDag);
#undef DEFOP;
};


template<typename RefAction>
void applyRijRef(FermionField1f &out, const FermionField1f &in, const int i, const int j, RefAction &action){
  FermionField2f tmp2f(in.Grid());
  tmp2f = Zero();
  PokeIndex<GparityFlavourIndex>(tmp2f, in, j);
  FermionField2f tmp2f_2(in.Grid());
  action.M(tmp2f,tmp2f_2);
  out = PeekIndex<GparityFlavourIndex>(tmp2f_2,i);
}

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  const int Ls=8;

  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});
  GridParallelRNG          RNG5(FGrid);  RNG5.SeedFixedIntegers(seeds5);
  GridParallelRNG          RNG4(UGrid);  RNG4.SeedFixedIntegers(seeds4);
  GridParallelRNG          RNG5rb(FrbGrid);  RNG5.SeedFixedIntegers(seeds5);

  LatticeGaugeField Umu(UGrid); 
  SU<Nc>::HotConfiguration(RNG4, Umu);

  Gamma X = Xmatrix();
  
  //Set up a regular MDWF action instance as well as X-conj and Xbar-conj versions
  RealD mass=0.01;
  RealD M5=1.8;
  RealD mob_b=1.5;
  GparityMobiusFermionD ::ImplParams params;
  std::vector<int> twists({1,1,1,0});
  params.twists = twists;
    
  GparityMobiusFermionD reg_action(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,mob_b,mob_b-1.,params);

  XconjugateMobiusFermionD::ImplParams xparams;
  xparams.twists = twists;
  xparams.boundary_phase = 1.0;
  
  XconjugateMobiusFermionD xconj_action(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,mob_b,mob_b-1.,xparams);

  xparams.boundary_phase = -1.0;
  XconjugateMobiusFermionD xbarconj_action(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,mob_b,mob_b-1.,xparams);
 

  //#######################################################################################################################################
  {
    std::cout << "Testing M on full grid" << std::endl;

    FermionField1f rand_sc(FGrid), rand_sc2(FGrid);
    gaussian(RNG5,rand_sc);   
    gaussian(RNG5,rand_sc2);

    FermionField2f rand_f0(FGrid); //v \delta_f,0
    rand_f0 = Zero();    
    PokeIndex<GparityFlavourIndex>(rand_f0, rand_sc, 0);

    FermionField2f rand_f1(FGrid); //v \delta_f,1
    rand_f1 = Zero();
    PokeIndex<GparityFlavourIndex>(rand_f1, rand_sc, 1);
    
    FermionField2f tmp(FGrid), tmp2(FGrid), tmp3(FGrid), tmp4(FGrid), tmp5(FGrid);
    FermionField1f tmpsc(FGrid), tmpsc2(FGrid), tmpsc3(FGrid), tmpsc4(FGrid);
    RealD nrm;

    std::cout << "Test the relationship between the upper and lower rows in flavor space of the G-parity Dirac operator" << std::endl;

    //Test  M00 v =  -X M11* X v = -X [M11 X v*]*
    reg_action.M(rand_f0,tmp);
    FermionField1f M00v = PeekIndex<GparityFlavourIndex>(tmp,0);
    FermionField1f M10v = PeekIndex<GparityFlavourIndex>(tmp,1);

    tmp = X*conjugate(rand_f1);    
    reg_action.M(tmp,tmp2);
    FermionField1f M11Xvconj = PeekIndex<GparityFlavourIndex>(tmp2,1);
    FermionField1f M01Xvconj = PeekIndex<GparityFlavourIndex>(tmp2,0);

    tmpsc = -(X*conjugate(M11Xvconj));
    nrm = norm2Diff(tmpsc, M00v);
    std::cout << "Test of M00 v =  -X M11* X v = -X [M11 X v*]* (expect 0): " << nrm << std::endl;
    assert(nrm < 1e-10);

    //Test  M10 v =  X M01* X v = X [M01 X v*]*
    tmpsc = X*conjugate(M01Xvconj);
    nrm = norm2Diff(tmpsc, M10v);
    std::cout << "Test of M10 v =  X M11* X v = X [M11 X v*]*  (expect 0): " << nrm << std::endl;
    assert(nrm < 1e-10);

    //Test the X-conjugate implementation
    std::cout << "Test the implementation of the X-conjugate action against the reference"<< std::endl;
    XconjugateDWF<GparityMobiusFermionD> xconj_action_ref(&reg_action);
    XconjugateDWF<GparityMobiusFermionD> xbarconj_action_ref(&reg_action,true);

    //Test upper boundary
    std::vector<int> L(4);
    for(int mu=0;mu<4;mu++)
      L[mu] = UGrid->GlobalDimensions()[mu];
    Coordinate site(5,0);
    typedef FermionField1f::vector_object SpinorV;
    typedef typename SpinorV::scalar_object SpinorS;

    tmpsc3 = Zero();
    for(int i=2;i<5;i++) site[i] = L[i-1]/2; //midpoint in y,z,t
    site[1] = 0; //field only on site on lower boundary

    SpinorS v;
    peekSite(v, rand_sc, site);    
    pokeSite(v, tmpsc3, site);

    xconj_action_ref.M(tmpsc3, tmpsc);
    xconj_action.M(tmpsc3, tmpsc2);
    nrm = norm2Diff(tmpsc, tmpsc2);
    std::cout << "Test of Xconjugate action M upper boundary only (expect 0): " << nrm << std::endl;
    assert(nrm < 1e-10);

    site[1] = L[0]-1;

    SpinorS r1,r2;
    peekSite(r1, tmpsc, site);
    peekSite(r2, tmpsc2, site);
    nrm = norm2Diff(r1,r2);
    std::cout << "Results L-1\nref:  " << r1 << "\nimpl: " << r2 << "\ndiff (expect 0): " << nrm << std::endl;
    assert(nrm < 1e-10);

    site[1] = 1;
    peekSite(r1, tmpsc, site);
    peekSite(r2, tmpsc2, site);
    nrm = norm2Diff(r1,r2);
    std::cout << "Results 1\nref:  " << r1 << "\nimpl: " << r2 << "\ndiff (expect 0): " << nrm << std::endl;
    assert(nrm < 1e-10);

    //Test lower boundary
    site[1] = L[0]-1;
    tmpsc3 = Zero();
    pokeSite(v, tmpsc3, site);


    xconj_action_ref.M(tmpsc3, tmpsc);
    xconj_action.M(tmpsc3, tmpsc2);
    nrm = norm2Diff(tmpsc,tmpsc2);
    std::cout << "Test of Xconjugate action M lower boundary only (expect 0): " << nrm << std::endl;
    assert(nrm < 1e-10);

    site[1] = 0;

    peekSite(r1, tmpsc, site);
    peekSite(r2, tmpsc2, site);
    nrm = norm2Diff(r1,r2);
    std::cout << "Results 0\nref:  " << r1 << "\nimpl: " << r2 << "\ndiff (expect 0): " << nrm << std::endl;
    assert(nrm < 1e-10);

    site[1] = L[0]-2;

    peekSite(r1, tmpsc, site);
    peekSite(r2, tmpsc2, site);
    nrm = norm2Diff(r1,r2);
    std::cout << "Results L-2\nref:  " << r1 << "\nimpl: " << r2 << "\ndiff (expect 0): " << nrm << std::endl;
    assert(nrm < 1e-10);

    //Test full
    xconj_action_ref.M(rand_sc, tmpsc);
    xconj_action.M(rand_sc, tmpsc2);
    nrm = norm2Diff(tmpsc,tmpsc2);
    std::cout << "Test of Xconjugate action M against reference, full field (expect 0): " << nrm << std::endl;
    assert(nrm < 1e-10);

    xconj_action_ref.Mdag(rand_sc, tmpsc);
    xconj_action.Mdag(rand_sc, tmpsc2);
    nrm = norm2Diff(tmpsc,tmpsc2);
    std::cout << "Test of Xconjugate action Mdag against reference, full field (expect 0): " << nrm << std::endl;
    assert(nrm < 1e-10);

    xbarconj_action_ref.M(rand_sc, tmpsc);
    xbarconj_action.M(rand_sc, tmpsc2);
    nrm = norm2Diff(tmpsc,tmpsc2);
    std::cout << "Test of Xbar-conjugate action M against reference, full field (expect 0): " << nrm << std::endl;
    assert(nrm < 1e-10);

    xbarconj_action_ref.Mdag(rand_sc, tmpsc);
    xbarconj_action.Mdag(rand_sc, tmpsc2);
    nrm = norm2Diff(tmpsc,tmpsc2);
    std::cout << "Test of Xbar-conjugate action Mdag against reference, full field (expect 0): " << nrm << std::endl;
    assert(nrm < 1e-10);

    RealD u, l;

    //Test the X-conjugate Dirac op acting on a field is the same as the G-parity Dirac op acting on the equivalent 2f field
    xconj_action.M(rand_sc, tmpsc);
    boostXconjToGparity(tmp, rand_sc);
    reg_action.M(tmp, tmp2);
    norm2DiffXconj(u,l,tmp2,tmpsc);
    std::cout << "Test X-conj Dop reproduces G-parity Dop acting on X-conjugate field, f=0 (expect 0): " << u << " f=1 (expect 0): " << l << std::endl;
    assert(l < 1e-10);
    assert(u < 1e-10);

    //Test the Xbar-conjugate Dirac op acting on a field is the same as the G-parity Dirac op acting on the equivalent 2f field
    xbarconj_action.M(rand_sc, tmpsc);
    boostXbarConjToGparity(tmp, rand_sc);
    reg_action.M(tmp, tmp2);
    norm2DiffXbarConj(u,l,tmp2,tmpsc);
    std::cout << "Test Xbar-conj Dop reproduces G-parity Dop acting on Xbar-conjugate field, f=0 (expect 0): " << u << " f=1 (expect 0): " << l << std::endl;
    assert(l < 1e-10);
    assert(u < 1e-10);

    //Test the X-conj 2f wrapper reproduces G-parity Dop acting on X-conjugate field 
    Xconjugate2fWrapper<XconjugateMobiusFermionD, GparityMobiusFermionD> xconj_2f_wrapper(&xconj_action,&reg_action);
    boostXconjToGparity(tmp, rand_sc);
    xconj_2f_wrapper.M(tmp, tmp3);
    reg_action.M(tmp, tmp2);
    nrm = norm2Diff(tmp3,tmp2);
    std::cout << "Test the X-conj 2f wrapper reproduces G-parity Dop acting on X-conjugate field (expect 0): " << nrm << std::endl;
    assert(nrm < 1e-10);

    //Show that U^dag * v with X-conj v is real
    boostXconjToGparity(tmp, rand_sc);
    applyUdag(tmp2,tmp);
    tmp3 = tmp2 - conjugate(tmp2);
    nrm = norm2(tmp3);
    std::cout << "Test U^dag v with v an X-conj vector results in a real vector (expect 0): " << nrm << std::endl;
    assert(nrm < 1e-10);

    //Check its form
    tmpsc = PeekIndex<GparityFlavourIndex>(tmp2,0);
    tmpsc2 = sqrt(2.)*(Xmatrix()*real(rand_sc));
    nrm = norm2Diff(tmpsc,tmpsc2);
    std::cout << "Test U^dag v upper component has expected form (expect 0): " << nrm << std::endl;
    assert(nrm < 1e-10);

    tmpsc = PeekIndex<GparityFlavourIndex>(tmp2,1);
    tmpsc2 = sqrt(2.)*imag(rand_sc);
    nrm = norm2Diff(tmpsc,tmpsc2);
    std::cout << "Test U^dag v lower component has expected form (expect 0): " << nrm << std::endl;
    assert(nrm < 1e-10);

    //Show that U is unitary
    gaussian(RNG5,tmp);
    applyUdag(tmp2,tmp);
    applyU(tmp3,tmp2);
    nrm = norm2Diff(tmp3,tmp);
    std::cout << "Test U U^dag = 1 (expect 0): " << nrm << std::endl;
    assert(nrm < 1e-10);
    
    applyU(tmp2,tmp);
    applyUdag(tmp3,tmp2);
    nrm = norm2Diff(tmp3,tmp);
    std::cout << "Test U^dag U = 1 (expect 0): " << nrm << std::endl;
    assert(nrm < 1e-10);

    {
      //Show -U^* U^dag = -U U^T = Xi
      gaussian(RNG5,tmp);

      tmp2 = conjugate(tmp);   // U^T v = (U^dag v^*)^*
      applyUdag(tmp3,tmp2);
      tmp2 = conjugate(tmp3);
      applyU(tmp3,tmp2); //U U^T v
      tmp3 = -tmp3;

      //Xi = -i sigma2 X
      GparityFlavour sigma2 = GparityFlavour(GparityFlavour::Algebra::SigmaY);
      tmp2 = ComplexD(0,-1)*(sigma2*(Xmatrix()*tmp));
      nrm = norm2Diff(tmp3,tmp2);
      std::cout << "Test Xi = -U U^T (expect 0): " << nrm << std::endl;
      assert(nrm < 1e-10);
    }


    {
      //Show R = U^dag M U is a real matrix through
      // (R v)^* = R v^*
      FermionField2f Rv_allstar(FGrid), R_vstar(FGrid), v(FGrid);
      gaussian(RNG5,v);
      applyU(tmp, v);
      reg_action.M(tmp, tmp2);
      applyUdag(Rv_allstar,tmp2); //Rv
      Rv_allstar = conjugate(Rv_allstar);

      tmp = conjugate(v);
      applyU(tmp2,tmp);
      reg_action.M(tmp2, tmp);
      applyUdag(R_vstar,tmp); //Rv

      nrm = norm2Diff(R_vstar,Rv_allstar);
      std::cout << "Test U^dag M U is a real matrix (expect 0): " << nrm << std::endl;
      assert(nrm < 1e-10);
    }

    {
      //Test relations between elements of M
      //M11 = -X M00^* X
      FermionField1f tmp(FGrid),tmp2(FGrid),tmp3(FGrid), v(FGrid), lhs(FGrid), rhs(FGrid);
      gaussian(RNG5,v);
      tmp = Xmatrix()*v;
      applyMijstar(tmp2,tmp,0,0,reg_action);
      rhs = -(Xmatrix()*tmp2);

      applyMij(lhs,v,1,1,reg_action);
      nrm = norm2Diff(lhs, rhs);
      std::cout << "Test M11 = -X M00* X (expect 0): " << nrm << std::endl;
      assert(nrm < 1e-10);

      //M10 = X M01* X
      tmp = Xmatrix()*v;  //X M01^* X v = (X M01 X v^* )^*
      applyMijstar(tmp2,tmp,0,1,reg_action);
      rhs = Xmatrix()*tmp2;  

      applyMij(lhs,v,1,0,reg_action);
      nrm = norm2Diff(lhs, rhs);
      std::cout << "Test M10 = X M01* X (expect 0): " << nrm << std::endl;
      assert(nrm < 1e-10);
    }

    {
      //Test expressions for elements of R
      //R11 = 0.5(M11 + M12 X + M12* X + M11*)
      FermionField1f tmp(FGrid),tmp2(FGrid),tmp3(FGrid), v(FGrid), lhs(FGrid), rhs(FGrid);
      gaussian(RNG5,v);
      applyMij(rhs,v,0,0,reg_action);

      tmp = Xmatrix()*v;
      applyMij(tmp2,tmp,0,1,reg_action);
      rhs = rhs + tmp2;

      tmp = Xmatrix()*v;
      applyMijstar(tmp2,tmp,0,1,reg_action);
      rhs = rhs + tmp2;
            
      applyMijstar(tmp2,v,0,0,reg_action);
      rhs = rhs + tmp2;

      rhs = rhs*0.5;

      applyRij(lhs,v,1,1,reg_action);
      
      nrm = norm2Diff(lhs, rhs);
      std::cout << "Test R11 = 0.5(M11 + M12 X + M12* X + M11*) (expect 0): " << nrm << std::endl;
      assert(nrm < 1e-10);

      //R00 = 0.5( -X M00 X   -X M01   -X M01*  -X M00* X )
      tmp = Xmatrix()*v;
      applyMij(tmp2,tmp,0,0,reg_action);
      rhs = -(Xmatrix()*tmp2);
      
      applyMij(tmp,v,0,1,reg_action);
      tmp2 = -(Xmatrix()*tmp);
      rhs = rhs + tmp2;

      applyMijstar(tmp,v,0,1,reg_action);
      tmp2 = -(Xmatrix()*tmp);
      rhs = rhs + tmp2;

      tmp = Xmatrix()*v;
      applyMijstar(tmp2,tmp,0,0,reg_action);
      tmp = -(Xmatrix()*tmp2);
      rhs = rhs + tmp;

      rhs = rhs * 0.5;
      
      applyRij(lhs,v,0,0,reg_action);

      nrm = norm2Diff(lhs, rhs);
      std::cout << "Test R00 = 0.5( -X M00 X   -X M01   -X M01*  -X M00* X ) (expect 0): " << nrm << std::endl;
      assert(nrm < 1e-10);

      //R01 = 0.5*( iX M00 + iX M01 X  -iX M01*X - iX M00* )
      applyMij(tmp,v,0,0,reg_action);
      rhs = Xmatrix()*tmp;
      
      tmp = Xmatrix()*v;
      applyMij(tmp2,tmp,0,1,reg_action);
      tmp = Xmatrix()*tmp2;
      rhs = rhs + tmp;

      tmp = Xmatrix()*v;
      applyMijstar(tmp2,tmp,0,1,reg_action);
      tmp = Xmatrix()*tmp2;
      rhs = rhs - tmp;

      applyMijstar(tmp,v,0,0,reg_action);
      tmp2 = Xmatrix()*tmp;
      rhs = rhs - tmp2;

      rhs = rhs * ComplexD(0,0.5);

      applyRij(lhs,v,0,1,reg_action);

      nrm = norm2Diff(lhs, rhs);
      std::cout << "Test R01 = 0.5*( iX M00 + iX M01 X  -iX M01*X - iX M00* ) (expect 0): " << nrm << std::endl;
      assert(nrm < 1e-10);

      //R10 = 0.5( i M00 X + i M01 -i M01* - iM00* X )
      tmp = Xmatrix() * v;
      applyMij(rhs,tmp,0,0,reg_action);
      
      applyMij(tmp,v,0,1,reg_action);
      rhs = rhs + tmp;

      applyMijstar(tmp,v,0,1,reg_action);
      rhs = rhs - tmp;
      
      tmp = Xmatrix() * v;
      applyMijstar(tmp2,tmp,0,0,reg_action);
      rhs = rhs - tmp2;
      
      rhs = rhs * ComplexD(0,0.5);

      applyRij(lhs,v,1,0,reg_action);

      nrm = norm2Diff(lhs, rhs);     
      std::cout << "Test R10 = 0.5( i M00 X + i M01 -i M01* - iM00* X ) (expect 0): " << nrm << std::endl;
      assert(nrm < 1e-10);
    }
    
    {
      //Test reference implementation
      RotatedGPrefAction<GparityMobiusFermionD> real_action_ref(&reg_action);
      FermionField1f tmp(FGrid),tmp2(FGrid),tmp3(FGrid), v(FGrid), lhs(FGrid), rhs(FGrid);
      gaussian(RNG5,v);
      applyRijRef(lhs,v,0,0,real_action_ref);
      applyRij(rhs,v,0,0,reg_action);
      nrm = norm2Diff(lhs, rhs);
      std::cout << "Test reference impl R00 (expect 0): " << nrm << std::endl;
      //assert(nrm < 1e-10);

      applyRijRef(lhs,v,0,1,real_action_ref);
      applyRij(rhs,v,0,1,reg_action);
      nrm = norm2Diff(lhs, rhs);
      std::cout << "Test reference impl R01 (expect 0): " << nrm << std::endl;
      //assert(nrm < 1e-10);

      applyRijRef(lhs,v,1,0,real_action_ref);
      applyRij(rhs,v,1,0,reg_action);
      nrm = norm2Diff(lhs, rhs);
      std::cout << "Test reference impl R10 (expect 0): " << nrm << std::endl;
      //assert(nrm < 1e-10);

      applyRijRef(lhs,v,1,1,real_action_ref);
      applyRij(rhs,v,1,1,reg_action);
      nrm = norm2Diff(lhs, rhs);
      std::cout << "Test reference impl R11 (expect 0): " << nrm << std::endl;
      //assert(nrm < 1e-10);
    }

    {
      //Test reference implementation of U^dag M U
      FermionField2f Rv_test(FGrid), Rv_ref(FGrid), v(FGrid);
      gaussian(RNG5,v);
      applyU(tmp, v);
      reg_action.M(tmp, tmp2);
      applyUdag(Rv_test,tmp2); //Rv
 
      RotatedGPrefAction<GparityMobiusFermionD> real_action_ref(&reg_action);
      real_action_ref.M(v,Rv_ref);
      
      nrm = norm2Diff(Rv_test,Rv_ref);
      std::cout << "Test reference implementation of real matrix M (expect 0): " << nrm << std::endl;
      assert(nrm < 1e-10);
    }
    exit(0);
    

  }

  //########################################################################
  {
    std::cout << "Test of red-black preconditioned operator" << std::endl;
    SchurDiagMooeeOperator<GparityMobiusFermionD,FermionField2f> reg_schurop(reg_action);
    SchurDiagMooeeOperator<XconjugateMobiusFermionD,FermionField1f> xconj_schurop(xconj_action);
    SchurDiagMooeeOperator<XconjugateMobiusFermionD,FermionField1f> xbarconj_schurop(xbarconj_action);
    XconjugateDWF<GparityMobiusFermionD> xconj_action_ref(&reg_action);

    FermionField1f rand_sc(FrbGrid), rand_sc2(FrbGrid); //v
    gaussian(RNG5rb,rand_sc);
    gaussian(RNG5rb,rand_sc2);

    FermionField2f rand_f0(FrbGrid); //v \delta_f,0
    rand_f0 = Zero();    
    PokeIndex<GparityFlavourIndex>(rand_f0, rand_sc, 0);

    FermionField2f rand_f1(FrbGrid); //v \delta_f,1
    rand_f1 = Zero();
    PokeIndex<GparityFlavourIndex>(rand_f1, rand_sc, 1);
    
    FermionField2f tmp(FrbGrid), tmp2(FrbGrid), tmp3(FrbGrid), tmp4(FrbGrid);
    FermionField1f tmpsc(FrbGrid), tmpsc2(FrbGrid), tmpsc3(FrbGrid);
    RealD nrm;

    std::cout << "Test the relationship between the upper and lower rows in flavor space of the G-parity preconditioned Dirac operator" << std::endl;    
    
    reg_schurop.Mpc(rand_f0,tmp);
    FermionField1f M00v = PeekIndex<GparityFlavourIndex>(tmp,0);
    FermionField1f M10v = PeekIndex<GparityFlavourIndex>(tmp,1);

    //Test  M00 v =  -X M11* X v = -X [M11 X v*]*
    tmp = X*conjugate(rand_f1);    

    reg_schurop.Mpc(tmp,tmp2);
    FermionField1f M11Xvconj = PeekIndex<GparityFlavourIndex>(tmp2,1);
    FermionField1f M01Xvconj = PeekIndex<GparityFlavourIndex>(tmp2,0);

    tmpsc = -(X*conjugate(M11Xvconj));
    nrm = norm2Diff(tmpsc, M00v);
    std::cout << "Test of M00 v =  -X M11* X v = -X [M11 X v*]*  (expect 0): " << nrm << std::endl;
    assert(nrm < 1e-10);

    //Test  M10 v =  X M01* X v = X [M01 X v*]*
    tmpsc = X*conjugate(M01Xvconj);
    nrm = norm2Diff(tmpsc,M10v);
    std::cout << "Test of M10 v =  X M11* X v = X [M11 X v*]*  (expect 0): " << nrm << std::endl;
    assert(nrm < 1e-10);

    std::cout << "Test the relationship between the upper and lower rows in flavor space of the G-parity preconditioned squared Dirac operator" << std::endl;    

    reg_schurop.HermOp(rand_f0,tmp);
    M00v = PeekIndex<GparityFlavourIndex>(tmp,0);
    M10v = PeekIndex<GparityFlavourIndex>(tmp,1);

    //Test  M00 v =  -X M11* X v = -X [M11 X v*]*
    tmp = X*conjugate(rand_f1);    

    reg_schurop.HermOp(tmp,tmp2);
    M11Xvconj = PeekIndex<GparityFlavourIndex>(tmp2,1);
    M01Xvconj = PeekIndex<GparityFlavourIndex>(tmp2,0);

    tmpsc = -(X*conjugate(M11Xvconj));
    nrm = norm2Diff(tmpsc,M00v);
    std::cout << "Test of M00 v =  -X M11* X v = -X [M11 X v*]*  (expect 0): " << nrm << std::endl;
    assert(nrm < 1e-10);

    //Test  M10 v =  X M01* X v = X [M01 X v*]*
    tmpsc = X*conjugate(M01Xvconj);
    nrm = norm2Diff(tmpsc,M10v);
    std::cout << "Test of M10 v =  X M11* X v = X [M11 X v*]*  (expect 0): " << nrm << std::endl;
    assert(nrm < 1e-10);

    //Test the X-conjugate implementation
    std::cout << "Test the implementation of the X-conjugate preconditioned action against the reference"<< std::endl;
    xconj_action_ref.Meooe(rand_sc, tmpsc);
    xconj_action.Meooe(rand_sc, tmpsc2);
    nrm = norm2Diff(tmpsc,tmpsc2);
    std::cout << "Test of X-conjugate action Meooe against reference, full field (expect 0): " << nrm << std::endl;
    assert(nrm < 1e-10);

    xconj_action_ref.MeooeDag(rand_sc, tmpsc);
    xconj_action.MeooeDag(rand_sc, tmpsc2);
    nrm = norm2Diff(tmpsc,tmpsc2);
    std::cout << "Test of X-conjugate action MeooeDag against reference, full field (expect 0): " << nrm << std::endl;
    assert(nrm < 1e-10);

    xconj_action_ref.Mooee(rand_sc, tmpsc);
    xconj_action.Mooee(rand_sc, tmpsc2);
    nrm = norm2Diff(tmpsc,tmpsc2);
    std::cout << "Test of X-conjugate action Mooee against reference, full field (expect 0): " << nrm << std::endl;
    assert(nrm < 1e-10);

    xconj_action_ref.MooeeDag(rand_sc, tmpsc);
    xconj_action.MooeeDag(rand_sc, tmpsc2);
    nrm = norm2Diff(tmpsc,tmpsc2);
    std::cout << "Test of X-conjugate action MooeeDag against reference, full field (expect 0): " << nrm << std::endl;
    assert(nrm < 1e-10);
        
    xconj_action_ref.MooeeInv(rand_sc, tmpsc);
    xconj_action.MooeeInv(rand_sc, tmpsc2);
    nrm = norm2Diff(tmpsc,tmpsc2);
    std::cout << "Test of X-conjugate action MooeeInv against reference, full field (expect 0): " << nrm << std::endl;
    assert(nrm < 1e-10);

    xconj_action_ref.MooeeInvDag(rand_sc, tmpsc);
    xconj_action.MooeeInvDag(rand_sc, tmpsc2);
    nrm = norm2Diff(tmpsc,tmpsc2);
    std::cout << "Test of X-conjugate action MooeeInvDag against reference, full field (expect 0): " << nrm << std::endl;
    assert(nrm < 1e-10);

    RealD u,l;

    //Test the X-conjugate Dirac op acting on a field is the same as the G-parity Dirac op acting on a field with explicit X-conjugate BCs
    xconj_schurop.HermOp(rand_sc, tmpsc);
    boostXconjToGparity(tmp, rand_sc);
    reg_schurop.HermOp(tmp, tmp2);
    norm2DiffXconj(u,l,tmp2,tmpsc);
    std::cout << "Test X-conj HermOp reproduces G-parity HermOp acting on X-conjugate field, f=0 (expect 0): " << u << " f=1 (expect 0): " << l << std::endl;
    assert(u < 1e-10);
    assert(l < 1e-10);
    
    //Test the X-conj 2f wrapper reproduces G-parity Dop acting on X-conjugate field
    Xconjugate2fWrapper<XconjugateMobiusFermionD,GparityMobiusFermionD> xconj_2f_wrapper(&xconj_action,&reg_action);
    SchurDiagMooeeOperator<Xconjugate2fWrapper<XconjugateMobiusFermionD,GparityMobiusFermionD>, FermionField2f> xconj_2f_schurop(xconj_2f_wrapper);
    boostXconjToGparity(tmp, rand_sc);
    reg_schurop.HermOp(tmp, tmp2);
    xconj_2f_schurop.HermOp(tmp,tmp3);
    nrm = norm2Diff(tmp3,tmp2);
    std::cout << "Test the X-conj 2f wrapper reproduces G-parity HermOp acting on X-conjugate field (expect 0): " << nrm << std::endl;
    assert(nrm < 1e-10);

    //Test the Xbar-conjugate Dirac op acting on a field is the same as the G-parity Dirac op acting on a field with explicit Xbar-conjugate BCs
    xbarconj_schurop.HermOp(rand_sc, tmpsc);
    boostXbarConjToGparity(tmp,rand_sc);
    reg_schurop.HermOp(tmp, tmp2);
    norm2DiffXbarConj(u,l,tmp2,tmpsc);
    std::cout << "Test Xbar-conj HermOp reproduces G-parity HermOp acting on Xbar-conjugate field, f=0 (expect 0): " << u << " f=1 (expect 0): " << l << std::endl;
    assert(u < 1e-10);
    assert(l < 1e-10);   

    //Test the Xbar-conj 2f wrapper reproduces G-parity Dop acting on Xbar-conjugate field
    Xconjugate2fWrapper<XconjugateMobiusFermionD,GparityMobiusFermionD> xbarconj_2f_wrapper(&xbarconj_action,&reg_action);
    SchurDiagMooeeOperator<Xconjugate2fWrapper<XconjugateMobiusFermionD,GparityMobiusFermionD>, FermionField2f> xbarconj_2f_schurop(xbarconj_2f_wrapper);
    boostXbarConjToGparity(tmp, rand_sc);
    reg_schurop.HermOp(tmp, tmp2);
    xbarconj_2f_schurop.HermOp(tmp,tmp3);
    nrm = norm2Diff(tmp3,tmp2);
    std::cout << "Test the Xbar-conj 2f wrapper reproduces G-parity HermOp acting on Xbar-conjugate field (expect 0): " << nrm << std::endl;
    assert(nrm < 1e-10);

    //Test reconstruction of G-parity Dop on arbitrary flavor vector using Xconj action
    PokeIndex<GparityFlavourIndex>(tmp, rand_sc, 0);
    PokeIndex<GparityFlavourIndex>(tmp, rand_sc2, 1);
    reg_schurop.HermOp(tmp, tmp2);
    
    FermionField1f rho(FrbGrid), xi(FrbGrid);
    rho = 0.5 * ( rand_sc + (X*conjugate(rand_sc2)) );
    xi = 0.5 * ( rand_sc - (X*conjugate(rand_sc2)) );
    xconj_schurop.HermOp(rho, tmpsc);
    xbarconj_schurop.HermOp(xi, tmpsc2);

    tmpsc3 = PeekIndex<GparityFlavourIndex>(tmp2,0) - tmpsc - tmpsc2;
    u = norm2(tmpsc3);

    tmpsc = -(X*conjugate(tmpsc));
    tmpsc2 = X*conjugate(tmpsc2);
    tmpsc3 = PeekIndex<GparityFlavourIndex>(tmp2,1) - tmpsc - tmpsc2;
    l = norm2(tmpsc3);

    std::cout << "Test reconstruction of GP HermOp on random 2f field from Xconj ops f=0 (expect 0): " << u << " f=1 (expect 0): " << l << std::endl;
    assert(u < 1e-10);
    assert(l < 1e-10);
  }
  std::cout << "All tests passed" << std::endl;
  
  Grid_finalize();
}

