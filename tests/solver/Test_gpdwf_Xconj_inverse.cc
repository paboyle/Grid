    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/solver/Test_gpdwf_Xconj_inverse.cc

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

typedef typename GparityMobiusFermionD::FermionField FermionField2f;
typedef typename XconjugateMobiusFermionD::FermionField FermionField1f;

void invertGparity(std::vector<FermionField2f> &out, GparityMobiusFermionD &action){
  GridBase* UGrid = action.GaugeGrid();
  GridBase* FGrid = action.FermionGrid();
  GridBase* FrbGrid = action.FermionRedBlackGrid();

  ConjugateGradient<FermionField2f> cg(1e-08,10000);
  SchurRedBlackDiagMooeeSolve<FermionField2f> solver(cg);

  typedef FermionField2f::scalar_object Spinor;

  LatticeInteger tcoor4d(UGrid);
  LatticeCoordinate(tcoor4d,3);

  FermionField2f tmp4d(UGrid);
  FermionField2f src4d(UGrid);
  FermionField2f src5d(FGrid), sol5d(FGrid);
  FermionField2f src5d_e(FrbGrid), src5d_o(FrbGrid), sol5d_o(FrbGrid);

  out = std::vector<FermionField2f>(24, UGrid);

  for(int f=0;f<2;f++){
    for(int s=0;s<4;s++){
      for(int c=0;c<3;c++){
	Spinor v = Zero();
	v(f)(s)(c) = 1.;

	tmp4d = v;
	src4d = Zero();
	src4d = where( tcoor4d == Integer(0), tmp4d, src4d); //unit vector on every site on timeslice 0, zero elsewhere

	action.ImportPhysicalFermionSource(src4d, src5d);
	solver.RedBlackSource(action, src5d, src5d_e, src5d_o);
	solver.RedBlackSolve(action, src5d_o, sol5d_o);
	solver.RedBlackSolution(action, sol5d_o, src5d_e, sol5d);
	action.ExportPhysicalFermionSolution(sol5d, out[c+3*(s + 4*f)]);
      }
    }
  }
}

//0.5*(1+iX) in
template<typename T>
T mulPplusLeft(const T &in){
  static Gamma C = Gamma(Gamma::Algebra::MinusGammaY) * Gamma(Gamma::Algebra::GammaT);
  static Gamma g5 = Gamma(Gamma::Algebra::Gamma5);
  static Gamma X = C*g5;
  static ComplexD _I(0,1);
  T out = 0.5*(in + _I * (X*in));
  return out;
}
//0.5*(1+iX) in
template<typename T>
T mulPminusLeft(const T &in){
  static Gamma C = Gamma(Gamma::Algebra::MinusGammaY) * Gamma(Gamma::Algebra::GammaT);
  static Gamma g5 = Gamma(Gamma::Algebra::Gamma5);
  static Gamma X = C*g5;
  static ComplexD _mI(0,-1);
  T out = 0.5*(in + _mI * (X*in));
  return out;
}

//0.5*(1+iX) in
template<typename T>
T mulPplusRight(const T &in){
  static Gamma C = Gamma(Gamma::Algebra::MinusGammaY) * Gamma(Gamma::Algebra::GammaT);
  static Gamma g5 = Gamma(Gamma::Algebra::Gamma5);
  static Gamma X = C*g5;
  static ComplexD _I(0,1);
  T out = 0.5*(in + _I * (in*X));
  return out;
}
//0.5*(1+iX) in
template<typename T>
T mulPminusRight(const T &in){
  static Gamma C = Gamma(Gamma::Algebra::MinusGammaY) * Gamma(Gamma::Algebra::GammaT);
  static Gamma g5 = Gamma(Gamma::Algebra::Gamma5);
  static Gamma X = C*g5;
  static ComplexD _mI(0,-1);
  T out = 0.5*(in + _mI * (in*X));
  return out;
}



typedef iScalar<iMatrix<iMatrix<typename FermionField1f::vector_type, Nc>, Ns> > vSCmatrix;
typedef iScalar<iVector<iVector<typename FermionField1f::vector_type, Nc>, Ns> > vSCvector;
typedef Lattice<vSCmatrix> SCmatrixField;

//Poke the spin-color vector field 'from' onto the sc,cc spin/color column of 'into'
void pokeSpinColorColumn(SCmatrixField &into, const FermionField1f &from, const int sc, const int cc){
  size_t Nsimd = FermionField1f::vector_type::Nsimd();
  autoView( from_v, from, AcceleratorRead);
  autoView( into_v, into, AcceleratorWrite);
  accelerator_for( i, from_v.size(), Nsimd, {
      auto site_from = from_v(i);
      auto site_into = into_v(i);
      for(int sr=0;sr<Ns;sr++)
	for(int cr=0;cr<Nc;cr++)
	  site_into()(sr,sc)(cr,cc) = site_from()(sr)(cr);
      coalescedWrite(into_v[i], site_into);
    });
}
void peekSpinColorColumn(FermionField1f &into, const SCmatrixField&from, const int sc, const int cc){
  size_t Nsimd = FermionField1f::vector_type::Nsimd();
  autoView( from_v, from, AcceleratorRead);
  autoView( into_v, into, AcceleratorWrite);
  accelerator_for( i, from_v.size(), Nsimd, {
      auto site_from = from_v(i);
      auto site_into = into_v(i);
      for(int sr=0;sr<Ns;sr++)
	for(int cr=0;cr<Nc;cr++)
	  site_into()(sr)(cr) = site_from()(sr,sc)(cr,cc);
      coalescedWrite(into_v[i], site_into);
    });
}


typedef iMatrix<iMatrix<iMatrix<typename FermionField1f::vector_type, Nc>, Ns>, Ngp> vSCFmatrix;
typedef iVector<iVector<iVector<typename FermionField1f::vector_type, Nc>, Ns>, Ngp> vSCFvector;
typedef Lattice<vSCFmatrix> SCFmatrixField;

//Poke the spin-color-flavor vector field 'from' onto the sc,cc,ff spin/color/flavor column of 'into'
void pokeSpinColorFlavorColumn(SCFmatrixField &into, const FermionField2f &from, const int sc, const int cc, const int fc){
  size_t Nsimd = FermionField2f::vector_type::Nsimd();
  autoView( from_v, from, AcceleratorRead);
  autoView( into_v, into, AcceleratorWrite);
  accelerator_for( i, from_v.size(), Nsimd, {
      auto site_from = from_v(i);
      auto site_into = into_v(i);
      for(int fr=0;fr<Ngp;fr++)
	for(int sr=0;sr<Ns;sr++)
	  for(int cr=0;cr<Nc;cr++)
	    site_into(fr,fc)(sr,sc)(cr,cc) = site_from(fr)(sr)(cr);
      coalescedWrite(into_v[i], site_into);
    });
}
void peekSpinColorFlavorColumn(FermionField2f &into, const SCFmatrixField&from, const int sc, const int cc, const int fc){
  size_t Nsimd = FermionField2f::vector_type::Nsimd();
  autoView( from_v, from, AcceleratorRead);
  autoView( into_v, into, AcceleratorWrite);
  accelerator_for( i, from_v.size(), Nsimd, {
      auto site_from = from_v(i);
      auto site_into = into_v(i);
      for(int fr=0;fr<Ngp;fr++)	
	for(int sr=0;sr<Ns;sr++)
	  for(int cr=0;cr<Nc;cr++)
	    site_into(fr)(sr)(cr) = site_from(fr,fc)(sr,sc)(cr,cc);
      coalescedWrite(into_v[i], site_into);
    });
}




void invertXconj2d(std::vector<FermionField2f> &out, GparityMobiusFermionD &action){
  GridBase* UGrid = action.GaugeGrid();
  GridBase* FGrid = action.FermionGrid();
  GridBase* FrbGrid = action.FermionRedBlackGrid();

  ConjugateGradient<FermionField2f> cg(1e-08,10000);
  SchurRedBlackDiagMooeeSolve<FermionField2f> solver(cg);

  typedef FermionField2f::scalar_object Spinor;
  typedef iScalar<iVector<iVector<Spinor::vector_type, Nc>, Ns> > FSpinor;

  LatticeInteger tcoor4d(UGrid);
  LatticeCoordinate(tcoor4d,3);

  FermionField2f tmp4d(UGrid);
  FermionField2f src4d(UGrid);
  FermionField2f src5d(FGrid), sol5d(FGrid);
  FermionField2f src5d_e(FrbGrid), src5d_o(FrbGrid), sol5d_o(FrbGrid);

  SCmatrixField tmp_vplus_f0_4d(UGrid);
  SCmatrixField tmp_vminus_f0_4d(UGrid);

  FermionField1f tmp1f(UGrid);

  static Gamma C = Gamma(Gamma::Algebra::MinusGammaY) * Gamma(Gamma::Algebra::GammaT);
  static Gamma g5 = Gamma(Gamma::Algebra::Gamma5);
  static Gamma X = C*g5;

  for(int s=0;s<4;s++){
    for(int c=0;c<3;c++){
      FSpinor vb = Zero();
      vb()(s)(c) = 1.;
      
      Spinor vplus;
      vplus(0) = ( mulPplusLeft(vb) )();
      vplus(1) = ( -(X*mulPminusLeft(vb)) )();

      Spinor vminus;
      vminus(0) = ( mulPminusLeft(vb) )();
      vminus(1) = ( -(X*mulPplusLeft(vb)) )();

      Spinor* vpm[2] = {&vplus, &vminus};
      SCmatrixField* solpm[2] = {&tmp_vplus_f0_4d, &tmp_vminus_f0_4d};

      for(int pm=0;pm<2;pm++){
	tmp4d = *vpm[pm];
	src4d = Zero();
	src4d = where( tcoor4d == Integer(0), tmp4d, src4d); //vector on every site on timeslice 0, zero elsewhere

	action.ImportPhysicalFermionSource(src4d, src5d);
	solver.RedBlackSource(action, src5d, src5d_e, src5d_o);
	solver.RedBlackSolve(action, src5d_o, sol5d_o);
	solver.RedBlackSolution(action, sol5d_o, src5d_e, sol5d);
	action.ExportPhysicalFermionSolution(sol5d, tmp4d);
	tmp1f = PeekIndex<GparityFlavourIndex>(tmp4d,0); //only need f0 component
	pokeSpinColorColumn(*solpm[pm], tmp1f, s,c);
      }
    }
  }

  SCmatrixField vplus_f0_Pplus = mulPplusRight(tmp_vplus_f0_4d);
  SCmatrixField vplus_f0_Pminus = mulPminusRight(tmp_vplus_f0_4d);

  SCmatrixField vminus_f0_Pplus = mulPplusRight(tmp_vminus_f0_4d);
  SCmatrixField vminus_f0_Pminus = mulPminusRight(tmp_vminus_f0_4d);

  SCmatrixField Minv11 = vplus_f0_Pplus + vminus_f0_Pminus;
  SCmatrixField Minv12 = vplus_f0_Pminus*X + vminus_f0_Pplus*X;
  SCmatrixField Minv22 = -(X*(conjugate(Minv11)*X));
  SCmatrixField Minv21 = X*(conjugate(Minv12)*X);
  
  out = std::vector<FermionField2f>(24, UGrid);

  SCmatrixField* Minv[2][2] = { {&Minv11,&Minv12}, {&Minv21,&Minv22} };

  for(int fc=0;fc<2;fc++){
    for(int sc=0;sc<4;sc++){
      for(int cc=0;cc<3;cc++){
	int out_idx = cc+3*(sc + 4*fc);

	for(int fr=0;fr<2;fr++){
	  peekSpinColorColumn(tmp1f, *Minv[fr][fc]  , sc, cc);

	  PokeIndex<GparityFlavourIndex>(out[out_idx], tmp1f, fr);
	}
      }
    }
  }
}




void invertXconj1d(std::vector<FermionField2f> &out, XconjugateMobiusFermionD &action){
  GridBase* UGrid = action.GaugeGrid();
  GridBase* FGrid = action.FermionGrid();
  GridBase* FrbGrid = action.FermionRedBlackGrid();

  ConjugateGradient<FermionField1f> cg(1e-08,10000);
  SchurRedBlackDiagMooeeSolve<FermionField1f> solver(cg);

  typedef FermionField1f::scalar_object Spinor;

  LatticeInteger tcoor4d(UGrid);
  LatticeCoordinate(tcoor4d,3);

  FermionField1f tmp4d(UGrid);
  FermionField1f src4d(UGrid);
  FermionField1f src5d(FGrid), sol5d(FGrid);
  FermionField1f src5d_e(FrbGrid), src5d_o(FrbGrid), sol5d_o(FrbGrid);

  SCmatrixField tmp_vplus_f0_4d(UGrid);
  SCmatrixField tmp_vminus_f0_4d(UGrid);

  static Gamma C = Gamma(Gamma::Algebra::MinusGammaY) * Gamma(Gamma::Algebra::GammaT);
  static Gamma g5 = Gamma(Gamma::Algebra::Gamma5);
  static Gamma X = C*g5;

  for(int s=0;s<4;s++){
    for(int c=0;c<3;c++){
      Spinor vb = Zero();
      vb()(s)(c) = 1.;
      
      Spinor vplus = mulPplusLeft(vb);
      Spinor vminus = mulPminusLeft(vb);

      Spinor* vpm[2] = {&vplus, &vminus};
      SCmatrixField* solpm[2] = {&tmp_vplus_f0_4d, &tmp_vminus_f0_4d};

      for(int pm=0;pm<2;pm++){
	tmp4d = *vpm[pm];
	src4d = Zero();
	src4d = where( tcoor4d == Integer(0), tmp4d, src4d); //vector on every site on timeslice 0, zero elsewhere

	action.ImportPhysicalFermionSource(src4d, src5d);
	solver.RedBlackSource(action, src5d, src5d_e, src5d_o);
	solver.RedBlackSolve(action, src5d_o, sol5d_o);
	solver.RedBlackSolution(action, sol5d_o, src5d_e, sol5d);
	action.ExportPhysicalFermionSolution(sol5d, tmp4d);
	pokeSpinColorColumn(*solpm[pm], tmp4d, s,c);
      }
    }
  }

  SCmatrixField vplus_f0_Pplus = mulPplusRight(tmp_vplus_f0_4d);
  SCmatrixField vplus_f0_Pminus = mulPminusRight(tmp_vplus_f0_4d);

  SCmatrixField vminus_f0_Pplus = mulPplusRight(tmp_vminus_f0_4d);
  SCmatrixField vminus_f0_Pminus = mulPminusRight(tmp_vminus_f0_4d);

  SCmatrixField Minv11 = vplus_f0_Pplus + vminus_f0_Pminus;
  SCmatrixField Minv12 = vplus_f0_Pminus*X + vminus_f0_Pplus*X;
  SCmatrixField Minv22 = -(X*(conjugate(Minv11)*X));
  SCmatrixField Minv21 = X*(conjugate(Minv12)*X);
  
  out = std::vector<FermionField2f>(24, UGrid);

  SCmatrixField* Minv[2][2] = { {&Minv11,&Minv12}, {&Minv21,&Minv22} };

  for(int fc=0;fc<2;fc++){
    for(int sc=0;sc<4;sc++){
      for(int cc=0;cc<3;cc++){
	int out_idx = cc+3*(sc + 4*fc);

	for(int fr=0;fr<2;fr++){
	  peekSpinColorColumn(tmp4d, *Minv[fr][fc]  , sc, cc);

	  PokeIndex<GparityFlavourIndex>(out[out_idx], tmp4d, fr);
	}
      }
    }
  }
}


template<typename T>
T mulURight(const T &in){
  static Gamma C = Gamma(Gamma::Algebra::MinusGammaY) * Gamma(Gamma::Algebra::GammaT);
  static Gamma g5 = Gamma(Gamma::Algebra::Gamma5);
  static Gamma X = C*g5;
  static GparityFlavour sigma3 = GparityFlavour(GparityFlavour::Algebra::SigmaZ);

  T out = 0.5*(in + in*X) + 0.5*(in*sigma3 - (in*X)*sigma3);
  return out;
}



void invertXconj1d_matrix(std::vector<FermionField2f> &out, XconjugateMobiusFermionD &action){
  GridBase* UGrid = action.GaugeGrid();
  GridBase* FGrid = action.FermionGrid();
  GridBase* FrbGrid = action.FermionRedBlackGrid();

  ConjugateGradient<FermionField1f> cg(1e-08,10000);
  SchurRedBlackDiagMooeeSolve<FermionField1f> solver(cg);

  typedef FermionField1f::scalar_object Spinor;

  LatticeInteger tcoor4d(UGrid);
  LatticeCoordinate(tcoor4d,3);

  FermionField1f tmp4d(UGrid);
  FermionField1f src4d(UGrid);
  FermionField1f src5d(FGrid), sol5d(FGrid);
  FermionField1f src5d_e(FrbGrid), src5d_o(FrbGrid), sol5d_o(FrbGrid);

  FermionField2f tmp2f(UGrid);

  SCFmatrixField V_4d(UGrid);

  static Gamma C = Gamma(Gamma::Algebra::MinusGammaY) * Gamma(Gamma::Algebra::GammaT);
  static Gamma g5 = Gamma(Gamma::Algebra::Gamma5);
  static Gamma X = C*g5;
  static GparityFlavour sigma1 = GparityFlavour(GparityFlavour::Algebra::SigmaX);
  

  for(int s=0;s<4;s++){
    for(int c=0;c<3;c++){
      Spinor vb = Zero();
      vb()(s)(c) = 1.;
      
      Spinor vplus = mulPplusLeft(vb);
      Spinor vminus = mulPminusLeft(vb);

      Spinor* vpm[2] = {&vplus, &vminus};

      for(int pm=0;pm<2;pm++){
	tmp4d = *vpm[pm];
	src4d = Zero();
	src4d = where( tcoor4d == Integer(0), tmp4d, src4d); //vector on every site on timeslice 0, zero elsewhere

	action.ImportPhysicalFermionSource(src4d, src5d);
	solver.RedBlackSource(action, src5d, src5d_e, src5d_o);
	solver.RedBlackSolve(action, src5d_o, sol5d_o);
	solver.RedBlackSolution(action, sol5d_o, src5d_e, sol5d);
	action.ExportPhysicalFermionSolution(sol5d, tmp4d);

	//Generate 2f X-conjugate output
	PokeIndex<GparityFlavourIndex>(tmp2f, tmp4d, 0);
	tmp4d = -(X*conjugate(tmp4d));
	PokeIndex<GparityFlavourIndex>(tmp2f, tmp4d, 1);

	pokeSpinColorFlavorColumn(V_4d, tmp2f, s,c,pm);
      }
    }
  }

  SCFmatrixField tmp = mulPplusRight(V_4d) + mulPminusRight(V_4d)*sigma1;  
  SCFmatrixField Minv = mulURight(tmp);

  out = std::vector<FermionField2f>(24, UGrid);

  for(int fc=0;fc<2;fc++){
    for(int sc=0;sc<4;sc++){
      for(int cc=0;cc<3;cc++){
	int out_idx = cc+3*(sc + 4*fc);
	peekSpinColorFlavorColumn(out[out_idx], Minv, sc, cc, fc);
      }
    }
  }
}


//Use the full 2x2 G-parity Dirac op inverted on a complex source with flavor structure
// |\phi   0     |
// | 0    \phi^* |
void invertGparityComplexSrc(std::vector<FermionField2f> &out, GparityMobiusFermionD &action, ComplexD phi){
  GridBase* UGrid = action.GaugeGrid();
  GridBase* FGrid = action.FermionGrid();
  GridBase* FrbGrid = action.FermionRedBlackGrid();

  ConjugateGradient<FermionField2f> cg(1e-08,10000);
  SchurRedBlackDiagMooeeSolve<FermionField2f> solver(cg);

  typedef FermionField2f::scalar_object Spinor;

  LatticeInteger tcoor4d(UGrid);
  LatticeCoordinate(tcoor4d,3);

  FermionField2f tmp4d(UGrid);
  FermionField2f src4d(UGrid);
  FermionField2f src5d(FGrid), sol5d(FGrid);
  FermionField2f src5d_e(FrbGrid), src5d_o(FrbGrid), sol5d_o(FrbGrid);

  out = std::vector<FermionField2f>(24, UGrid);

  for(int f=0;f<2;f++){
    for(int s=0;s<4;s++){
      for(int c=0;c<3;c++){
	Spinor v = Zero();
	v(f)(s)(c) = f==0 ? phi : conjugate(phi);

	tmp4d = v;
	src4d = Zero();
	src4d = where( tcoor4d == Integer(0), tmp4d, src4d); //unit vector on every site on timeslice 0, zero elsewhere

	action.ImportPhysicalFermionSource(src4d, src5d);
	solver.RedBlackSource(action, src5d, src5d_e, src5d_o);
	solver.RedBlackSolve(action, src5d_o, sol5d_o);
	solver.RedBlackSolution(action, sol5d_o, src5d_e, sol5d);
	action.ExportPhysicalFermionSolution(sol5d, out[c+3*(s + 4*f)]);
      }
    }
  }
}

//Do the same as the above but with the X-conjugate Dirac op
void invertXconj1dComplexSrc(std::vector<FermionField2f> &out, XconjugateMobiusFermionD &action, ComplexD phi){
  GridBase* UGrid = action.GaugeGrid();
  GridBase* FGrid = action.FermionGrid();
  GridBase* FrbGrid = action.FermionRedBlackGrid();

  ConjugateGradient<FermionField1f> cg(1e-08,10000);
  SchurRedBlackDiagMooeeSolve<FermionField1f> solver(cg);

  typedef FermionField1f::scalar_object Spinor;

  LatticeInteger tcoor4d(UGrid);
  LatticeCoordinate(tcoor4d,3);

  FermionField1f tmp4d(UGrid);
  FermionField1f src4d(UGrid);
  FermionField1f src5d(FGrid), sol5d(FGrid);
  FermionField1f src5d_e(FrbGrid), src5d_o(FrbGrid), sol5d_o(FrbGrid);

  FermionField2f tmp2f(UGrid);

  SCFmatrixField V_4d(UGrid);

  static Gamma C = Gamma(Gamma::Algebra::MinusGammaY) * Gamma(Gamma::Algebra::GammaT);
  static Gamma g5 = Gamma(Gamma::Algebra::Gamma5);
  static Gamma X = C*g5;
  static GparityFlavour sigma1 = GparityFlavour(GparityFlavour::Algebra::SigmaX);
  

  for(int s=0;s<4;s++){
    for(int c=0;c<3;c++){
      Spinor vb = Zero();
      vb()(s)(c) = phi;
      
      Spinor vplus = mulPplusLeft(vb);
      Spinor vminus = mulPminusLeft(vb);

      Spinor* vpm[2] = {&vplus, &vminus};

      for(int pm=0;pm<2;pm++){
	tmp4d = *vpm[pm];
	src4d = Zero();
	src4d = where( tcoor4d == Integer(0), tmp4d, src4d); //vector on every site on timeslice 0, zero elsewhere

	action.ImportPhysicalFermionSource(src4d, src5d);
	solver.RedBlackSource(action, src5d, src5d_e, src5d_o);
	solver.RedBlackSolve(action, src5d_o, sol5d_o);
	solver.RedBlackSolution(action, sol5d_o, src5d_e, sol5d);
	action.ExportPhysicalFermionSolution(sol5d, tmp4d);

	//Generate 2f X-conjugate output
	PokeIndex<GparityFlavourIndex>(tmp2f, tmp4d, 0);
	tmp4d = -(X*conjugate(tmp4d));
	PokeIndex<GparityFlavourIndex>(tmp2f, tmp4d, 1);

	pokeSpinColorFlavorColumn(V_4d, tmp2f, s,c,pm);
      }
    }
  }

  SCFmatrixField tmp = mulPplusRight(V_4d) + mulPminusRight(V_4d)*sigma1;  
  SCFmatrixField Minv = mulURight(tmp);

  out = std::vector<FermionField2f>(24, UGrid);

  for(int fc=0;fc<2;fc++){
    for(int sc=0;sc<4;sc++){
      for(int cc=0;cc<3;cc++){
	int out_idx = cc+3*(sc + 4*fc);
	peekSpinColorFlavorColumn(out[out_idx], Minv, sc, cc, fc);
      }
    }
  }
}



int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  RealD tol = 1e-10;
  for(int i=1;i<argc;i++){
    std::string sarg(argv[i]);
    if(sarg == "-tol"){
      std::stringstream ss; ss << argv[i+1]; ss >> tol;
      std::cout << "Set tolerance to " << tol << std::endl;
    }
  }
  
  const int Ls=4;

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

  Gamma C = Gamma(Gamma::Algebra::MinusGammaY) * Gamma(Gamma::Algebra::GammaT);
  Gamma g5 = Gamma(Gamma::Algebra::Gamma5);
  Gamma X = C*g5;
  
  //Set up a regular MDWF action instance as well as X-conj and Xbar-conj versions
  RealD mass=0.04;
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

  FermionField2f tmp2f(UGrid);
  std::cout.precision(12);

  {
    std::cout << "Tests with real source" << std::endl;
    std::cout << "-------------------------------------------------" << std::endl;
  
    std::vector<FermionField2f> sol_gp;
    invertGparity(sol_gp, reg_action);

    std::vector<FermionField2f> sol_Xconj2d;
    invertXconj2d(sol_Xconj2d, reg_action);

    std::vector<FermionField2f> sol_Xconj1d;
    invertXconj1d(sol_Xconj1d, xconj_action);

    std::vector<FermionField2f> sol_Xconj1d_matrix;
    invertXconj1d_matrix(sol_Xconj1d_matrix, xconj_action);

    std::cout << "Comparisons for real source" << std::endl;
    for(int fc=0;fc<2;fc++){
      for(int sc=0;sc<4;sc++){
	for(int cc=0;cc<3;cc++){
	  int out_idx = cc+3*(sc + 4*fc);
	  tmp2f = sol_gp[out_idx] - sol_Xconj2d[out_idx];
	  RealD nrm2d = norm2(tmp2f);
	  tmp2f = sol_gp[out_idx] - sol_Xconj1d[out_idx];
	  RealD nrm1d = norm2(tmp2f);
	  tmp2f = sol_gp[out_idx] - sol_Xconj1d_matrix[out_idx];
	  RealD nrm1dm = norm2(tmp2f);

	  std::cout << fc << " " << sc << " " << cc << " " << nrm2d << " " << nrm1d << " " << nrm1dm << std::endl;
	  assert(nrm2d < tol);
	  assert(nrm1d < tol);
	  assert(nrm1dm < tol);
	}
      }
    }
  }

  {
    std::cout << "Tests with complex source" << std::endl;
    std::cout << "-------------------------------------------------" << std::endl;
  
    ComplexD phi(1.234,  -6.444); //"random" phase

    std::vector<FermionField2f> sol_gp;
    invertGparityComplexSrc(sol_gp, reg_action, phi);

    std::vector<FermionField2f> sol_Xconj1d_matrix;
    invertXconj1dComplexSrc(sol_Xconj1d_matrix, xconj_action, phi);

    std::cout << "Comparisons for complex source" << std::endl;
    for(int fc=0;fc<2;fc++){
      for(int sc=0;sc<4;sc++){
	for(int cc=0;cc<3;cc++){
	  int out_idx = cc+3*(sc + 4*fc);
	  tmp2f = sol_gp[out_idx] - sol_Xconj1d_matrix[out_idx];
	  RealD nrm1dm = norm2(tmp2f);

	  std::cout << fc << " " << sc << " " << cc << " " << nrm1dm << std::endl;
	  assert(nrm1dm < tol);
	}
      }
    }
  }

  std::cout << "All tests passed" << std::endl;
  
  Grid_finalize();
}
