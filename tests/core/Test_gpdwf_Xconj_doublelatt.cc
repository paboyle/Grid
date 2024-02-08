    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/core/Test_gpdwf_Xconj_doublelatt.cc

    Copyright (C) 2015

Author: Christopher Kelly <ckelly@bnl.gov>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: Peter Boyle <peterboyle@Peters-MacBook-Pro-2.local>
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

/**
   This test uses the relationship between G-parity BCs, X-conjugate BCs and antiperiodic BCs on a doubled lattice
   to perform a complete test of the implementations
 **/

#include <Grid/Grid.h>

using namespace std;
using namespace Grid;
 ;

//typedef GparityDomainWallFermionD GparityDiracOp;
//typedef DomainWallFermionD StandardDiracOp;
//typedef XconjugateDomainWallFermionD XconjDiracOp;
//#define DOP_PARAMS

typedef GparityMobiusFermionD GparityDiracOp;
typedef MobiusFermionD StandardDiracOp;
typedef XconjugateMobiusFermionD XconjDiracOp;
#define DOP_PARAMS ,1.5, 0.5


typedef typename GparityDiracOp::FermionField GparityFermionField;
typedef typename GparityDiracOp::GaugeField GparityGaugeField;
typedef typename GparityFermionField::vector_type vComplexType;

typedef typename StandardDiracOp::FermionField StandardFermionField;
typedef typename StandardDiracOp::GaugeField StandardGaugeField;

enum{ same_vComplex = std::is_same<vComplexType, typename StandardFermionField::vector_type>::value };
static_assert(same_vComplex == 1, "Dirac Operators must have same underlying SIMD complex type");

const Gamma & Xmatrix(){
  static Gamma C = Gamma(Gamma::Algebra::MinusGammaY) * Gamma(Gamma::Algebra::GammaT);
  static Gamma g5 = Gamma(Gamma::Algebra::Gamma5);
  static Gamma X = C*g5;
  return X;
}

void boostXconjToGparity(GparityFermionField &out, const StandardFermionField &in){
  PokeIndex<GparityFlavourIndex>(out, in, 0);
  StandardFermionField tmp = -(Xmatrix()*conjugate(in));
  PokeIndex<GparityFlavourIndex>(out, tmp, 1);
}
void boostXconjToDoubleLatt(StandardFermionField &out, const StandardFermionField &in, const int nu, const int L){
  assert(out.Grid() != in.Grid());
  assert(!out.Grid()->_isCheckerBoarded);
  assert(!in.Grid()->_isCheckerBoarded);

  LatticeInteger    xcoor_2L_5(out.Grid());
  LatticeCoordinate(xcoor_2L_5,1+nu); //note '1+nu'! This is because for 5D fields the s-direction is direction 0
  
  StandardFermionField tmp_2L(out.Grid());
  Replicate(in, out);
  tmp_2L = -(Xmatrix()*conjugate(out));
  out = where( xcoor_2L_5 >= Integer(L), tmp_2L, out );
}
void boostGparityToDoubleLatt(StandardFermionField &out, const GparityFermionField &in, const int nu, const int L){
  assert(out.Grid() != in.Grid());
  assert(!out.Grid()->_isCheckerBoarded);
  assert(!in.Grid()->_isCheckerBoarded);

  LatticeInteger    xcoor_2L_5(out.Grid());
  LatticeCoordinate(xcoor_2L_5,1+nu); //note '1+nu'! This is because for 5D fields the s-direction is direction 0

  StandardFermionField tmp_L(in.Grid());
  tmp_L = PeekIndex<GparityFlavourIndex>(in,0);
  Replicate(tmp_L, out);

  tmp_L = PeekIndex<GparityFlavourIndex>(in,1);
  StandardFermionField tmp_2L(out.Grid());
  Replicate(tmp_L, tmp_2L);

  out = where( xcoor_2L_5 >= Integer(L), tmp_2L, out );
}

//Break up a larger lattice into smaller lattices
//output vector is resized to the number of subdivisions, with lexicographic block indexing x+rx*(y+ry*(z+rz*t))
//where rx,ry,rz,rt are the grid size ratios
template<class vobj>
void SubDivide(std::vector<Lattice<vobj> > &out, const Lattice<vobj> &in, GridBase* out_grid)
{
  out.resize(0,out_grid);
  typedef typename vobj::scalar_object sobj;

  GridBase *og = out_grid;
  GridBase *ig = in.Grid();

  int nd = ig->_ndimension;
  assert(og->_ndimension==nd);

  subdivides(og,ig); 

  int nout = 1;
  Coordinate ratio(nd);
  for(int d=0;d<nd;d++){
    ratio[d] = ig->_fdimensions[d]/og->_fdimensions[d];
    nout *= ratio[d];
  }
  out.resize(nout, out_grid);
  
  Coordinate ocoor(nd);
  Coordinate icoor(nd);
  for(int g=0;g<ig->gSites();g++){

    int oidx = 0;
    ig->GlobalIndexToGlobalCoor(g,icoor);
    for(int d=nd-1;d>=0;d--){
      ocoor[d] = icoor[d]%og->_gdimensions[d];
      oidx = icoor[d] / og->_gdimensions[d] + ratio[d]*oidx; //lex mapping  x+rx*(y+ry*(z+rz*t))
    }
    
    sobj tmp;
    peekSite(tmp,in,icoor);
    pokeSite(tmp,out[oidx],ocoor);
  }

}



int main (int argc, char ** argv)
{
  int nu = 0;
  int tbc_aprd = 0; //use antiperiodic BCs in the time direction?
  
  Grid_init(&argc,&argv);

  for(int i=1;i<argc;i++){
    if(std::string(argv[i]) == "--Gparity-dir"){
      std::stringstream ss; ss << argv[i+1]; ss >> nu;
      std::cout << GridLogMessage << "Set Gparity direction to " << nu << std::endl;
    }else if(std::string(argv[i]) == "--Tbc-APRD"){
      tbc_aprd = 1;
      std::cout << GridLogMessage << "Using antiperiodic BCs in the time direction" << std::endl;
    }
  }

  std::cout << GridLogMessage<< "*****************************************************************" <<std::endl;
  std::cout << GridLogMessage<< "* Kernel options --dslash-generic, --dslash-unroll, --dslash-asm" <<std::endl;
  std::cout << GridLogMessage<< "*****************************************************************" <<std::endl;
  std::cout << GridLogMessage<< "*****************************************************************" <<std::endl;
  std::cout << GridLogMessage<< "* Testing X-conjugate and Gparity Dirac operator                  "<<std::endl;
  std::cout << GridLogMessage<< "* Vectorising space-time by "<<vComplexType::Nsimd()<<std::endl;
#ifdef GRID_OMP
  if ( WilsonKernelsStatic::Comms == WilsonKernelsStatic::CommsAndCompute ) std::cout << GridLogMessage<< "* Using Overlapped Comms/Compute" <<std::endl;
  if ( WilsonKernelsStatic::Comms == WilsonKernelsStatic::CommsThenCompute) std::cout << GridLogMessage<< "* Using sequential comms compute" <<std::endl;
#endif
  if ( WilsonKernelsStatic::Opt == WilsonKernelsStatic::OptGeneric   ) std::cout << GridLogMessage<< "* Using GENERIC Nc WilsonKernels" <<std::endl;
  if ( WilsonKernelsStatic::Opt == WilsonKernelsStatic::OptHandUnroll) std::cout << GridLogMessage<< "* Using UNROLLED Nc=3       WilsonKernels" <<std::endl;
  if ( WilsonKernelsStatic::Opt == WilsonKernelsStatic::OptInlineAsm ) std::cout << GridLogMessage<< "* Using Asm Nc=3   WilsonKernels" <<std::endl;
  std::cout << GridLogMessage<< "*****************************************************************" <<std::endl;

  const int Ls=4;
  Coordinate latt_L = GridDefaultLatt();
  Coordinate latt_2L(latt_L); latt_2L[nu] = 2*latt_L[nu];
  int L = latt_L[nu];

  Coordinate simd_layout = GridDefaultSimd(Nd,vComplexType::Nsimd());
  Coordinate mpi_layout  = GridDefaultMpi(); //node layout

  GridCartesian         * UGrid_L   = SpaceTimeGrid::makeFourDimGrid(latt_L, simd_layout, mpi_layout);
  GridRedBlackCartesian * UrbGrid_L = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid_L);
  GridCartesian         * FGrid_L   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid_L);
  GridRedBlackCartesian * FrbGrid_L = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid_L);


  GridCartesian         * UGrid_2L   = SpaceTimeGrid::makeFourDimGrid(latt_2L, simd_layout, mpi_layout);
  GridRedBlackCartesian * UrbGrid_2L = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid_2L);
  GridCartesian         * FGrid_2L   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid_2L);
  GridRedBlackCartesian * FrbGrid_2L = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid_2L);

  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});
  GridParallelRNG          RNG5_L(FGrid_L);  RNG5_L.SeedFixedIntegers(seeds5);
  GridParallelRNG          RNG4_L(UGrid_L);  RNG4_L.SeedFixedIntegers(seeds4);

  GparityGaugeField Umu_L(UGrid_L); //note type is same as StandardGaugeField
  SU<Nc>::HotConfiguration(RNG4_L,Umu_L);
  
  ////////////////////////
  //Copy-conjugate the gauge field
  ////////////////////////
  StandardGaugeField Umu_2L(UGrid_2L); 
  Replicate(Umu_L,Umu_2L); //output on right, grr!

  LatticeInteger xcoor_2L_4(UGrid_2L);
  LatticeCoordinate(xcoor_2L_4,nu);
  StandardGaugeField Umu_2L_conj = conjugate(Umu_2L);
  Umu_2L = where( xcoor_2L_4 >= Integer(L), Umu_2L_conj, Umu_2L ); //don't need Cshift as replicate already duplicates U onto the doubled lattice
 
  ////////////////////////
  //Generate two one-flavor fields with implicit X-conjugate BCs
  ////////////////////////

  StandardFermionField src_X_L(FGrid_L), src2_X_L(FGrid_L);
  random(RNG5_L,src_X_L);
  random(RNG5_L,src2_X_L);
   
  ////////////////////////
  //Create the corresponding two-flavor G-parity fields
  ////////////////////////
  GparityFermionField src_GP_L(FGrid_L), src2_GP_L(FGrid_L);
  boostXconjToGparity(src_GP_L,src_X_L);
  boostXconjToGparity(src2_GP_L,src2_X_L);

  ////////////////////////
  //Create the corresponding doubled-lattice fields
  ////////////////////////
  StandardFermionField src_std_2L(FGrid_2L), src2_std_2L(FGrid_2L);
  boostXconjToDoubleLatt(src_std_2L, src_X_L, nu, L);
  boostXconjToDoubleLatt(src2_std_2L, src2_X_L, nu, L);

  ////////////////////////
  //Create the Dirac operators
  ////////////////////////
  RealD mass=0.01;
  RealD M5=1.8;

  //Standard Dirac op on doubled lattice
  AcceleratorVector<Complex,4> bc_std(Nd, 1.0);
  if(tbc_aprd) bc_std[Nd-1] = -1.; //antiperiodic time BC
  bc_std[nu] = -1; //antiperiodic in G-parity direction
  StandardDiracOp::ImplParams std_params(bc_std);
  StandardDiracOp Dstd_2L(Umu_2L,*FGrid_2L,*FrbGrid_2L,*UGrid_2L,*UrbGrid_2L,mass,M5 DOP_PARAMS, std_params);

  SchurDiagMooeeOperator<StandardDiracOp,StandardFermionField> HermOp_std_2L(Dstd_2L);
  ConjugateGradient<StandardFermionField> CG_std(1.0e-8,10000);

  //G-parity Dirac op on single lattice
  std::vector<int> twists(Nd,0);
  twists[nu] = 1;
  if(tbc_aprd) twists[Nd-1] = 1;

  GparityDiracOp::ImplParams gp_params;
  gp_params.twists = twists;
  GparityDiracOp DGP_L(Umu_L,*FGrid_L,*FrbGrid_L,*UGrid_L,*UrbGrid_L,mass,M5 DOP_PARAMS,gp_params);
  
  SchurDiagMooeeOperator<GparityDiracOp,GparityFermionField> HermOp_GP_L(DGP_L);
  ConjugateGradient<GparityFermionField> CG_GP(1.0e-8,10000);

  //X-conj Dirac op on single lattice
  XconjDiracOp::ImplParams x_params;
  x_params.twists = twists;
  XconjDiracOp DX_L(Umu_L,*FGrid_L,*FrbGrid_L,*UGrid_L,*UrbGrid_L,mass,M5 DOP_PARAMS,x_params);
  
  SchurDiagMooeeOperator<XconjDiracOp,StandardFermionField> HermOp_X_L(DX_L);

  /////////////////////
  //Solve inverse with 1st source
  /////////////////////
  StandardFermionField    src_o_std_2L(FrbGrid_2L);
  StandardFermionField result_o_std_2L(FrbGrid_2L);
  pickCheckerboard(Odd,src_o_std_2L,src_std_2L);
  result_o_std_2L=Zero();
  CG_std(HermOp_std_2L,src_o_std_2L,result_o_std_2L);
  StandardFermionField result_std_2L(FGrid_2L);
  result_std_2L = Zero();
  setCheckerboard(result_std_2L,result_o_std_2L); 

  StandardFermionField    src_o_X_L(FrbGrid_L);
  StandardFermionField result_o_X_L(FrbGrid_L);
  pickCheckerboard(Odd,src_o_X_L,src_X_L);
  result_o_X_L=Zero();
  CG_std(HermOp_X_L,src_o_X_L,result_o_X_L);
  StandardFermionField result_X_L(FGrid_L);
  result_X_L = Zero();
  setCheckerboard(result_X_L,result_o_X_L); 

  GparityFermionField    src_o_GP_L(FrbGrid_L);
  GparityFermionField result_o_GP_L(FrbGrid_L);
  pickCheckerboard(Odd,src_o_GP_L,src_GP_L);
  result_o_GP_L=Zero();
  CG_GP(HermOp_GP_L,src_o_GP_L,result_o_GP_L);
  GparityFermionField result_GP_L(FGrid_L);
  result_GP_L = Zero();
  setCheckerboard(result_GP_L,result_o_GP_L); 

  std::cout << "Norms:  std:" << norm2(result_o_std_2L) << " X:" << 2*norm2(result_o_X_L) << " GP:" << norm2(result_o_GP_L) << std::endl;

  /////////////////////
  //Convert X-conj and GP solutions to 2L formulation for direct comparison
  ////////////////////
  std::cout << "Computing differences between solutions in 2L formulation" << std::endl;
  StandardFermionField result_X_2L(FGrid_2L);
  boostXconjToDoubleLatt(result_X_2L, result_X_L, nu, L);
  
  StandardFermionField result_GP_2L(FGrid_2L);
  boostGparityToDoubleLatt(result_GP_2L, result_GP_L, nu, L);
     
  StandardFermionField diff_2L(FGrid_2L);
  diff_2L = result_X_2L - result_std_2L;
  std::cout << "X vs std: " << norm2(diff_2L) << std::endl;
  assert(norm2(diff_2L) < 1e-8);
  diff_2L = result_GP_2L - result_std_2L;
  std::cout << "GP vs std: " << norm2(diff_2L) << std::endl;
  assert(norm2(diff_2L) < 1e-8);

  ////////////////////////////////
  //Repeat experiments with a G-parity source, back convert to X-conjugate via doubled lattice as in doc
  ////////////////////////////////
  std::cout << "Using a random two-flavor G-parity field" << std::endl;
  
  random(RNG5_L,src_GP_L);
  boostGparityToDoubleLatt(src_std_2L, src_GP_L, nu, L);
  std::vector<StandardFermionField> src_std_2L_split;
  SubDivide(src_std_2L_split, src_std_2L, FGrid_L);
  assert(src_std_2L_split.size() == 2);
  
  //v[0]  = a+ib
  //v[1]  = -X[a-ib]*
  StandardFermionField tmp_X_L(FGrid_L), tmp2_X_L(FGrid_L);
  tmp_X_L = Xmatrix()*conjugate(src_std_2L_split[1]);
  src_X_L = 0.5*( src_std_2L_split[0] + tmp_X_L ); //a
  src2_X_L = ComplexD(0,-0.5)*( src_std_2L_split[0] - tmp_X_L ); //b
  
  /////////////////////
  //Solve inverse with new sources
  /////////////////////
  pickCheckerboard(Odd,src_o_std_2L,src_std_2L);
  result_o_std_2L=Zero();
  CG_std(HermOp_std_2L,src_o_std_2L,result_o_std_2L);
  result_std_2L = Zero();
  setCheckerboard(result_std_2L,result_o_std_2L); 

  //note 2 inversions for Xconj
  pickCheckerboard(Odd,src_o_X_L,src_X_L);
  result_o_X_L=Zero();
  CG_std(HermOp_X_L,src_o_X_L,result_o_X_L);
  result_X_L = Zero();
  setCheckerboard(result_X_L,result_o_X_L); 

  StandardFermionField    src2_o_X_L(FrbGrid_L);
  StandardFermionField result2_o_X_L(FrbGrid_L);
  pickCheckerboard(Odd,src2_o_X_L,src2_X_L);
  result2_o_X_L=Zero();
  CG_std(HermOp_X_L,src2_o_X_L,result2_o_X_L);
  StandardFermionField result2_X_L(FGrid_L);
  result2_X_L = Zero();
  setCheckerboard(result2_X_L,result2_o_X_L); 

  //Gparity
  pickCheckerboard(Odd,src_o_GP_L,src_GP_L);
  result_o_GP_L=Zero();
  CG_GP(HermOp_GP_L,src_o_GP_L,result_o_GP_L);
  result_GP_L = Zero();
  setCheckerboard(result_GP_L,result_o_GP_L); 

  ///////////////////
  //Convert to doubled latt
  ///////////////////
  std::cout << "Computing differences between solutions in 2L formulation" << std::endl;
  boostXconjToDoubleLatt(result_X_2L, result_X_L, nu, L);
  //Note we need to reconstruct the full solution
  StandardFermionField tmp_X_2L(FGrid_2L);
  boostXconjToDoubleLatt(tmp_X_2L, result2_X_L, nu, L);
  result_X_2L = result_X_2L + ComplexD(0,1)*tmp_X_2L;

  boostGparityToDoubleLatt(result_GP_2L, result_GP_L, nu, L);
     
  diff_2L = result_X_2L - result_std_2L;
  std::cout << "X vs std: " << norm2(diff_2L) << std::endl;
  assert(norm2(diff_2L) < 1e-8);
  diff_2L = result_GP_2L - result_std_2L;
  std::cout << "GP vs std: " << norm2(diff_2L) << std::endl;
  assert(norm2(diff_2L) < 1e-8);

  Grid_finalize();
}
