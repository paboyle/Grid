    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: tests/hmc/Test_action_dwf_gparity2fvs1f.cc

    Copyright (C) 2015

    Author: Christopher Kelly <ckelly@bnl.gov>
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
#include <Grid/Grid.h>

using namespace Grid;



template<typename FermionField2f, typename FermionField1f>
void copy2fTo1fFermionField(FermionField1f &out, const FermionField2f &in, int gpdir){
  auto f0_halfgrid = PeekIndex<GparityFlavourIndex>(in,0); //on 2f Grid
  FermionField1f f0_fullgrid_dbl(out.Grid());
  Replicate(f0_halfgrid, f0_fullgrid_dbl); //double it up to live on the 1f Grid

  auto f1_halfgrid = PeekIndex<GparityFlavourIndex>(in,1);
  FermionField1f f1_fullgrid_dbl(out.Grid());
  Replicate(f1_halfgrid, f1_fullgrid_dbl);
  
  const Coordinate &dim_2f = in.Grid()->GlobalDimensions();
  const Coordinate &dim_1f = out.Grid()->GlobalDimensions();

  //We have to be careful for 5d fields; the s-direction is placed before the x,y,z,t and so we need to shift gpdir by 1
  std::cout << "gpdir " << gpdir << std::endl;

  gpdir+=1;
  std::cout << "gpdir for 5D fields " << gpdir << std::endl;

  std::cout << "dim_2f " << dim_2f << std::endl;
  std::cout << "dim_1f " << dim_1f << std::endl;
  
  assert(dim_1f[gpdir] == 2*dim_2f[gpdir]);

  LatticeInteger xcoor_1f(out.Grid()); //5d lattice integer
  LatticeCoordinate(xcoor_1f,gpdir);

  Integer L = dim_2f[gpdir];

  out = where(xcoor_1f < L, f0_fullgrid_dbl, f1_fullgrid_dbl);
}

//Both have the same field type
void copy2fTo1fGaugeField(LatticeGaugeField &out, const LatticeGaugeField &in, int gpdir){
  LatticeGaugeField U_dbl(out.Grid());
  Replicate(in, U_dbl);
  
  LatticeGaugeField Uconj_dbl = conjugate( U_dbl );

  const Coordinate &dim_2f = in.Grid()->GlobalDimensions();
  
  LatticeInteger xcoor_1f(out.Grid());
  LatticeCoordinate(xcoor_1f,gpdir);

  Integer L = dim_2f[gpdir];
  
  out = where(xcoor_1f < L, U_dbl, Uconj_dbl);
}


std::ostream & operator<<(std::ostream &os, const Coordinate &x){
  os << "(";
  for(int i=0;i<x.size();i++) os << x[i] <<  (i<x.size()-1 ? " " : "");
  os << ")";
  return os;
}


int main(int argc, char **argv) {
  using namespace Grid;
  
  Grid_init(&argc, &argv);
  int threads = GridThread::GetThreads();

  std::cout << GridLogMessage << "Grid is setup to use " << threads << " threads" << std::endl;

  int Ls = 16;

  Coordinate latt_2f = GridDefaultLatt();
  Coordinate simd_layout = GridDefaultSimd(Nd, vComplexD::Nsimd());
  Coordinate mpi_layout = GridDefaultMpi();

  int mu = 0; //Gparity direction

  Coordinate latt_1f = latt_2f;
  latt_1f[mu] *= 2;

  GridCartesian         * UGrid_1f   = SpaceTimeGrid::makeFourDimGrid(latt_1f, simd_layout, mpi_layout);
  GridRedBlackCartesian * UrbGrid_1f = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid_1f);
  GridCartesian         * FGrid_1f   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid_1f);
  GridRedBlackCartesian * FrbGrid_1f = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid_1f);


  GridCartesian         * UGrid_2f   = SpaceTimeGrid::makeFourDimGrid(latt_2f, simd_layout, mpi_layout);
  GridRedBlackCartesian * UrbGrid_2f = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid_2f);
  GridCartesian         * FGrid_2f   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid_2f);
  GridRedBlackCartesian * FrbGrid_2f = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid_2f);


  std::cout << "SIMD layout " << simd_layout << std::endl;
  std::cout << "MPI layout " << mpi_layout << std::endl;
  std::cout << "2f dimensions " << latt_2f << std::endl;
  std::cout << "1f dimensions " << latt_1f << std::endl;

  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});
  GridParallelRNG          RNG5_2f(FGrid_2f);  RNG5_2f.SeedFixedIntegers(seeds5);
  GridParallelRNG          RNG4_2f(UGrid_2f);  RNG4_2f.SeedFixedIntegers(seeds4);

  std::cout << "Generating hot 2f gauge configuration" << std::endl;
  LatticeGaugeField Umu_2f(UGrid_2f);
  SU<Nc>::HotConfiguration(RNG4_2f,Umu_2f);

  std::cout << "Copying 2f->1f gauge field" << std::endl;
  LatticeGaugeField Umu_1f(UGrid_1f);
  copy2fTo1fGaugeField(Umu_1f, Umu_2f, mu);  

  typedef GparityWilsonImplR FermionImplPolicy2f;
  typedef GparityDomainWallFermionD FermionAction2f;
  typedef typename FermionAction2f::FermionField FermionField2f;
  
  typedef WilsonImplR FermionImplPolicy1f;
  typedef DomainWallFermionD FermionAction1f;
  typedef typename FermionAction1f::FermionField FermionField1f;

  std::cout << "Generating eta 2f" << std::endl;
  FermionField2f eta_2f(FGrid_2f);
  gaussian(RNG5_2f, eta_2f);

  RealD scale = std::sqrt(0.5);
  eta_2f=eta_2f*scale;

  std::cout << "Copying 2f->1f eta" << std::endl;
  FermionField1f eta_1f(FGrid_1f);
  copy2fTo1fFermionField(eta_1f, eta_2f, mu);
  
  Real beta         = 2.13;
  Real light_mass   = 0.01;
  Real strange_mass = 0.032;
  Real pv_mass      = 1.0;
  RealD M5  = 1.8;

  //Setup the Dirac operators
  std::cout << "Initializing Dirac operators" << std::endl;
  
  FermionAction2f::ImplParams Params_2f;
  Params_2f.twists[mu] = 1;
  Params_2f.twists[Nd-1] = 1; //APBC in time direction

  //note 'Num' and 'Den' here refer to the determinant ratio, not the operator ratio in the pseudofermion action where the two are inverted
  //to my mind the Pauli Villars and 'denominator' are synonymous but the Grid convention has this as the 'Numerator' operator in the RHMC implementation
  FermionAction2f NumOp_2f(Umu_2f,*FGrid_2f,*FrbGrid_2f,*UGrid_2f, *UrbGrid_2f, light_mass,M5,Params_2f); 
  FermionAction2f DenOp_2f(Umu_2f,*FGrid_2f,*FrbGrid_2f,*UGrid_2f, *UrbGrid_2f, pv_mass, M5,Params_2f);

  FermionAction1f::ImplParams Params_1f;
  Params_1f.boundary_phases[mu] = -1; //antiperiodic in doubled lattice in GP direction
  Params_1f.boundary_phases[Nd-1] = -1;
  
  FermionAction1f NumOp_1f(Umu_1f,*FGrid_1f,*FrbGrid_1f,*UGrid_1f, *UrbGrid_1f, light_mass,M5,Params_1f);
  FermionAction1f DenOp_1f(Umu_1f,*FGrid_1f,*FrbGrid_1f,*UGrid_1f, *UrbGrid_1f, pv_mass, M5,Params_1f);

  //Test the replication routines by running a CG on eta
  double StoppingCondition = 1e-10;
  double MaxCGIterations = 30000;
  ConjugateGradient<FermionField2f>  CG_2f(StoppingCondition,MaxCGIterations);
  ConjugateGradient<FermionField1f>  CG_1f(StoppingCondition,MaxCGIterations);

  NumOp_1f.ImportGauge(Umu_1f);
  NumOp_2f.ImportGauge(Umu_2f);

  FermionField1f test_1f(FGrid_1f);
  FermionField2f test_2f(FGrid_2f);
  
  MdagMLinearOperator<FermionAction1f, FermionField1f> Linop_1f(NumOp_1f);
  MdagMLinearOperator<FermionAction2f, FermionField2f> Linop_2f(NumOp_2f);
  
  CG_1f(Linop_1f, eta_1f, test_1f);
  CG_2f(Linop_2f, eta_2f, test_2f);
  RealD test_1f_norm = norm2(test_1f);
  RealD test_2f_norm = norm2(test_2f);

  std::cout << "Verification of replication routines: " << test_1f_norm << " " << test_2f_norm << " " << test_1f_norm - test_2f_norm << std::endl;


#if 1
  typedef GeneralEvenOddRatioRationalPseudoFermionAction<FermionImplPolicy2f> Action2f;
  typedef GeneralEvenOddRatioRationalPseudoFermionAction<FermionImplPolicy1f> Action1f;

  RationalActionParams rational_params;
  rational_params.inv_pow = 2;
  rational_params.lo = 1e-5;
  rational_params.hi = 32;
  rational_params.md_degree = 16;
  rational_params.action_degree = 16;

  Action2f action_2f(DenOp_2f, NumOp_2f, rational_params);
  Action1f action_1f(DenOp_1f, NumOp_1f, rational_params);
#else
  typedef TwoFlavourEvenOddRatioPseudoFermionAction<FermionImplPolicy2f> Action2f;
  typedef TwoFlavourEvenOddRatioPseudoFermionAction<FermionImplPolicy1f> Action1f;

  Action2f action_2f(DenOp_2f, NumOp_2f, CG_2f, CG_2f);
  Action1f action_1f(DenOp_1f, NumOp_1f, CG_1f, CG_1f);
#endif


  std::cout << "Action refresh" << std::endl;
  action_2f.refresh(Umu_2f, eta_2f);
  action_1f.refresh(Umu_1f, eta_1f);

  std::cout << "Action compute post heatbath" << std::endl;
  RealD S_2f = action_2f.S(Umu_2f);
  RealD S_1f = action_1f.S(Umu_1f);

  std::cout << "Action comparison post heatbath" << std::endl;
  std::cout << S_2f << " " << S_1f << " " << S_2f-S_1f << std::endl;

  //Change the gauge field between refresh and action eval else the matrix and inverse matrices all cancel and we just get |eta|^2
  SU<Nc>::HotConfiguration(RNG4_2f,Umu_2f);
  copy2fTo1fGaugeField(Umu_1f, Umu_2f, mu);  

  //Now compute the action with the new gauge field
  std::cout << "Action compute post gauge field update" << std::endl;
  S_2f = action_2f.S(Umu_2f);
  S_1f = action_1f.S(Umu_1f);

  std::cout << "Action comparison post gauge field update" << std::endl;
  std::cout << S_2f << " " << S_1f << " " << S_2f-S_1f << std::endl;

  Grid_finalize();
} // main


