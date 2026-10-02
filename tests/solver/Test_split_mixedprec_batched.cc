/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./tests/solver/Test_split_mixedprec_batched.cc

    Copyright (C) 2026

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
/////////////////////////////////////////////////////////////////////////////////////////////
// MixedPrecisionConjugateGradientBatched with inner solves on split-communicator
// partitions (--batched-solver-split), against the same solver unsplit.
//
//   mpirun -n 2 ./Test_split_mixedprec_batched --grid 8.8.8.8 --mpi 1.1.1.2 --batched-solver-split 1.1.1.1
//
// With no --batched-solver-split the partitions are single ranks (1.1.1.1).
// NBATCH is deliberately not a multiple of the partition count, so the last group is
// zero-padded.
/////////////////////////////////////////////////////////////////////////////////////////////
#include <Grid/Grid.h>

using namespace std;
using namespace Grid;

const int   NBATCH    = 3;
const RealD TOLERANCE = 1.0e-8;

template<class Field>
RealD RelativeDifference(const Field &a,const Field &b)
{
  Field diff(a.Grid());
  diff = a - b;
  return std::sqrt(norm2(diff)/norm2(b));
}

/////////////////////////////////////////////////////////////////////////////////////////////
// The solver assumes vector p of a group lands in the partition whose
// GridSplitVectorIndex is p: fill vector p with the constant p+1, split, and check
// every partition sees its own index.
/////////////////////////////////////////////////////////////////////////////////////////////
void CheckPartitionOrder(GridCartesian *UGrid,const Coordinate &layout)
{
  GridCartesian SGrid(UGrid->FullDimensions(),UGrid->_simd_layout,layout,*UGrid);
  int P     = UGrid->ProcessorCount()/SGrid.ProcessorCount();
  int index = GridSplitVectorIndex(UGrid,&SGrid);

  std::vector<LatticeComplexF> full(P,UGrid);
  for(int p=0;p<P;p++){
    full[p] = ComplexF(p+1,0.0);
  }
  LatticeComplexF split(&SGrid);
  Grid_split(full,split);

  LatticeComplexF expect(&SGrid);
  expect = ComplexF(index+1,0.0);
  RealD err = norm2(split - expect);

  std::cout << GridLogMessage << "Partition order: vector index " << index << " error " << err << std::endl;
  GRID_ASSERT(err == 0.0);
}

/////////////////////////////////////////////////////////////////////////////////////////////
// For one pair of double/float linear operators:
//  (1) red-black-aware split/unsplit round trip is exact;
//  (2) the split clone's HermOp agrees with the full-grid HermOp;
//  (3) the batched solver with split inner solves agrees with the unsplit solver;
//  (4) the split solution has true residual at the requested tolerance.
/////////////////////////////////////////////////////////////////////////////////////////////
template<class FieldD,class FieldF>
void CheckSplitBatched(const std::string &name,
                       LinearOperatorBase<FieldD> &Linop_d,
                       LinearOperatorBase<FieldF> &Linop_f,
                       GridBase *grid_f,
                       std::vector<FieldD> &src,
                       const Coordinate &layout)
{
  std::cout << GridLogMessage << "==================================================" << std::endl;
  std::cout << GridLogMessage << name << std::endl;
  std::cout << GridLogMessage << "==================================================" << std::endl;

  int cb = src[0].Checkerboard();

  SplitOperator<FieldF> *split = Linop_f.SplitClone(layout);
  GRID_ASSERT(split != nullptr);
  int P = split->Partitions;

  std::vector<FieldF> full(P,grid_f);
  std::vector<FieldF> back(P,grid_f);
  std::vector<FieldF> Mfull(P,grid_f);
  for(int p=0;p<P;p++){
    full[p].Checkerboard()  = cb;
    back[p].Checkerboard()  = cb;
    Mfull[p].Checkerboard() = cb;
    precisionChange(full[p],src[p%src.size()]);
  }

  // (1) round trip
  FieldF s_in(split->FieldGrid);
  FieldF s_out(split->FieldGrid);
  Grid_split(full,s_in);
  Grid_unsplit(back,s_in);
  for(int p=0;p<P;p++){
    RealD err = norm2(back[p] - full[p]);
    std::cout << GridLogMessage << name << ": split/unsplit round trip " << p << " error " << err << std::endl;
    GRID_ASSERT(err == 0.0);
  }

  // (2) operator
  split->Linop->HermOp(s_in,s_out);
  Grid_unsplit(back,s_out);
  for(int p=0;p<P;p++){
    Linop_f.HermOp(full[p],Mfull[p]);
    RealD err = RelativeDifference(back[p],Mfull[p]);
    std::cout << GridLogMessage << name << ": split HermOp vs full HermOp " << p << " relative difference " << err << std::endl;
    GRID_ASSERT(err < 1.0e-6);
  }
  delete split;

  // (3) solver, unsplit then split
  int NBatch = src.size();
  std::vector<FieldD> sol_ref(NBatch,src[0].Grid());
  std::vector<FieldD> sol_split(NBatch,src[0].Grid());
  for(int i=0;i<NBatch;i++){
    sol_ref[i].Checkerboard()   = cb;
    sol_split[i].Checkerboard() = cb;
    sol_ref[i]   = Zero();
    sol_split[i] = Zero();
  }

  MixedPrecisionConjugateGradientBatched<FieldD,FieldF> mCG(TOLERANCE,10000,50,1000,grid_f,Linop_f,Linop_d);

  std::cout << GridLogMessage << name << ": unsplit batched solve" << std::endl;
  mCG.BatchedSplit     = Coordinate();
  mCG.BatchedSplitNode = false;
  mCG(src,sol_ref);

  std::cout << GridLogMessage << name << ": split batched solve, partition layout " << layout << std::endl;
  mCG.BatchedSplit = layout;
  mCG(src,sol_split);

  FieldD Msol(src[0].Grid());
  Msol.Checkerboard() = cb;
  for(int i=0;i<NBatch;i++){
    RealD diff = RelativeDifference(sol_split[i],sol_ref[i]);
    Linop_d.HermOp(sol_split[i],Msol);
    RealD resid = RelativeDifference(Msol,src[i]);
    std::cout << GridLogMessage << name << ": rhs " << i
              << " split vs unsplit " << diff
              << " true residual " << resid << std::endl;
    GRID_ASSERT(diff  < 1.0e-6);
    GRID_ASSERT(resid < 10.0*TOLERANCE);
  }
}

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  const int Ls = 8;

  Coordinate layout = GridDefaultBatchedSolverSplit();
  if ( layout.size() == 0 ) {
    layout = Coordinate(Nd,1);
  }

  GridCartesian         *UGrid_d   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplexD::Nsimd()), GridDefaultMpi());
  GridRedBlackCartesian *UrbGrid_d = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid_d);
  GridCartesian         *FGrid_d   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid_d);
  GridRedBlackCartesian *FrbGrid_d = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid_d);

  GridCartesian         *UGrid_f   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplexF::Nsimd()), GridDefaultMpi());
  GridRedBlackCartesian *UrbGrid_f = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid_f);
  GridCartesian         *FGrid_f   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid_f);
  GridRedBlackCartesian *FrbGrid_f = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid_f);

  CheckPartitionOrder(UGrid_f,layout);

  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});
  GridParallelRNG RNG4(UGrid_d);
  GridParallelRNG RNG5(FGrid_d);
  RNG4.SeedFixedIntegers(seeds4);
  RNG5.SeedFixedIntegers(seeds5);

  LatticeGaugeFieldD Umu_d(UGrid_d);
  LatticeGaugeFieldF Umu_f(UGrid_f);
  SU<Nc>::HotConfiguration(RNG4,Umu_d);
  precisionChange(Umu_f,Umu_d);

  // Antiperiodic in time for some actions, so boundary phases must travel with the links
  WilsonImplParams antiperiodic;
  antiperiodic.boundary_phases[Nd-1] = -1.0;

  // Sources: odd checkerboard (Schur) and full lattice (MdagM), 4d and 5d
  std::vector<LatticeFermionD> src4_o(NBATCH,UrbGrid_d);
  std::vector<LatticeFermionD> src5_o(NBATCH,FrbGrid_d);
  std::vector<LatticeFermionD> src5(NBATCH,FGrid_d);
  LatticeFermionD tmp4(UGrid_d);
  LatticeFermionD tmp5(FGrid_d);
  for(int i=0;i<NBATCH;i++){
    random(RNG4,tmp4);
    random(RNG5,tmp5);
    pickCheckerboard(Odd,src4_o[i],tmp4);
    pickCheckerboard(Odd,src5_o[i],tmp5);
    src5[i] = tmp5;
  }

  //////////////////////////////////////////
  // Wilson
  //////////////////////////////////////////
  {
    RealD mass = 0.1;
    WilsonFermionD Dd(Umu_d,*UGrid_d,*UrbGrid_d,mass);
    WilsonFermionF Df(Umu_f,*UGrid_f,*UrbGrid_f,mass);
    SchurDiagMooeeOperator<WilsonFermionD,LatticeFermionD> Ld(Dd);
    SchurDiagMooeeOperator<WilsonFermionF,LatticeFermionF> Lf(Df);
    CheckSplitBatched("Wilson SchurDiagMooee",Ld,Lf,UrbGrid_f,src4_o,layout);
  }

  //////////////////////////////////////////
  // Wilson clover
  //////////////////////////////////////////
  {
    RealD mass  = 0.1;
    RealD csw_r = 1.0;
    RealD csw_t = 1.0;
    WilsonCloverFermionD Dd(Umu_d,*UGrid_d,*UrbGrid_d,mass,csw_r,csw_t);
    WilsonCloverFermionF Df(Umu_f,*UGrid_f,*UrbGrid_f,mass,csw_r,csw_t);
    SchurDiagMooeeOperator<WilsonCloverFermionD,LatticeFermionD> Ld(Dd);
    SchurDiagMooeeOperator<WilsonCloverFermionF,LatticeFermionF> Lf(Df);
    CheckSplitBatched("WilsonClover SchurDiagMooee",Ld,Lf,UrbGrid_f,src4_o,layout);
  }

  //////////////////////////////////////////
  // Compact Wilson clover, antiperiodic
  //////////////////////////////////////////
  {
    RealD mass  = 0.1;
    RealD csw_r = 1.0;
    RealD csw_t = 1.0;
    RealD cF    = 1.0;
    WilsonAnisotropyCoefficients anis;
    CompactWilsonCloverFermionD Dd(Umu_d,*UGrid_d,*UrbGrid_d,mass,csw_r,csw_t,cF,anis,antiperiodic);
    CompactWilsonCloverFermionF Df(Umu_f,*UGrid_f,*UrbGrid_f,mass,csw_r,csw_t,cF,anis,antiperiodic);
    SchurDiagMooeeOperator<CompactWilsonCloverFermionD,LatticeFermionD> Ld(Dd);
    SchurDiagMooeeOperator<CompactWilsonCloverFermionF,LatticeFermionF> Lf(Df);
    CheckSplitBatched("CompactWilsonClover SchurDiagMooee antiperiodic",Ld,Lf,UrbGrid_f,src4_o,layout);
  }

  //////////////////////////////////////////
  // Domain wall
  //////////////////////////////////////////
  {
    RealD mass = 0.1;
    RealD M5   = 1.8;
    DomainWallFermionD Dd(Umu_d,*FGrid_d,*FrbGrid_d,*UGrid_d,*UrbGrid_d,mass,M5);
    DomainWallFermionF Df(Umu_f,*FGrid_f,*FrbGrid_f,*UGrid_f,*UrbGrid_f,mass,M5);
    SchurDiagMooeeOperator<DomainWallFermionD,LatticeFermionD> Ld(Dd);
    SchurDiagMooeeOperator<DomainWallFermionF,LatticeFermionF> Lf(Df);
    CheckSplitBatched("DomainWall SchurDiagMooee",Ld,Lf,FrbGrid_f,src5_o,layout);
  }

  //////////////////////////////////////////
  // Mobius, antiperiodic, with unequal
  // masses; all wrapper kinds
  //////////////////////////////////////////
  {
    RealD mass = 0.1;
    RealD M5   = 1.8;
    RealD b    = 1.5;
    RealD c    = 0.5;
    MobiusFermionD Dd(Umu_d,*FGrid_d,*FrbGrid_d,*UGrid_d,*UrbGrid_d,mass,M5,b,c,antiperiodic);
    MobiusFermionF Df(Umu_f,*FGrid_f,*FrbGrid_f,*UGrid_f,*UrbGrid_f,mass,M5,b,c,antiperiodic);
    Dd.SetMass(0.1,0.12);
    Df.SetMass(0.1,0.12);
    {
      SchurDiagMooeeOperator<MobiusFermionD,LatticeFermionD> Ld(Dd);
      SchurDiagMooeeOperator<MobiusFermionF,LatticeFermionF> Lf(Df);
      CheckSplitBatched("Mobius SchurDiagMooee antiperiodic",Ld,Lf,FrbGrid_f,src5_o,layout);
    }
    {
      SchurDiagOneOperator<MobiusFermionD,LatticeFermionD> Ld(Dd);
      SchurDiagOneOperator<MobiusFermionF,LatticeFermionF> Lf(Df);
      CheckSplitBatched("Mobius SchurDiagOne antiperiodic",Ld,Lf,FrbGrid_f,src5_o,layout);
    }
    {
      SchurDiagTwoOperator<MobiusFermionD,LatticeFermionD> Ld(Dd);
      SchurDiagTwoOperator<MobiusFermionF,LatticeFermionF> Lf(Df);
      CheckSplitBatched("Mobius SchurDiagTwo antiperiodic",Ld,Lf,FrbGrid_f,src5_o,layout);
    }
    {
      MdagMLinearOperator<MobiusFermionD,LatticeFermionD> Ld(Dd);
      MdagMLinearOperator<MobiusFermionF,LatticeFermionF> Lf(Df);
      CheckSplitBatched("Mobius MdagM antiperiodic",Ld,Lf,FGrid_f,src5,layout);
    }
  }

  //////////////////////////////////////////
  // ZMobius
  //////////////////////////////////////////
  {
    RealD mass = 0.1;
    RealD M5   = 1.8;
    RealD b    = 1.0;
    RealD c    = 0.0;
    std::vector<ComplexD> gamma(Ls);
    for(int s=0;s<Ls;s++){
      gamma[s] = ComplexD(1.0+0.05*s, (s%2) ? 0.02 : -0.02);
    }
    ZMobiusFermionD Dd(Umu_d,*FGrid_d,*FrbGrid_d,*UGrid_d,*UrbGrid_d,mass,M5,gamma,b,c);
    ZMobiusFermionF Df(Umu_f,*FGrid_f,*FrbGrid_f,*UGrid_f,*UrbGrid_f,mass,M5,gamma,b,c);
    SchurDiagMooeeOperator<ZMobiusFermionD,LatticeFermionD> Ld(Dd);
    SchurDiagMooeeOperator<ZMobiusFermionF,LatticeFermionF> Lf(Df);
    CheckSplitBatched("ZMobius SchurDiagMooee",Ld,Lf,FrbGrid_f,src5_o,layout);
  }

  //////////////////////////////////////////
  // A derived operator must not inherit
  // its parent's SplitClone
  //////////////////////////////////////////
  {
    RealD mass = 0.1;
    RealD mu   = 0.1;
    WilsonTMFermionF Df(Umu_f,*UGrid_f,*UrbGrid_f,mass,mu);
    SchurDiagMooeeOperator<WilsonTMFermionF,LatticeFermionF> Lf(Df);
    SplitOperator<LatticeFermionF> *split = Lf.SplitClone(layout);
    std::cout << GridLogMessage << "WilsonTM SplitClone refused: " << (split == nullptr) << std::endl;
    GRID_ASSERT(split == nullptr);
  }

  std::cout << GridLogMessage << "Test_split_mixedprec_batched: all checks passed" << std::endl;

  Grid_finalize();
}
