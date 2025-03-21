/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_general_coarse_hdcg.cc

    Copyright (C) 2023

Author: Peter Boyle <pboyle@bnl.gov>

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
#include <Grid/algorithms/iterative/ImplicitlyRestartedBlockLanczos.h>
#include <Grid/algorithms/iterative/ImplicitlyRestartedBlockLanczosCoarse.h>
#include <Grid/algorithms/iterative/AdefMrhs.h>

using namespace std;
using namespace Grid;

template<class aggregation>
void SaveFineEvecs(aggregation &Agg,std::string file)
{
#ifdef HAVE_LIME
  emptyUserRecord record;
  ScidacWriter WR(Agg[0].Grid()->IsBoss());
  WR.open(file);
  for(int b=0;b<Agg.size();b++){
    WR.writeScidacFieldRecord(Agg[b],record,0,Grid::BinaryIO::BINARYIO_LEXICOGRAPHIC);
  }
  WR.close();
#endif
}
template<class aggregation>
void SaveBasis(aggregation &Agg,std::string file)
{
#ifdef HAVE_LIME
  emptyUserRecord record;
  ScidacWriter WR(Agg.FineGrid->IsBoss());
  WR.open(file);
  for(int b=0;b<Agg.subspace.size();b++){
    WR.writeScidacFieldRecord(Agg.subspace[b],record,0,Grid::BinaryIO::BINARYIO_LEXICOGRAPHIC);
    //    WR.writeScidacFieldRecord(Agg.subspace[b],record);
  }
  WR.close();
#endif
}
template<class aggregation>
void LoadBasis(aggregation &Agg, std::string file)
{
#ifdef HAVE_LIME
  emptyUserRecord record;
  ScidacReader RD ;
  RD.open(file);
  for(int b=0;b<Agg.subspace.size();b++){
    RD.readScidacFieldRecord(Agg.subspace[b],record,Grid::BinaryIO::BINARYIO_LEXICOGRAPHIC);
    //    RD.readScidacFieldRecord(Agg.subspace[b],record,0);
  }    
  RD.close();
#endif
}
template<class aggregation>
void LoadFineEvecs(aggregation &Agg, std::string file,LatticeFermionF & conv_tmp)
{
#ifdef HAVE_LIME
  emptyUserRecord record;
  ScidacReader RD ;
  RD.open(file);
  for(int b=0;b<Agg.size();b++){
    RD.readScidacFieldRecord(conv_tmp,record,Grid::BinaryIO::BINARYIO_LEXICOGRAPHIC);
    precisionChange(Agg[b],conv_tmp);
  }    
  RD.close();
#endif
}
template<class CoarseVector>
void SaveEigenvectors(std::vector<RealD>            &eval,
		      std::vector<CoarseVector>     &evec,
		      std::string evec_file,
		      std::string eval_file)
{
#ifdef HAVE_LIME
  emptyUserRecord record;
  ScidacWriter WR(evec[0].Grid()->IsBoss());
  WR.open(evec_file);
  for(int b=0;b<evec.size();b++){
    WR.writeScidacFieldRecord(evec[b],record,0,0);
  }
  WR.close();
  XmlWriter WRx(eval_file);
  write(WRx,"evals",eval);
#endif
}
template<class CoarseVector>
void LoadEigenvectors(std::vector<RealD>            &eval,
		      std::vector<CoarseVector>     &evec,
		      std::string evec_file,
		      std::string eval_file)
{
#ifdef HAVE_LIME
    XmlReader RDx(eval_file);
    read(RDx,"evals",eval);
    emptyUserRecord record;

    Grid::ScidacReader RD ;
    RD.open(evec_file);
    assert(evec.size()==eval.size());
    for(int k=0;k<eval.size();k++) {
      RD.readScidacFieldRecord(evec[k],record);
    }
    RD.close();
#endif
}

// Want Op in CoarsenOp to call MatPcDagMatPc
template<class Field>
class HermOpAdaptor : public LinearOperatorBase<Field>
{
  LinearOperatorBase<Field> & wrapped;
public:
  HermOpAdaptor(LinearOperatorBase<Field> &wrapme) : wrapped(wrapme)  {};
  void Op     (const Field &in, Field &out)   { wrapped.HermOp(in,out);  }
  void HermOp(const Field &in, Field &out)    { wrapped.HermOp(in,out); }
  void AdjOp     (const Field &in, Field &out){ wrapped.HermOp(in,out);  }
  void OpDiag (const Field &in, Field &out)                  {    assert(0);  }
  void OpDir  (const Field &in, Field &out,int dir,int disp) {    assert(0);  }
  void OpDirAll  (const Field &in, std::vector<Field> &out)  {    assert(0);  };
  void HermOpAndNorm(const Field &in, Field &out,RealD &n1,RealD &n2){    assert(0);  }
};

template<class Field> class CGSmoother : public LinearFunction<Field>
{
public:
  using LinearFunction<Field>::operator();
  typedef LinearOperatorBase<Field> FineOperator;
  FineOperator   & _SmootherOperator;
  int iters;
  CGSmoother(int _iters, FineOperator &SmootherOperator) :
    _SmootherOperator(SmootherOperator),
    iters(_iters)
  {
    std::cout << GridLogMessage<<" Mirs smoother order "<<iters<<std::endl;
  };
  void operator() (const Field &in, Field &out) 
  {
    ConjugateGradient<Field>  CG(0.0,iters,false); // non-converge is just fine in a smoother

    out=Zero();

    CG(_SmootherOperator,in,out);
  }
};


int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  const int Ls=24;
  const int nbasis = 62;
  const int cb = 0 ;
  RealD mass=0.00078;
  RealD M5=1.8;
  RealD b=1.5;
  RealD c=0.5;

  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(),
								   GridDefaultSimd(Nd,vComplex::Nsimd()),
								   GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  // Construct a coarsened grid with 4^4 cell
  Coordinate Block({4,4,6,4});
  Coordinate clatt = GridDefaultLatt();
  for(int d=0;d<clatt.size();d++){
    clatt[d] = clatt[d]/Block[d];
  }

  //////////////////////////////////////////
  // Double precision grids 
  //////////////////////////////////////////
  GridCartesian *Coarse4d =  SpaceTimeGrid::makeFourDimGrid(clatt,
							    GridDefaultSimd(Nd,vComplex::Nsimd()),
							    GridDefaultMpi());;
  GridCartesian *Coarse5d =  SpaceTimeGrid::makeFiveDimGrid(1,Coarse4d);

  //////////////////////////////////////////
  // Single precision grids -- lanczos + smoother
  //////////////////////////////////////////
  GridCartesian         * UGridF   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(),
								   GridDefaultSimd(Nd,vComplexF::Nsimd()),
								   GridDefaultMpi());
  GridRedBlackCartesian * UrbGridF = SpaceTimeGrid::makeFourDimRedBlackGrid(UGridF);
  GridCartesian         * FGridF   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGridF);
  GridRedBlackCartesian * FrbGridF = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGridF);
  ///////////////////////// RNGs /////////////////////////////////
  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});
  std::vector<int> cseeds({5,6,7,8});

  GridParallelRNG          RNG5(FGrid);   RNG5.SeedFixedIntegers(seeds5);
  GridParallelRNG          RNG4(UGrid);   RNG4.SeedFixedIntegers(seeds4);
  GridParallelRNG          CRNG(Coarse5d);CRNG.SeedFixedIntegers(cseeds);

  ///////////////////////// Configuration /////////////////////////////////
  LatticeGaugeField Umu(UGrid);

  FieldMetaData header;
  std::string file("ckpoint_lat.1000");
  NerscIO::readConfiguration(Umu,header,file);

  //////////////////////// Fermion action //////////////////////////////////
  MobiusFermionD Ddwf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,b,c);

  SchurDiagMooeeOperator<MobiusFermionD, LatticeFermion> HermOpEO(Ddwf);
  
  const int Fine_Nstop = 200;
  const int Fine_Nk = 100;
  const int Fine_Np = 100;
  const int Fine_Nm = Fine_Nk+Fine_Np;

  typedef LatticeFermion FermionField;
  std::vector<RealD>        Fine_eval;
  std::vector<FermionField> Fine_evec;

  LatticeFermionF conv_tmp(FrbGridF);
  Fine_eval.resize(Fine_Nstop);
  Fine_evec.resize(Fine_Nstop,FrbGrid);
  std::string evec_file("/lustre/orion/phy157/proj-shared/phy157_dwf/paboyle/Subspace.phys48.evecF");
  LoadFineEvecs(Fine_evec,evec_file,conv_tmp);
  
  typedef HermOpAdaptor<LatticeFermionD> HermFineMatrix;
  HermFineMatrix FineHermOp(HermOpEO);

  ////////////////////////////////////////////////////////////
  ///////////// Coarse basis and Little Dirac Operator ///////
  ////////////////////////////////////////////////////////////
  typedef GeneralCoarsenedMatrix<vSpinColourVector,vTComplex,nbasis> LittleDiracOperator;
  typedef LittleDiracOperator::CoarseVector CoarseVector;

  NextToNextToNextToNearestStencilGeometry5D geom(Coarse5d);

  typedef Aggregation<vSpinColourVector,vTComplex,nbasis> Subspace;
  Subspace Aggregates(Coarse5d,FrbGrid,cb);

  ////////////////////////////////////////////////////////////
  // Need to check about red-black grid coarsening
  ////////////////////////////////////////////////////////////
  //  std::string subspace_file("/lustre/orion/phy157/proj-shared/phy157_dwf/paboyle/Subspace.phys48.mixed.2500.60");
  //  //  std::string subspace_file("/lustre/orion/phy157/proj-shared/phy157_dwf/paboyle/Subspace.phys48.new.62");
  //  std::string refine_file("/lustre/orion/phy157/proj-shared/phy157_dwf/paboyle/Subspace.phys48.evec");
  std::string refine_file("/lustre/orion/phy157/proj-shared/phy157_dwf/paboyle/Refine.phys48.mixed.2500.60");
  //  std::string ldop_file("/lustre/orion/phy157/proj-shared/phy157_dwf/paboyle/LittleDiracOp.phys48.mixed.60");
  //  std::string evec_file("/lustre/orion/phy157/proj-shared/phy157_dwf/paboyle/evecs.scidac");
  //  std::string eval_file("/lustre/orion/phy157/proj-shared/phy157_dwf/paboyle/eval.xml");
  bool load_agg=true;
  bool load_refine=true;

  //////////////////////////////////////////
  // Block projector for coarse/fine
  //////////////////////////////////////////
  MultiRHSBlockProject<LatticeFermionD> MrhsProjector;


  /////////////////////////////////////////////////
  // Mirs smoother
  /////////////////////////////////////////////////
  int ord=8;
  RealD lo=2.0;
  RealD MirsShift = lo;
  ShiftedHermOpLinearOperator<LatticeFermionD> ShiftedFineHermOp(HermOpEO,MirsShift);
  CGSmoother<LatticeFermionD> CGsmooth(ord,ShiftedFineHermOp) ;
  
  LoadBasis(Aggregates,refine_file);
  Aggregates.Orthogonalise();

  std::cout << "**************************************"<<std::endl;
  std::cout << " Using filtered subspace"<<std::endl;
  std::cout << "**************************************"<<std::endl;
  MrhsProjector.Allocate(nbasis,FrbGrid,Coarse5d);
  MrhsProjector.ImportBasis(Aggregates.subspace);

  FermionField Ftmp(FrbGrid);
  std::vector<FermionField> Fine_ev(1,FrbGrid);
  std::vector<FermionField> Fine_ev_compressed(1,FrbGrid);
  std::vector<CoarseVector>  c_evec(1,Coarse5d);
  for(int ev=0;ev<Fine_evec.size();ev++){
    Fine_ev[0] = Fine_evec[ev];
    MrhsProjector.blockProject(Fine_ev,c_evec);
    MrhsProjector.blockPromote(Fine_ev_compressed,c_evec);
    Ftmp = Fine_ev_compressed[0];
    RealD div = 1.0/ sqrt(norm2(Ftmp));
    Ftmp = Ftmp * div;
    std::cout << GridLogMessage<<" "<<ev<<" uncomp "<< norm2(Fine_ev[0])  <<std::endl;
    std::cout << GridLogMessage<<" "<<ev<<" comp   "<< norm2(Ftmp)  <<std::endl;
    Ftmp = Fine_ev[0] - Ftmp;
    std::cout << GridLogMessage<<" "<<ev<<" diff "<< norm2(Ftmp)  <<std::endl;
    CGsmooth(Fine_ev_compressed[0],Ftmp);
    Ftmp = Ftmp *lo;
    std::cout << GridLogMessage<<" "<<ev<<" smoothed "<< norm2(Ftmp)  <<std::endl;
    div = 1.0/ sqrt(norm2(Ftmp));
    Ftmp=Ftmp*div;
    Ftmp = Fine_ev[0]-Ftmp;
    std::cout << GridLogMessage<<" "<<ev<<" diff "<< norm2(Ftmp)  <<std::endl;
  }

  std::cout << "**************************************"<<std::endl;
  std::cout << " Using eigenvector subspace "<<std::endl;
  std::cout << "**************************************"<<std::endl;
  for(int i=0;i<Aggregates.subspace.size();i++){
    Aggregates.subspace[i] = Fine_evec[i];
  }
  Aggregates.Orthogonalise();
  MrhsProjector.ImportBasis(Aggregates.subspace);
  for(int ev=0;ev<Fine_evec.size();ev++){
    Fine_ev[0] = Fine_evec[ev];
    MrhsProjector.blockProject(Fine_ev,c_evec);
    MrhsProjector.blockPromote(Fine_ev_compressed,c_evec);
    Ftmp = Fine_ev_compressed[0];
    RealD div = 1.0/ sqrt(norm2(Ftmp));
    Ftmp = Ftmp * div;
    std::cout << GridLogMessage<<" "<<ev<<" uncomp "<< norm2(Fine_ev[0])  <<std::endl;
    std::cout << GridLogMessage<<" "<<ev<<" comp   "<< norm2(Ftmp)  <<std::endl;
    Ftmp = Fine_ev[0] - Ftmp;
    std::cout << GridLogMessage<<" "<<ev<<" diff "<< norm2(Ftmp)  <<std::endl;
    CGsmooth(Fine_ev_compressed[0],Ftmp);
    Ftmp = Ftmp *lo;
    std::cout << GridLogMessage<<" "<<ev<<" smoothed "<< norm2(Ftmp)  <<std::endl;
    div = 1.0/ sqrt(norm2(Ftmp));
    Ftmp=Ftmp*div;
    Ftmp = Fine_ev[0]-Ftmp;
    std::cout << GridLogMessage<<" "<<ev<<" diff "<< norm2(Ftmp)  <<std::endl;
  }

  // Standard CG
  Grid_finalize();
  return 0;
}
