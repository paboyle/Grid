    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_dwf_compressed_lanczos_reorg.cc

    Copyright (C) 2017

Author: Leans heavily on Christoph Lehner's code
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
/*
 *  Reimplement the badly named "multigrid" lanczos as compressed Lanczos using the features 
 *  in Grid that were intended to be used to support blocked Aggregates, from
 */
#include <Grid/Grid.h>
#include <Grid/algorithms/iterative/ImplicitlyRestartedLanczos.h>
#include <Grid/algorithms/iterative/LocalCoherenceLanczos.h>

using namespace std;
using namespace Grid;

template<class Fobj,class CComplex,int nbasis>
class LocalCoherenceLanczosScidac : public LocalCoherenceLanczos<Fobj,CComplex,nbasis>
{ 
public:
  typedef iVector<CComplex,nbasis >           CoarseSiteVector;
  typedef Lattice<CoarseSiteVector>           CoarseField;
  typedef Lattice<CComplex>   CoarseScalar; // used for inner products on fine field
  typedef Lattice<Fobj>          FineField;

  LocalCoherenceLanczosScidac(GridBase *FineGrid,GridBase *CoarseGrid,
			      LinearOperatorBase<FineField> &FineOp,
			      int checkerboard) 
    // Base constructor
    : LocalCoherenceLanczos<Fobj,CComplex,nbasis>(FineGrid,CoarseGrid,FineOp,checkerboard) 
  {};

  void checkpointFine(std::string evecs_file,std::string evals_file)
  {
#ifdef HAVE_LIME
    assert(this->subspace.size()==nbasis);
    emptyUserRecord record;
    Grid::ScidacWriter WR(this->_FineGrid->IsBoss());
    WR.open(evecs_file);
    for(int k=0;k<nbasis;k++) {
      WR.writeScidacFieldRecord(this->subspace[k],record);
    }
    WR.close();
    
    XmlWriter WRx(evals_file);
    write(WRx,"evals",this->evals_fine);
#else
    assert(0);
#endif
  }

  void checkpointFineRestore(std::string evecs_file,std::string evals_file)
  {
#ifdef HAVE_LIME
    this->evals_fine.resize(nbasis);
    this->subspace.resize(nbasis,this->_FineGrid);
    
    std::cout << GridLogIRL<< "checkpointFineRestore:  Reading evals from "<<evals_file<<std::endl;
    XmlReader RDx(evals_file);
    read(RDx,"evals",this->evals_fine);
    
    assert(this->evals_fine.size()==nbasis);
    
    std::cout << GridLogIRL<< "checkpointFineRestore:  Reading evecs from "<<evecs_file<<std::endl;
    emptyUserRecord record;
    Grid::ScidacReader RD ;
    RD.open(evecs_file);
    for(int k=0;k<nbasis;k++) {
      this->subspace[k].checkerboard=this->_checkerboard;
      RD.readScidacFieldRecord(this->subspace[k],record);
      
    }
    RD.close();
#else
    assert(0);
#endif 
  }

  void checkpointCoarse(std::string evecs_file,std::string evals_file)
  {
#ifdef HAVE_LIME
    int n = this->evec_coarse.size();
    emptyUserRecord record;
    Grid::ScidacWriter WR(this->_CoarseGrid->IsBoss());
    WR.open(evecs_file);
    for(int k=0;k<n;k++) {
      WR.writeScidacFieldRecord(this->evec_coarse[k],record);
    }
    WR.close();
    
    XmlWriter WRx(evals_file);
    write(WRx,"evals",this->evals_coarse);
#else
    assert(0);
#endif
  }

  void checkpointCoarseRestore(std::string evecs_file,std::string evals_file,int nvec)
  {
#ifdef HAVE_LIME
    std::cout << "resizing coarse vecs to " << nvec<< std::endl;
    this->evals_coarse.resize(nvec);
    this->evec_coarse.resize(nvec,this->_CoarseGrid);
    std::cout << GridLogIRL<< "checkpointCoarseRestore:  Reading evals from "<<evals_file<<std::endl;
    XmlReader RDx(evals_file);
    read(RDx,"evals",this->evals_coarse);

    assert(this->evals_coarse.size()==nvec);
    emptyUserRecord record;
    std::cout << GridLogIRL<< "checkpointCoarseRestore:  Reading evecs from "<<evecs_file<<std::endl;
    Grid::ScidacReader RD ;
    RD.open(evecs_file);
    for(int k=0;k<nvec;k++) {
      RD.readScidacFieldRecord(this->evec_coarse[k],record);
    }
    RD.close();
#else 
    assert(0);
#endif
  }
};

int main (int argc, char ** argv) {

  Grid_init(&argc,&argv);
  GridLogIRL.TimingMode(1);

  LocalCoherenceLanczosParams Params;
  {
    Params.omega.resize(10);
    Params.blockSize.resize(5);
    XmlWriter writer("Params_template.xml");
    write(writer,"Params",Params);
    std::cout << GridLogMessage << " Written Params_template.xml" <<std::endl;
  }
  
  { 
    XmlReader reader(std::string("./Params.xml"));
    read(reader, "Params", Params);
  }

  int     Ls = (int)Params.omega.size();
  RealD mass = Params.mass;
  RealD M5   = Params.M5;
  std::vector<int> blockSize = Params.blockSize;
  std::vector<int> latt({32,32,32,32});
  uint64_t     vol = Ls*latt[0]*latt[1]*latt[2]*latt[3];
  //  double   mat_flop= 2.0*1320.0*vol;    
  // Grids
  GridCartesian         * UGrid     = SpaceTimeGrid::makeFourDimGrid(latt,
								     GridDefaultSimd(Nd,vComplex::Nsimd()),
								     GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid   = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid     = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid   = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  std::vector<int> fineLatt     = latt;
  int dims=fineLatt.size();
  assert(blockSize.size()==dims+1);
  std::vector<int> coarseLatt(dims);
  std::vector<int> coarseLatt5d ;

  for (int d=0;d<coarseLatt.size();d++){
    coarseLatt[d] = fineLatt[d]/blockSize[d];    assert(coarseLatt[d]*blockSize[d]==fineLatt[d]);
  }

  std::cout << GridLogMessage<< " 5d coarse lattice is ";
  for (int i=0;i<coarseLatt.size();i++){
    std::cout << coarseLatt[i]<<"x";
  } 
  int cLs = Ls/blockSize[dims]; assert(cLs*blockSize[dims]==Ls);
  std::cout << cLs<<std::endl;
  
  GridCartesian         * CoarseGrid4    = SpaceTimeGrid::makeFourDimGrid(coarseLatt, GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * CoarseGrid4rb  = SpaceTimeGrid::makeFourDimRedBlackGrid(CoarseGrid4);
  GridCartesian         * CoarseGrid5    = SpaceTimeGrid::makeFiveDimGrid(cLs,CoarseGrid4);

  // Gauge field
  std::vector<int> seeds4({1,2,3,4});
  GridParallelRNG          RNG4(UGrid);  RNG4.SeedFixedIntegers(seeds4);
  LatticeGaugeField Umu(UGrid);
  SU<Nc>::HotConfiguration(RNG4,Umu);
  //  FieldMetaData header;
  //  NerscIO::readConfiguration(Umu,header,Params.config);

  std::cout << GridLogMessage << "Lattice dimensions: " << latt << "   Ls: " << Ls << std::endl;

  // ZMobius EO Operator
  ZMobiusFermionR Ddwf(Umu, *FGrid, *FrbGrid, *UGrid, *UrbGrid, mass, M5, Params.omega,1.,0.);
  SchurDiagTwoOperator<ZMobiusFermionR,LatticeFermion> HermOp(Ddwf);

  // Eigenvector storage
  LanczosParams fine  =Params.FineParams;  
  LanczosParams coarse=Params.CoarseParams;  

  const int Ns1 = fine.Nstop;   //const int Ns2 = coarse.Nstop;
  const int Nk1 = fine.Nk;      //const int Nk2 = coarse.Nk;
  const int Nm1 = fine.Nm;      const int Nm2 = coarse.Nm;

  std::cout << GridLogMessage << "Keep " << fine.Nstop   << " fine   vectors" << std::endl;
  std::cout << GridLogMessage << "Keep " << coarse.Nstop << " coarse vectors" << std::endl;
  assert(Nm2 >= Nm1);

  const int nbasis= 60;
  assert(nbasis==Ns1);
  LocalCoherenceLanczosScidac<vSpinColourVector,vTComplex,nbasis> _LocalCoherenceLanczos(FrbGrid,CoarseGrid5,HermOp,Odd);
  std::cout << GridLogMessage << "Constructed LocalCoherenceLanczos" << std::endl;

  assert( (Params.doFine)||(Params.doFineRead));

  if ( Params.doFine ) { 
    std::cout << GridLogMessage << "Performing fine grid IRL Nstop "<< Ns1 << " Nk "<<Nk1<<" Nm "<<Nm1<< std::endl;
    double t0=-usecond();
    _LocalCoherenceLanczos.calcFine(fine.Cheby,
		 fine.Nstop,fine.Nk,fine.Nm,
		 fine.resid,fine.MaxIt, 
		 fine.betastp,fine.MinRes);
    t0+=usecond();

    double t1=-usecond();
    if ( Params.saveEvecs ) {
      std::cout << GridLogIRL<<"Checkpointing Fine evecs"<<std::endl;
      _LocalCoherenceLanczos.checkpointFine(std::string("evecs.scidac"),std::string("evals.xml"));
    }
    t1+=usecond();

    std::cout << GridLogMessage << "Computation time is " << (t0)/1.0e6 <<" seconds"<<std::endl;
    if ( Params.saveEvecs )  std::cout << GridLogMessage << "I/O         time is " << (t1)/1.0e6 <<" seconds"<<std::endl;
    std::cout << GridLogMessage << "Time to solution is " << (t1+t0)/1.0e6 <<" seconds"<<std::endl;
    std::cout << GridLogMessage << "Done"<<std::endl;
  }

  Grid_finalize();
}

