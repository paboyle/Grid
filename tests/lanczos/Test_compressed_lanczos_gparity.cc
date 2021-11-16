    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_compressed_lanczos_gparity.cc

    Copyright (C) 2017

Author: Christopher Kelly <ckelly@bnl.gov>
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

// template<class VectorInt>
// void GridCmdOptionIntVector(const std::string &str, VectorInt & vec)
// {
//   vec.resize(0);
//   std::stringstream ss(str);
//   int i;
//   while (ss >> i){
//     vec.push_back(i);
//     if(std::ispunct(ss.peek()))
//       ss.ignore();
//   }
//   return;
// }


//For the CPS configurations we have to manually seed the RNG and deal with an incorrect factor of 2 in the plaquette metadata
void readConfiguration(LatticeGaugeFieldD &U,
		       const std::string &config,
		       bool is_cps_cfg = false){

  if(is_cps_cfg) NerscIO::exitOnReadPlaquetteMismatch() = false;

  typedef GaugeStatistics<ConjugateGimplD> GaugeStats;
     
  FieldMetaData header;
  NerscIO::readConfiguration<GaugeStats>(U, header, config);

  if(is_cps_cfg) NerscIO::exitOnReadPlaquetteMismatch() = true;
}

//Lanczos parameters in CPS conventions
struct CPSLanczosParams : Serializable {
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(CPSLanczosParams,
				  RealD, alpha,
				  RealD, beta,
				  int, ch_ord,
				  int, N_use,
				  int, N_get,
				  int, N_true_get,
				  RealD, stop_rsd,
				  int, maxits);

  //Translations
  ChebyParams getChebyParams() const{
    ChebyParams out;
    out.alpha = beta*beta; //aka lo
    out.beta = alpha*alpha; //aka hi
    out.Npoly = ch_ord+1;
    return out;
  }
  int Nstop() const{ return N_true_get; }
  int Nm() const{ return N_use; }
  int Nk() const{ return N_get; }
};

//Maybe this class should be in the main library?
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
  }

  void checkpointFineRestore(std::string evecs_file,std::string evals_file)
  {
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
      this->subspace[k].Checkerboard()=this->_checkerboard;
      RD.readScidacFieldRecord(this->subspace[k],record);
      
    }
    RD.close();
  }

  void checkpointCoarse(std::string evecs_file,std::string evals_file)
  {
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
  }

  void checkpointCoarseRestore(std::string evecs_file,std::string evals_file,int nvec)
  {
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
  }
};


//Note:  because we rely upon physical properties we must use a "real" gauge configuration
int main (int argc, char ** argv) {
  Grid_init(&argc,&argv);
  GridLogIRL.TimingMode(1);

  std::vector<int> blockSize = {2,2,2,2,2};
  std::vector<int> GparityDirs = {1,1,1}; //1 for each GP direction

  int Ls = 12;
  RealD mass = 0.01;
  RealD M5 = 1.8;
  bool is_cps_cfg = false;

  CPSLanczosParams fine, coarse;

  fine.alpha = 2;
  fine.beta = 0.1;
  fine.ch_ord = 100;
  fine.N_use = 70;
  fine.N_get = 60;
  fine.N_true_get = 60;
  fine.stop_rsd = 1e-8;
  fine.maxits = 10000;

  coarse.alpha = 2;
  coarse.beta = 0.1;
  coarse.ch_ord = 100;
  coarse.N_use = 200;
  coarse.N_get = 190;
  coarse.N_true_get = 190;
  coarse.stop_rsd = 1e-8;
  coarse.maxits = 10000;

  double coarse_relax_tol = 1e5;
  int smoother_ord = 20;
  
  if(argc < 3){
    std::cout << GridLogMessage << "Usage: <exe> <config> <gparity dirs> <options>" << std::endl;
    std::cout << GridLogMessage << "<gparity dirs> should have the format a.b.c where a,b,c are 0,1 depending on whether there are G-parity BCs in that direction" << std::endl;
    std::cout << GridLogMessage << "Options:" << std::endl;
    std::cout << GridLogMessage << "--Ls <value> : Set Ls (default 12)" << std::endl;
    std::cout << GridLogMessage << "--mass <value> : Set the mass (default 0.01)" << std::endl;
    std::cout << GridLogMessage << "--block <value> : Set the block size. Format should be a.b.c.d.e where a-e are the block extents  (default 2.2.2.2.2)" << std::endl;
    std::cout << GridLogMessage << "--is_cps_cfg : Indicate that the configuration was generated with CPS where until recently the stored plaquette was wrong by a factor of 2" << std::endl;
    std::cout << GridLogMessage << "--write_irl_templ: Write a template for the parameters file of the Lanczos to \"irl_templ.xml\"" << std::endl;
    std::cout << GridLogMessage << "--read_irl_fine <filename>: Real the parameters file for the fine Lanczos" << std::endl;
    std::cout << GridLogMessage << "--read_irl_coarse <filename>: Real the parameters file for the coarse Lanczos" << std::endl;
    std::cout << GridLogMessage << "--write_fine <filename stub>: Write fine evecs/evals to filename starting with the stub" << std::endl;
    std::cout << GridLogMessage << "--read_fine <filename stub>: Read fine evecs/evals from filename starting with the stub" << std::endl;
    std::cout << GridLogMessage << "--write_coarse <filename stub>: Write coarse evecs/evals to filename starting with the stub" << std::endl;
    std::cout << GridLogMessage << "--read_coarse <filename stub>: Read coarse evecs/evals from filename starting with the stub" << std::endl;
    std::cout << GridLogMessage << "--smoother_ord :  Set the Chebyshev order of the smoother (default 20)" << std::endl;
    std::cout << GridLogMessage << "--coarse_relax_tol : Set the relaxation parameter for evaluating the residual of the reconstructed eigenvectors outside of the basis (default 1e5)" << std::endl;
    Grid_finalize();
    return 1;
  }
  std::string config = argv[1];
  GridCmdOptionIntVector(argv[2], GparityDirs);
  assert(GparityDirs.size() == 3);

  bool write_fine = false;
  std::string write_fine_file;

  bool read_fine = false;
  std::string read_fine_file;

  bool write_coarse = false;
  std::string write_coarse_file;

  bool read_coarse = false;
  std::string read_coarse_file;
 
  for(int i=3;i<argc;i++){
    std::string sarg = argv[i];
    if(sarg == "--Ls"){
      Ls = std::stoi(argv[i+1]);
      std::cout << GridLogMessage << "Set Ls to " << Ls << std::endl;
    }else if(sarg == "--mass"){
      std::istringstream ss(argv[i+1]); ss >> mass;
      std::cout << GridLogMessage << "Set quark mass to " << mass << std::endl;
    }else if(sarg == "--block"){
      GridCmdOptionIntVector(argv[i+1], blockSize);
      assert(blockSize.size() == 5);
      std::cout << GridLogMessage << "Set block size to ";
      for(int q=0;q<5;q++) std::cout << blockSize[q] << " ";
      std::cout << std::endl;      
    }else if(sarg == "--is_cps_cfg"){
      is_cps_cfg = true;
    }else if(sarg == "--write_irl_templ"){
      XmlWriter writer("irl_templ.xml");
      write(writer,"Params",fine);
      Grid_finalize();
      return 0;
    }else if(sarg == "--read_irl_fine"){
      std::cout << GridLogMessage << "Reading fine IRL params from " << argv[i+1] << std::endl;
      XmlReader reader(argv[i+1]);
      read(reader, "Params", fine);
    }else if(sarg == "--read_irl_coarse"){
      std::cout << GridLogMessage << "Reading coarse IRL params from " << argv[i+1] << std::endl;
      XmlReader reader(argv[i+1]);
      read(reader, "Params", coarse);
    }else if(sarg == "--write_fine"){
      write_fine = true;
      write_fine_file = argv[i+1];
    }else if(sarg == "--read_fine"){
      read_fine = true;
      read_fine_file = argv[i+1];
    }else if(sarg == "--write_coarse"){
      write_coarse = true;
      write_coarse_file = argv[i+1];
    }else if(sarg == "--read_coarse"){
      read_coarse = true;
      read_coarse_file = argv[i+1];
    }else if(sarg == "--smoother_ord"){
      std::istringstream ss(argv[i+1]); ss >> smoother_ord;
      std::cout << GridLogMessage << "Set smoother order to " << smoother_ord << std::endl;
    }else if(sarg == "--coarse_relax_tol"){
      std::istringstream ss(argv[i+1]); ss >> coarse_relax_tol;
      std::cout << GridLogMessage << "Set coarse IRL relaxation parameter to " << coarse_relax_tol << std::endl;
    }      
  }
  
  //Fine grids
  GridCartesian         * UGrid     = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(),  GridDefaultSimd(Nd,vComplex::Nsimd()),   GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid   = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid     = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid   = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  //Setup G-parity BCs
  assert(Nd == 4);
  std::vector<int> dirs4(4);
  for(int i=0;i<3;i++) dirs4[i] = GparityDirs[i];
  dirs4[3] = 0; //periodic gauge BC in time
  
  std::cout << GridLogMessage << "Gauge BCs: " << dirs4 << std::endl;
  ConjugateGimplD::setDirections(dirs4); //gauge BC

  GparityWilsonImplD::ImplParams Params;
  for(int i=0;i<Nd-1;i++) Params.twists[i] = GparityDirs[i]; //G-parity directions
  Params.twists[Nd-1] = 1; //APBC in time direction
  std::cout << GridLogMessage << "Fermion BCs: " << Params.twists << std::endl;
  
  //Read the gauge field
  LatticeGaugeField Umu(UGrid);  
  readConfiguration(Umu, config, is_cps_cfg);

  //Setup the coarse grids  
  auto fineLatt     = GridDefaultLatt();
  Coordinate coarseLatt(4);
  for (int d=0;d<4;d++){
    coarseLatt[d] = fineLatt[d]/blockSize[d];    assert(coarseLatt[d]*blockSize[d]==fineLatt[d]);
  }

  std::cout << GridLogMessage<< " 5d coarse lattice is ";
  for (int i=0;i<4;i++){
    std::cout << coarseLatt[i]<<"x";
  } 
  int cLs = Ls/blockSize[4]; assert(cLs*blockSize[4]==Ls);
  std::cout << cLs<<std::endl;
  
  GridCartesian         * CoarseGrid4    = SpaceTimeGrid::makeFourDimGrid(coarseLatt, GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * CoarseGrid4rb  = SpaceTimeGrid::makeFourDimRedBlackGrid(CoarseGrid4);
  GridCartesian         * CoarseGrid5    = SpaceTimeGrid::makeFiveDimGrid(cLs,CoarseGrid4);

  //Dirac operator
  GparityDomainWallFermionD action(Umu, *FGrid, *FrbGrid, *UGrid, *UrbGrid, mass, M5, Params);
  typedef GparityDomainWallFermionD::FermionField FermionField;
  
  SchurDiagTwoOperator<GparityDomainWallFermionD,FermionField> SchurOp(action);

  typedef GparityWilsonImplD::SiteSpinor SiteSpinor;

  std::cout << GridLogMessage << "Keep " << fine.N_true_get   << " fine   vectors" << std::endl;
  std::cout << GridLogMessage << "Keep " << coarse.N_true_get << " coarse vectors" << std::endl;
  assert(coarse.N_true_get >= fine.N_true_get);

  const int nbasis= 60;
  assert(nbasis<=fine.N_true_get);
  LocalCoherenceLanczosScidac<SiteSpinor,vTComplex,nbasis> _LocalCoherenceLanczos(FrbGrid,CoarseGrid5,SchurOp,Odd);
  std::cout << GridLogMessage << "Constructed LocalCoherenceLanczos" << std::endl;

  //Compute and/or read fine evecs
  if(read_fine){
    _LocalCoherenceLanczos.checkpointFineRestore(read_fine_file + "_evecs.scidac", read_fine_file + "_evals.xml");
  }else{
    std::cout << GridLogMessage << "Performing fine grid IRL" << std::endl;
    std::cout << GridLogMessage << "Using Chebyshev alpha=" << fine.alpha << " beta=" << fine.beta << " ord=" << fine.ch_ord << std::endl;
    _LocalCoherenceLanczos.calcFine(fine.getChebyParams(),
				    fine.Nstop(),fine.Nk(),fine.Nm(),
				    fine.stop_rsd,fine.maxits,0,0);
    if(write_fine){
      std::cout << GridLogIRL<<"Checkpointing Fine evecs"<<std::endl;
      _LocalCoherenceLanczos.checkpointFine(write_fine_file + "_evecs.scidac", write_fine_file + "_evals.xml");
    }
  }
  
  //Block orthonormalise (this should be part of calcFine?)
  std::cout << GridLogIRL<<"Orthogonalising"<<std::endl;
  _LocalCoherenceLanczos.Orthogonalise();
  std::cout << GridLogIRL<<"Orthogonaled"<<std::endl;

  ChebyParams smoother = fine.getChebyParams();
  smoother.Npoly = smoother_ord+1;

  if(read_coarse){
    _LocalCoherenceLanczos.checkpointCoarseRestore(read_coarse_file + "_evecs.scidac", read_coarse_file + "_evals.xml",coarse.Nstop());

  }else{
    std::cout << GridLogMessage << "Performing coarse grid IRL" << std::endl;
    std::cout << GridLogMessage << "Using Chebyshev alpha=" << coarse.alpha << " beta=" << coarse.beta << " ord=" << coarse.ch_ord << std::endl;	
    _LocalCoherenceLanczos.calcCoarse(coarse.getChebyParams(), smoother, coarse_relax_tol,
				      coarse.Nstop(), coarse.Nk() ,coarse.Nm(),
				      coarse.stop_rsd, coarse.maxits, 
				      0,0);

    if(write_coarse){
      std::cout << GridLogIRL<<"Checkpointing Coarse evecs"<<std::endl;
      _LocalCoherenceLanczos.checkpointCoarse(write_coarse_file + "_evecs.scidac", write_coarse_file + "_evals.xml");
    }

  }

  //Test the eigenvectors
  FermionField evec(FrbGrid);
  FermionField tmp(FrbGrid);
  RealD eval;
  
  for(int i=0;i<coarse.N_true_get;i++){    
    _LocalCoherenceLanczos.getFineEvecEval(evec, eval, i);
    SchurOp.HermOp(evec, tmp);
    tmp = tmp - eval*evec;
    std::cout << GridLogMessage << "Eval " << eval << " resid " << sqrt(norm2(tmp)) << std::endl;
  }

  Grid_finalize();
}

