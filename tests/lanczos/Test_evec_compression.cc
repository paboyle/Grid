    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_evec_compression.cc

    Copyright (C) 2017

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
/*
 *
 * This test generates eigenvectors using the Lanczos algorithm then attempts to use local coherence compression
 * to express those vectors in terms of a basis formed from a subset. This test is useful for finding the optimal
 * blocking and basis size for performing a Local Coherence Lanczos
 */
#include <Grid/Grid.h>
#include <Grid/algorithms/iterative/ImplicitlyRestartedLanczos.h>
#include <Grid/algorithms/iterative/LocalCoherenceLanczos.h>

using namespace std;
using namespace Grid;

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


template<class Fobj,class CComplex,int nbasis>
class LocalCoherenceCompressor{
public:
  typedef iVector<CComplex,nbasis >           CoarseSiteVector;
  typedef Lattice<CComplex>                   CoarseScalar; // used for inner products on fine field
  typedef Lattice<CoarseSiteVector>           CoarseField;
  typedef Lattice<Fobj>                       FineField;
  
  void compress(std::vector<FineField> &basis,
		std::vector<CoarseField> &compressed_evecs,
		const std::vector<FineField> &evecs_in,
		GridBase *FineGrid,
		GridBase *CoarseGrid){
    int nevecs = evecs_in.size();
    assert(nevecs > nbasis);
    
    //Construct the basis
    basis.resize(nbasis, FineGrid);
    for(int b=0;b<nbasis;b++) basis[b] = evecs_in[b];

    //Block othornormalize basis
    CoarseScalar InnerProd(CoarseGrid);
    std::cout << GridLogMessage <<" Gramm-Schmidt pass 1"<<std::endl;
    blockOrthogonalise(InnerProd,basis);
    std::cout << GridLogMessage <<" Gramm-Schmidt pass 2"<<std::endl;
    blockOrthogonalise(InnerProd,basis);

    //The coarse grid representation is the field of vectors of block inner products
    std::cout << GridLogMessage << "Compressing eigevectors" << std::endl;
    compressed_evecs.resize(nevecs, CoarseGrid);
    for(int i=0;i<nevecs;i++) blockProject(compressed_evecs[i], evecs_in[i], basis);
    std::cout << GridLogMessage << "Compression complete" << std::endl;
  }

  void uncompress(FineField &evec, const int i, const std::vector<FineField> &basis, const std::vector<CoarseField> &compressed_evecs) const{
    blockPromote(compressed_evecs[i],evec,basis);  
  }

  //Test uncompressed eigenvectors of Linop.HermOp to precision 'base_tolerance' for i<nbasis and 'base_tolerance*relax' for i>=nbasis
  bool testCompression(LinearOperatorBase<FineField> &Linop,
		       const std::vector<FineField> &basis, const std::vector<CoarseField> &compressed_evecs, const std::vector<RealD> &evals,
		       const RealD base_tolerance, const RealD relax){
    GridBase* FineGrid = basis[0].Grid();
    GridBase* CoarseGrid = compressed_evecs[0].Grid();

    bool fail = false;
    FineField evec(FineGrid), Mevec(FineGrid);
    for(int i=0;i<compressed_evecs.size();i++){
      std::cout << GridLogMessage << "Uncompressing evec " << i << std::endl;
      uncompress(evec, i, basis, compressed_evecs);

      std::cout << GridLogMessage << "Computing residual for evec " << i << std::endl;
      std::cout << GridLogMessage << "Linop" << std::endl;
      Linop.HermOp(evec, Mevec);
      std::cout << GridLogMessage << "Linalg" << std::endl;
      Mevec = Mevec - evals[i]*evec;

      std::cout << GridLogMessage << "Resid" << std::endl;
      RealD tol = base_tolerance * (i<nbasis ? 1. : relax);
      RealD res = sqrt(norm2(Mevec));
      std::cout << GridLogMessage << "Evec idx " << i << " res " << res << " tol " << tol << std::endl;
      if(res > tol) fail = true;
    }
    return fail;
  }
};

template<class Fobj,class CComplex,int nbasis>
void compareBlockPromoteTimings(const std::vector<Lattice<Fobj> > &basis, const std::vector<Lattice<iVector<CComplex,nbasis > > > &compressed_evecs){
  typedef iVector<CComplex,nbasis >           CoarseSiteVector;
  typedef Lattice<CComplex>                   CoarseScalar; 
  typedef Lattice<CoarseSiteVector>           CoarseField;
  typedef Lattice<Fobj>                       FineField;

  GridStopWatch timer;
  
  GridBase* FineGrid = basis[0].Grid();
  GridBase* CoarseGrid = compressed_evecs[0].Grid();

  FineField v1(FineGrid), v2(FineGrid);

  //Start with a cold start
  for(int i=0;i<basis.size();i++){
    autoView( b_ , basis[i], CpuWrite);
  }
  for(int i=0;i<compressed_evecs.size();i++){
    autoView( b_ , compressed_evecs[i], CpuWrite);
  }
  {
    autoView( b_, v1, CpuWrite );
  }

  timer.Start();
  blockPromote(compressed_evecs[0],v1,basis);  
  timer.Stop();
  std::cout << GridLogMessage << "Time for cold blockPromote v1 " << timer.Elapsed() << std::endl;

  //Test to ensure it is actually doing a cold start by repeating
  for(int i=0;i<basis.size();i++){
    autoView( b_ , basis[i], CpuWrite);
  }
  for(int i=0;i<compressed_evecs.size();i++){
    autoView( b_ , compressed_evecs[i], CpuWrite);
  }
  {
    autoView( b_, v1, CpuWrite );
  }

  timer.Reset();
  timer.Start();
  blockPromote(compressed_evecs[0],v1,basis);  
  timer.Stop();
  std::cout << GridLogMessage << "Time for cold blockPromote v1 repeat (should be the same as above) " << timer.Elapsed() << std::endl;
}


  

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

  CPSLanczosParams fine;

  fine.alpha = 2;
  fine.beta = 0.1;
  fine.ch_ord = 100;
  fine.N_use = 70;
  fine.N_get = 60;
  fine.N_true_get = 60;
  fine.stop_rsd = 1e-8;
  fine.maxits = 10000;

  double coarse_relax_tol = 1e5;
 
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
    std::cout << GridLogMessage << "--write_fine <filename stub>: Write fine evecs/evals to filename starting with the stub" << std::endl;
    std::cout << GridLogMessage << "--read_fine <filename stub>: Read fine evecs/evals from filename starting with the stub" << std::endl;
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
    }else if(sarg == "--write_fine"){
      write_fine = true;
      write_fine_file = argv[i+1];
    }else if(sarg == "--read_fine"){
      read_fine = true;
      read_fine_file = argv[i+1];
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
  const int nbasis= 60;
  typedef vTComplex CComplex; 
  typedef iVector<CComplex,nbasis >           CoarseSiteVector;
  typedef Lattice<CComplex>                   CoarseScalar;
  typedef Lattice<CoarseSiteVector>           CoarseField;
 
  //Dirac operator
  GparityDomainWallFermionD action(Umu, *FGrid, *FrbGrid, *UGrid, *UrbGrid, mass, M5, Params);
  typedef GparityDomainWallFermionD::FermionField FermionField;
  
  SchurDiagTwoOperator<GparityDomainWallFermionD,FermionField> SchurOp(action);

  typedef GparityWilsonImplD::SiteSpinor SiteSpinor;

  //Do the fine Lanczos
  std::vector<RealD> evals;
  std::vector<FermionField> evecs;

  if(read_fine){
    evals.resize(fine.N_true_get);
    evecs.resize(fine.N_true_get, FrbGrid);

    std::string evals_file = read_fine_file + "_evals.xml";
    std::string evecs_file = read_fine_file + "_evecs.scidac";
    
    std::cout << GridLogIRL<< "Reading evals from "<<evals_file<<std::endl;
    XmlReader RDx(evals_file);
    read(RDx,"evals",evals);
    
    assert(evals.size()==fine.N_true_get);
    
    std::cout << GridLogIRL<< "Reading evecs from "<<evecs_file<<std::endl;
    emptyUserRecord record;
    Grid::ScidacReader RD ;
    RD.open(evecs_file);
    for(int k=0;k<fine.N_true_get;k++) {
      evecs[k].Checkerboard()=Odd;
      RD.readScidacFieldRecord(evecs[k],record);
      
    }
    RD.close();
  }else{ 
    int Nstop = fine.Nstop();
    int Nm = fine.Nm();
    int Nk = fine.Nk();
    RealD resid = fine.stop_rsd;
    int MaxIt = fine.maxits;
    
    assert(nbasis<=Nm);    
    Chebyshev<FermionField>      Cheby(fine.getChebyParams());
    FunctionHermOp<FermionField> ChebyOp(Cheby,SchurOp);
    PlainHermOp<FermionField>    Op(SchurOp);

    evals.resize(Nm);
    evecs.resize(Nm,FrbGrid);
    
    ImplicitlyRestartedLanczos<FermionField> IRL(ChebyOp,Op,Nstop,Nk,Nm,resid,MaxIt,0,0);

    FermionField src(FrbGrid); 
    typedef typename FermionField::scalar_type Scalar;
    src=Scalar(1.0); 
    src.Checkerboard() = Odd;

    int Nconv;
    IRL.calc(evals, evecs,src,Nconv,false);

    if(write_fine){
      std::string evals_file = write_fine_file + "_evals.xml";
      std::string evecs_file = write_fine_file + "_evecs.scidac";

      std::cout << GridLogIRL<< "Writing evecs to "<<evecs_file<<std::endl;

      emptyUserRecord record;
      Grid::ScidacWriter WR(FrbGrid->IsBoss());
      WR.open(evecs_file);
      for(int k=0;k<evecs.size();k++) {
	WR.writeScidacFieldRecord(evecs[k],record);
      }
      WR.close();

      std::cout << GridLogIRL<< "Writing evals to "<<evals_file<<std::endl;
      
      XmlWriter WRx(evals_file);
      write(WRx,"evals",evals);
    }    
  }
    
  //Do the compression
  LocalCoherenceCompressor<SiteSpinor,vTComplex,nbasis> compressor;
  std::vector<FermionField> basis(nbasis,FrbGrid);
  std::vector<CoarseField> compressed_evecs(evecs.size(),CoarseGrid5);
  
  compressor.compress(basis, compressed_evecs, evecs, FrbGrid, CoarseGrid5);

  compareBlockPromoteTimings(basis, compressed_evecs);


  
  //Test the result
  assert( compressor.testCompression(SchurOp, basis, compressed_evecs, evals, fine.stop_rsd, coarse_relax_tol) );   

  Grid_finalize();
}
