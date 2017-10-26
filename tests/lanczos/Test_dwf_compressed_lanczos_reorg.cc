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

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

struct LanczosParams : Serializable {
 public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(LanczosParams,
				  ChebyParams, Cheby,/*Chebyshev*/
				  int, Nstop,    /*Vecs in Lanczos must converge Nstop < Nk < Nm*/
				  int, Nk,       /*Vecs in Lanczos seek converge*/
				  int, Nm,       /*Total vecs in Lanczos include restart*/
				  RealD, resid,  /*residual*/
 				  int, MaxIt, 
				  RealD, betastp,  /* ? */
				  int, MinRes);    // Must restart
};

struct CompressedLanczosParams : Serializable {
 public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(CompressedLanczosParams,
				  LanczosParams, FineParams,
				  LanczosParams, CoarseParams,
				  ChebyParams,   Smoother,
				  RealD        , coarse_relax_tol,
				  std::vector<int>, blockSize,
				  std::string, config,
				  std::vector < std::complex<double>  >, omega,
				  RealD, mass,
				  RealD, M5
				  );
};

// Duplicate functionality; ProjectedFunctionHermOp could be used with the trivial function
template<class Fobj,class CComplex,int nbasis>
class ProjectedHermOp : public LinearFunction<Lattice<iVector<CComplex,nbasis > > > {
public:
  typedef iVector<CComplex,nbasis >           CoarseSiteVector;
  typedef Lattice<CoarseSiteVector>           CoarseField;
  typedef Lattice<CComplex>   CoarseScalar; // used for inner products on fine field
  typedef Lattice<Fobj>          FineField;

  LinearOperatorBase<FineField> &_Linop;
  Aggregation<Fobj,CComplex,nbasis> &_Aggregate;

  ProjectedHermOp(LinearOperatorBase<FineField>& linop,  Aggregation<Fobj,CComplex,nbasis> &aggregate) : 
    _Linop(linop),
    _Aggregate(aggregate)  {  };

  void operator()(const CoarseField& in, CoarseField& out) {

    GridBase *FineGrid = _Aggregate.FineGrid;
    FineField fin(FineGrid);
    FineField fout(FineGrid);

    _Aggregate.PromoteFromSubspace(in,fin);    std::cout<<GridLogIRL<<"ProjectedHermop : Promote to fine"<<std::endl;
    _Linop.HermOp(fin,fout);                   std::cout<<GridLogIRL<<"ProjectedHermop : HermOp (fine) "<<std::endl;
    _Aggregate.ProjectToSubspace(out,fout);    std::cout<<GridLogIRL<<"ProjectedHermop : Project to coarse "<<std::endl;
  }
};

template<class Fobj,class CComplex,int nbasis>
class ProjectedFunctionHermOp : public LinearFunction<Lattice<iVector<CComplex,nbasis > > > {
public:
  typedef iVector<CComplex,nbasis >           CoarseSiteVector;
  typedef Lattice<CoarseSiteVector>           CoarseField;
  typedef Lattice<CComplex>   CoarseScalar; // used for inner products on fine field
  typedef Lattice<Fobj>          FineField;


  OperatorFunction<FineField>   & _poly;
  LinearOperatorBase<FineField> &_Linop;
  Aggregation<Fobj,CComplex,nbasis> &_Aggregate;

  ProjectedFunctionHermOp(OperatorFunction<FineField> & poly,LinearOperatorBase<FineField>& linop, 
			  Aggregation<Fobj,CComplex,nbasis> &aggregate) : 
    _poly(poly),
    _Linop(linop),
    _Aggregate(aggregate)  {  };

  void operator()(const CoarseField& in, CoarseField& out) {

    GridBase *FineGrid = _Aggregate.FineGrid;

    FineField fin(FineGrid) ;fin.checkerboard  =_Aggregate.checkerboard;
    FineField fout(FineGrid);fout.checkerboard =_Aggregate.checkerboard;
    
    _Aggregate.PromoteFromSubspace(in,fin);    std::cout<<GridLogIRL<<"ProjectedFunctionHermop : Promote to fine"<<std::endl;
    _poly(_Linop,fin,fout);                    std::cout<<GridLogIRL<<"ProjectedFunctionHermop : Poly "<<std::endl;
    _Aggregate.ProjectToSubspace(out,fout);    std::cout<<GridLogIRL<<"ProjectedFunctionHermop : Project to coarse "<<std::endl;
  }
};

template<class Fobj,class CComplex,int nbasis>
class ImplicitlyRestartedLanczosSmoothedTester  : public ImplicitlyRestartedLanczosTester<Lattice<iVector<CComplex,nbasis > > >
{
 public:
  typedef iVector<CComplex,nbasis >           CoarseSiteVector;
  typedef Lattice<CoarseSiteVector>           CoarseField;
  typedef Lattice<CComplex>   CoarseScalar; // used for inner products on fine field
  typedef Lattice<Fobj>          FineField;

  LinearFunction<CoarseField> & _Poly;
  OperatorFunction<FineField>   & _smoother;
  LinearOperatorBase<FineField> &_Linop;
  Aggregation<Fobj,CComplex,nbasis> &_Aggregate;
  RealD                             _coarse_relax_tol;
  ImplicitlyRestartedLanczosSmoothedTester(LinearFunction<CoarseField>   &Poly,
					   OperatorFunction<FineField>   &smoother,
					   LinearOperatorBase<FineField> &Linop,
					   Aggregation<Fobj,CComplex,nbasis> &Aggregate,
					   RealD coarse_relax_tol=5.0e3) 
    : _smoother(smoother), _Linop(Linop),_Aggregate(Aggregate), _Poly(Poly), _coarse_relax_tol(coarse_relax_tol)  {    };

  int TestConvergence(int j,RealD eresid,CoarseField &B, RealD &eval,RealD evalMaxApprox)
  {
    CoarseField v(B);
    RealD eval_poly = eval;
    // Apply operator
    _Poly(B,v);

    RealD vnum = real(innerProduct(B,v)); // HermOp.
    RealD vden = norm2(B);
    RealD vv0  = norm2(v);
    eval   = vnum/vden;
    v -= eval*B;

    RealD vv = norm2(v) / ::pow(evalMaxApprox,2.0);

    std::cout.precision(13);
    std::cout<<GridLogIRL  << "[" << std::setw(3)<<j<<"] "
	     <<"eval = "<<std::setw(25)<< eval << " (" << eval_poly << ")"
	     <<" |H B[i] - eval[i]B[i]|^2 / evalMaxApprox^2 " << std::setw(25) << vv
	     <<std::endl;

    int conv=0;
    if( (vv<eresid*eresid) ) conv = 1;
    return conv;
  }
  int ReconstructEval(int j,RealD eresid,CoarseField &B, RealD &eval,RealD evalMaxApprox)
  {
    GridBase *FineGrid = _Aggregate.FineGrid;

    int checkerboard   = _Aggregate.checkerboard;

    FineField fB(FineGrid);fB.checkerboard =checkerboard;
    FineField fv(FineGrid);fv.checkerboard =checkerboard;

    _Aggregate.PromoteFromSubspace(B,fv);
    _smoother(_Linop,fv,fB); 

    RealD eval_poly = eval;
    _Linop.HermOp(fB,fv);

    RealD vnum = real(innerProduct(fB,fv)); // HermOp.
    RealD vden = norm2(fB);
    RealD vv0  = norm2(fv);
    eval   = vnum/vden;
    fv -= eval*fB;
    RealD vv = norm2(fv) / ::pow(evalMaxApprox,2.0);

    std::cout.precision(13);
    std::cout<<GridLogIRL  << "[" << std::setw(3)<<j<<"] "
	     <<"eval = "<<std::setw(25)<< eval << " (" << eval_poly << ")"
	     <<" |H B[i] - eval[i]B[i]|^2 / evalMaxApprox^2 " << std::setw(25) << vv
	     <<std::endl;
    if ( j > nbasis ) eresid = eresid*_coarse_relax_tol;
    if( (vv<eresid*eresid) ) return 1;
    return 0;
  }
};


////////////////////////////////////////////
// Make serializable Lanczos params
////////////////////////////////////////////
template<class Fobj,class CComplex,int nbasis>
class CoarseFineIRL 
{
public:
  typedef iVector<CComplex,nbasis >           CoarseSiteVector;
  typedef Lattice<CComplex>                   CoarseScalar; // used for inner products on fine field
  typedef Lattice<CoarseSiteVector>           CoarseField;
  typedef Lattice<Fobj>                       FineField;

private:
  GridBase *_CoarseGrid;
  GridBase *_FineGrid;
  int _checkerboard;
  LinearOperatorBase<FineField>                 & _FineOp;
  
  // FIXME replace Aggregation with vector of fine; the code reuse is too small for
  // the hassle and complexity of cross coupling.
  Aggregation<Fobj,CComplex,nbasis>               _Aggregate;  
  std::vector<RealD>                              evals_fine;
  std::vector<RealD>                              evals_coarse; 
  std::vector<CoarseField>                        evec_coarse;
public:
  CoarseFineIRL(GridBase *FineGrid,
		GridBase *CoarseGrid,
		LinearOperatorBase<FineField> &FineOp,
		int checkerboard) :
    _CoarseGrid(CoarseGrid),
    _FineGrid(FineGrid),
    _Aggregate(CoarseGrid,FineGrid,checkerboard),
    _FineOp(FineOp),
    _checkerboard(checkerboard)
  {
    evals_fine.resize(0);
    evals_coarse.resize(0);
  };
  void Orthogonalise(void ) { _Aggregate.Orthogonalise(); }

  template<typename T>  static RealD normalise(T& v) 
  {
    RealD nn = norm2(v);
    nn = ::sqrt(nn);
    v = v * (1.0/nn);
    return nn;
  }

  void testFine(void)
  {
    int Nk = nbasis;
    _Aggregate.subspace.resize(Nk,_FineGrid);
    _Aggregate.subspace[0]=1.0;
    _Aggregate.subspace[0].checkerboard=_checkerboard;
    normalise(_Aggregate.subspace[0]);
    PlainHermOp<FineField>    Op(_FineOp);
    for(int k=1;k<Nk;k++){
      _Aggregate.subspace[k].checkerboard=_checkerboard;
      Op(_Aggregate.subspace[k-1],_Aggregate.subspace[k]);
      normalise(_Aggregate.subspace[k]);
    }
  }


  void checkpointFine(std::string evecs_file,std::string evals_file)
  {
    assert(_Aggregate.subspace.size()==nbasis);
    emptyUserRecord record;
    {
      ScidacWriter WR;
      WR.open(evecs_file);
      for(int k=0;k<nbasis;k++) {
	WR.writeScidacFieldRecord(_Aggregate.subspace[k],record);
      }
      WR.close();
    }
    {
      XmlWriter WR(evals_file);
      write(WR,"evals",evals_fine);
    }
  }
  void checkpointCoarse(std::string evecs_file,std::string evals_file)
  {
    int n = evec_coarse.size();
    emptyUserRecord record;
    {
      ScidacWriter WR;
      WR.open(evecs_file);
      for(int k=0;k<n;k++) {
	WR.writeScidacFieldRecord(evec_coarse[k],record);
      }
      WR.close();
    }
    {
      XmlWriter WR(evals_file);
      write(WR,"evals",evals_coarse);
    }
  }

  void checkpointFineRestore(std::string evecs_file,std::string evals_file)
  {
    {
      XmlReader RD(evals_file);
      read(RD,"evals",evals_fine);
    }
    assert(evals_fine.size()==nbasis);

    emptyUserRecord record;
    {
      ScidacReader RD ;
      RD.open(evecs_file);
      for(int k=0;k<nbasis;k++) {
	RD.readScidacFieldRecord(_Aggregate.subspace[k],record);
      }
      RD.close();
    }
  }

  void calcFine(ChebyParams cheby_parms,int Nstop,int Nk,int Nm,RealD resid, 
		RealD MaxIt, RealD betastp, int MinRes)
  {
    assert(nbasis<=Nm);
    Chebyshev<FineField>      Cheby(cheby_parms);
    FunctionHermOp<FineField> ChebyOp(Cheby,_FineOp);
    PlainHermOp<FineField>    Op(_FineOp);

    evals_fine.resize(Nm);
    _Aggregate.subspace.resize(Nm,_FineGrid);

    ImplicitlyRestartedLanczos<FineField> IRL(ChebyOp,Op,Nstop,Nk,Nm,resid,MaxIt,betastp,MinRes);

    FineField src(_FineGrid); src=1.0; src.checkerboard = _checkerboard;

    int Nconv;
    IRL.calc(evals_fine,_Aggregate.subspace,src,Nconv,false);
    
    // Shrink down to number saved
    assert(Nstop>=nbasis);
    assert(Nconv>=nbasis);
    evals_fine.resize(nbasis);
    _Aggregate.subspace.resize(nbasis,_FineGrid);
  }
  void calcCoarse(ChebyParams cheby_op,ChebyParams cheby_smooth,RealD relax,
		  int Nstop, int Nk, int Nm,RealD resid, 
		  RealD MaxIt, RealD betastp, int MinRes)
  {
    Chebyshev<FineField>                          Cheby(cheby_op);
    ProjectedHermOp<Fobj,CComplex,nbasis>         Op(_FineOp,_Aggregate);
    ProjectedFunctionHermOp<Fobj,CComplex,nbasis> ChebyOp (Cheby,_FineOp,_Aggregate);
    //////////////////////////////////////////////////////////////////////////////////////////////////
    // create a smoother and see if we can get a cheap convergence test and smooth inside the IRL
    //////////////////////////////////////////////////////////////////////////////////////////////////

    Chebyshev<FineField>                                           ChebySmooth(cheby_smooth);
    ImplicitlyRestartedLanczosSmoothedTester<Fobj,CComplex,nbasis> ChebySmoothTester(ChebyOp,ChebySmooth,_FineOp,_Aggregate,relax);

    evals_coarse.resize(Nm);
    evec_coarse.resize(Nm,_CoarseGrid);

    CoarseField src(_CoarseGrid);     src=1.0; 

    ImplicitlyRestartedLanczos<CoarseField> IRL(ChebyOp,ChebyOp,ChebySmoothTester,Nstop,Nk,Nm,resid,MaxIt,betastp,MinRes);
    int Nconv=0;
    IRL.calc(evals_coarse,evec_coarse,src,Nconv,false);
    assert(Nconv>=Nstop);

    for (int i=0;i<Nstop;i++){
      std::cout << i << " Coarse eval = " << evals_coarse[i]  << std::endl;
    }
  }
};


int main (int argc, char ** argv) {

  Grid_init(&argc,&argv);
  GridLogIRL.TimingMode(1);

  CompressedLanczosParams Params;
  {
    Params.omega.resize(10);
    Params.blockSize.resize(5);
    XmlWriter writer("Params_template.xml");
    write(writer,"Params",Params);
    std::cout << GridLogMessage << " Written Params_template.xml" <<std::endl;
  }
  
  { 
    XmlReader reader("./Params.xml");
    read(reader, "Params", Params);
  }

  int     Ls = (int)Params.omega.size();
  RealD mass = Params.mass;
  RealD M5   = Params.M5;
  std::vector<int> blockSize = Params.blockSize;

  // Grids
  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid   = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid     = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid   = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  std::vector<int> fineLatt     = GridDefaultLatt();
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
  GridRedBlackCartesian * CoarseGrid5rb  = SpaceTimeGrid::makeFourDimRedBlackGrid(CoarseGrid5);

  // Gauge field
  LatticeGaugeField Umu(UGrid);
  FieldMetaData header;
  NerscIO::readConfiguration(Umu,header,Params.config);
  std::cout << GridLogMessage << "Lattice dimensions: " << GridDefaultLatt() << "   Ls: " << Ls << std::endl;

  // ZMobius EO Operator
  ZMobiusFermionR Ddwf(Umu, *FGrid, *FrbGrid, *UGrid, *UrbGrid, mass, M5, Params.omega,1.,0.);
  SchurDiagTwoOperator<ZMobiusFermionR,LatticeFermion> HermOp(Ddwf);

  // Eigenvector storage
  LanczosParams fine  =Params.FineParams;  
  LanczosParams coarse=Params.CoarseParams;  

  const int Ns1 = fine.Nstop;   const int Ns2 = coarse.Nstop;
  const int Nk1 = fine.Nk;      const int Nk2 = coarse.Nk;
  const int Nm1 = fine.Nm;      const int Nm2 = coarse.Nm;

  std::cout << GridLogMessage << "Keep " << fine.Nstop   << " fine   vectors" << std::endl;
  std::cout << GridLogMessage << "Keep " << coarse.Nstop << " coarse vectors" << std::endl;
  assert(Nm2 >= Nm1);

  const int nbasis= 60;
  assert(nbasis==Ns1);
  CoarseFineIRL<vSpinColourVector,vTComplex,nbasis> IRL(FrbGrid,CoarseGrid5rb,HermOp,Odd);
  std::cout << GridLogMessage << "Constructed CoarseFine IRL" << std::endl;

  int do_fine   = 1;
  int do_coarse = 0;
  int do_smooth = 0;
  if ( do_fine ) { 
    std::cout << GridLogMessage << "Performing fine grid IRL Nstop "<< Ns1 << " Nk "<<Nk1<<" Nm "<<Nm1<< std::endl;
    IRL.calcFine(fine.Cheby,
		 fine.Nstop,fine.Nk,fine.Nm,
		 fine.resid,fine.MaxIt, 
		 fine.betastp,fine.MinRes);

    std::cout << GridLogIRL<<"Checkpointing Fine evecs"<<std::endl;
    IRL.checkpointFine(std::string("evecs.scidac"),std::string("evals.xml"));
  } else { 
    //    IRL.testFine();
    IRL.checkpointFineRestore(std::string("evecs.scidac"),std::string("evals.xml"));
  }
  
  std::cout << GridLogMessage << "Orthogonalising " << nbasis<<" Nm "<<Nm2<< std::endl;
  IRL.Orthogonalise();

  std::cout << GridLogMessage << "Performing coarse grid IRL Nstop "<< Ns2<< " Nk "<<Nk2<<" Nm "<<Nm2<< std::endl;
  IRL.calcCoarse(coarse.Cheby,Params.Smoother,Params.coarse_relax_tol,
		 coarse.Nstop, coarse.Nk,coarse.Nm,
		 coarse.resid, coarse.MaxIt, 
		 coarse.betastp,coarse.MinRes);


  std::cout << GridLogIRL<<"Checkpointing coarse evecs"<<std::endl;
  IRL.checkpointCoarse(std::string("evecs.coarse.scidac"),std::string("evals.coarse.xml"));

  Grid_finalize();
}

