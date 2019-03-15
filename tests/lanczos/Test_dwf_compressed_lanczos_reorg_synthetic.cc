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

    _Aggregate.PromoteFromSubspace(in,fin);
    _Linop.HermOp(fin,fout);
    _Aggregate.ProjectToSubspace(out,fout);
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

    _Aggregate.PromoteFromSubspace(in,fin);
    _poly(_Linop,fin,fout);
    _Aggregate.ProjectToSubspace(out,fout);
  }
};

// Make serializable Lanczos params

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
  Aggregation<Fobj,CComplex,nbasis>               _Aggregate; 

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
  {};

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
      Op(_Aggregate.subspace[k-1],_Aggregate.subspace[k]);
      normalise(_Aggregate.subspace[k]);
      std::cout << GridLogMessage << "testFine subspace "<<k<<" " <<norm2(_Aggregate.subspace[k])<<std::endl;
    }
    for(int k=0;k<Nk;k++){
      std::cout << GridLogMessage << "testFine subspace "<<k<<"  cb " <<_Aggregate.subspace[k].checkerboard<<std::endl;
    }
    _Aggregate.Orthogonalise();
  }

  void calcFine(RealD alpha, RealD beta,int Npoly,int Nm,RealD resid, 
		RealD MaxIt, RealD betastp, int MinRes)
  {
    assert(nbasis<=Nm);
    Chebyshev<FineField>      Cheby(alpha,beta,Npoly);
    FunctionHermOp<FineField> ChebyOp(Cheby,_FineOp);
    PlainHermOp<FineField>    Op(_FineOp);

    int Nk = nbasis;

    std::vector<RealD>          eval(Nm);

    FineField src(_FineGrid); src=1.0; src.checkerboard = _checkerboard;

    ImplicitlyRestartedLanczos<FineField> IRL(ChebyOp,Op,Nk,Nk,Nm,resid,MaxIt,betastp,MinRes);
    _Aggregate.subspace.resize(Nm,_FineGrid);
    IRL.calc(eval,_Aggregate.subspace,src,Nk,false);
    _Aggregate.subspace.resize(Nk,_FineGrid);
    for(int k=0;k<Nk;k++){
      std::cout << GridLogMessage << "testFine subspace "<<k<<"  cb " <<_Aggregate.subspace[k].checkerboard<<std::endl;
    }
    _Aggregate.Orthogonalise();
  }
  void calcCoarse(RealD alpha, RealD beta,int Npoly,
		  int Nk, int Nm,RealD resid, 
		  RealD MaxIt, RealD betastp, int MinRes)
  {
    Chebyshev<FineField>                          Cheby(alpha,beta,Npoly);
    ProjectedHermOp<Fobj,CComplex,nbasis>         Op(_FineOp,_Aggregate);
    ProjectedFunctionHermOp<Fobj,CComplex,nbasis> ChebyOp(Cheby,_FineOp,_Aggregate);

    std::vector<RealD>          eval(Nm);
    std::vector<CoarseField>    evec(Nm,_CoarseGrid);

    CoarseField src(_CoarseGrid);     src=1.0; 

    ImplicitlyRestartedLanczos<CoarseField> IRL(ChebyOp,ChebyOp,Nk,Nk,Nm,resid,MaxIt,betastp,MinRes);
    IRL.calc(eval,evec,src,Nk,false);

    // We got the evalues of the Cheby operator;
    // Reconstruct eigenvalues of original operator via Chebyshev inverse
    for (int i=0;i<Nk;i++){

      RealD eval_guess;
      if (i==0) eval_guess = 0;
      else      eval[i-1]  = 0;

      RealD eval_poly  = eval[i];
      RealD eval_op    = Cheby.approxInv(eval_poly,eval_guess,100,1e-10);
      std::cout << i << " Reconstructed eval = " << eval_op << " from quess " << eval_op << std::endl;
      eval[i] = eval_op;
    }

  }
};

struct LanczosParams : Serializable {
 public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(LanczosParams,
				  RealD, alpha, 
				  RealD, beta,
				  int, Npoly,
				  int, Nk,
				  int, Nm,
				  RealD, resid, 
				  int, MaxIt, 
				  RealD, betastp, 
				  int, MinRes);
};
struct CompressedLanczosParams : Serializable {
 public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(CompressedLanczosParams,
				  LanczosParams, FineParams,
				  LanczosParams, CoarseParams,
				  std::vector<int>, blockSize,
				  std::string, config,
				  std::vector < std::complex<double>  >, omega,
				  RealD, mass,
				  RealD, M5
				  );
};

int main (int argc, char ** argv) {

  Grid_init(&argc,&argv);

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
  //  FieldMetaData header;
  //  NerscIO::readConfiguration(Umu,header,Params.config);
  {
    std::vector<int> seeds4({1,2,3,4});
    GridParallelRNG          RNG4(UGrid);  RNG4.SeedFixedIntegers(seeds4);
    SU3::HotConfiguration(RNG4, Umu);
  }
  std::cout << GridLogMessage << "Lattice dimensions: " << GridDefaultLatt() << "   Ls: " << Ls << std::endl;

  // ZMobius EO Operator
  ZMobiusFermionR Ddwf(Umu, *FGrid, *FrbGrid, *UGrid, *UrbGrid, mass, M5, Params.omega,1.,0.);
  SchurDiagTwoOperator<ZMobiusFermionR,LatticeFermion> HermOp(Ddwf);

  // Eigenvector storage
  LanczosParams fine  =Params.FineParams;  
  LanczosParams coarse=Params.CoarseParams;  
  const int Nm1 = fine.Nm;
  const int Nm2 = coarse.Nm;

  std::cout << GridLogMessage << "Keep " << fine.Nk   << " full vectors" << std::endl;
  std::cout << GridLogMessage << "Keep " << coarse.Nk << " total vectors" << std::endl;
  assert(Nm2 >= Nm1);

  const int nbasis= 70;
  CoarseFineIRL<vSpinColourVector,vTComplex,nbasis> IRL(FrbGrid,CoarseGrid5rb,HermOp,Odd);

  std::cout << GridLogMessage << "Constructed CoarseFine IRL" << std::endl;

  std::cout << GridLogMessage << "Performing fine grid IRL Nk "<< nbasis<<" Nm "<<Nm1<< std::endl;
  IRL.testFine();
  //  IRL.calcFine(fine.alpha,fine.beta,fine.Npoly,fine.Nm,fine.resid, fine.MaxIt, fine.betastp, fine.MinRes);

  std::cout << GridLogMessage << "Performing coarse grid (poly) IRL " << nbasis<<" Nm "<<Nm2<< std::endl;
  IRL.calcCoarse(coarse.alpha,coarse.beta,coarse.Npoly,coarse.Nk,coarse.Nm,coarse.resid, coarse.MaxIt, coarse.betastp, coarse.MinRes);

  //  IRL.smoothedCoarseEigenvalues();

  Grid_finalize();
}

