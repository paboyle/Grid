#include <Grid.h>
#include <algorithms/iterative/PrecGeneralisedConjugateResidual.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

template<class Fobj,class CComplex,int nbasis>
class MultiGridPreconditioner : public LinearFunction< Lattice<Fobj> > {
public:

  typedef Aggregation<Fobj,CComplex,nbasis> Aggregates;
  typedef CoarsenedMatrix<Fobj,CComplex,nbasis> CoarseOperator;

  typedef typename Aggregation<Fobj,CComplex,nbasis>::siteVector     siteVector;
  typedef typename Aggregation<Fobj,CComplex,nbasis>::CoarseScalar CoarseScalar;
  typedef typename Aggregation<Fobj,CComplex,nbasis>::CoarseVector CoarseVector;
  typedef typename Aggregation<Fobj,CComplex,nbasis>::CoarseMatrix CoarseMatrix;
  typedef typename Aggregation<Fobj,CComplex,nbasis>::FineField    FineField;
  typedef LinearOperatorBase<FineField>                            FineOperator;

  Aggregates     & _Aggregates;
  CoarseOperator & _CoarseOperator;
  FineOperator   & _FineOperator;

  // Constructor
  MultiGridPreconditioner(Aggregates &Agg, CoarseOperator &Coarse, FineOperator &Fine) 
    : _Aggregates(Agg),
      _CoarseOperator(Coarse),
      _FineOperator(Fine)
  {
  }

  void operator()(const FineField &in, FineField & out) {

    FineField Min(in._grid);

    CoarseVector Csrc(_CoarseOperator.Grid());
    CoarseVector Ctmp(_CoarseOperator.Grid());
    CoarseVector Csol(_CoarseOperator.Grid());

    // Monitor completeness of low mode space
    _Aggregates.ProjectToSubspace  (Csrc,in);
    _Aggregates.PromoteFromSubspace(Csrc,out);
    std::cout<<"Completeness: "<<std::sqrt(norm2(out)/norm2(in))<<std::endl;

    ConjugateResidual<FineField>    MCR(1.0e-2,1000);
    ConjugateGradient<CoarseVector>  CG(1.0e-2,10000);

    // [PTM+Q] in = [1 - Q A] M in + Q in = Min + Q [ in -A Min]
    // Smoothing step, followed by coarse grid correction

    MCR(_FineOperator,in,Min);
    _FineOperator.Op(Min,out);
    out = in -out; // out = in - A Min

    MdagMLinearOperator<CoarseOperator,CoarseVector>     MdagMOp(_CoarseOperator);
    HermitianLinearOperator<CoarseOperator,CoarseVector> HermOp(_CoarseOperator);
    Csol=zero;
    _Aggregates.ProjectToSubspace  (Csrc,out);
    HermOp.AdjOp(Csrc,Ctmp);// Normal equations
    CG(MdagMOp  ,Ctmp,Csol);
    _Aggregates.PromoteFromSubspace(Csol,out);

    out = Min + out;;
  }

};

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  const int Ls=8;

  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplexF::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);

  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  ///////////////////////////////////////////////////
  // Construct a coarsened grid; utility for this?
  ///////////////////////////////////////////////////
  std::vector<int> clatt = GridDefaultLatt();
  for(int d=0;d<clatt.size();d++){
    clatt[d] = clatt[d]/4;
  }
  GridCartesian *Coarse4d =  SpaceTimeGrid::makeFourDimGrid(clatt, GridDefaultSimd(Nd,vComplexF::Nsimd()),GridDefaultMpi());;
  GridCartesian *Coarse5d =  SpaceTimeGrid::makeFiveDimGrid(1,Coarse4d);

  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});
  std::vector<int> cseeds({5,6,7,8});
  GridParallelRNG          RNG5(FGrid);   RNG5.SeedFixedIntegers(seeds5);
  GridParallelRNG          RNG4(UGrid);   RNG4.SeedFixedIntegers(seeds4);
  GridParallelRNG          CRNG(Coarse5d);CRNG.SeedFixedIntegers(cseeds);

  LatticeFermion    src(FGrid); gaussian(RNG5,src);
  LatticeFermion result(FGrid); result=zero;
  LatticeFermion    ref(FGrid); ref=zero;
  LatticeFermion    tmp(FGrid);
  LatticeFermion    err(FGrid);
  LatticeGaugeField Umu(UGrid); 

  NerscField header;
  std::string file("./ckpoint_lat.4000");
  readNerscConfiguration(Umu,header,file);

  //  SU3::ColdConfiguration(RNG4,Umu);
  //  SU3::TepidConfiguration(RNG4,Umu);
  //  SU3::HotConfiguration(RNG4,Umu);
  //  Umu=zero;

  RealD mass=0.04;
  RealD M5=1.8;

  std::cout << "**************************************************"<< std::endl;
  std::cout << "Building g5R5 hermitian DWF operator" <<std::endl;
  std::cout << "**************************************************"<< std::endl;
  DomainWallFermion Ddwf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);

  const int nbasis = 4;

  typedef Aggregation<vSpinColourVector,vTComplex,nbasis>              Subspace;
  typedef CoarsenedMatrix<vSpinColourVector,vTComplex,nbasis>          CoarseOperator;
  typedef CoarseOperator::CoarseVector                                 CoarseVector;

  std::cout << "**************************************************"<< std::endl;
  std::cout << "Calling Aggregation class to build subspace" <<std::endl;
  std::cout << "**************************************************"<< std::endl;
  MdagMLinearOperator<DomainWallFermion,LatticeFermion> HermDefOp(Ddwf);
  Subspace Aggregates(Coarse5d,FGrid);
  Aggregates.CreateSubspace(RNG5,HermDefOp);

  std::cout << "**************************************************"<< std::endl;
  std::cout << "Building coarse representation of Indef operator" <<std::endl;
  std::cout << "**************************************************"<< std::endl;
  Gamma5R5HermitianLinearOperator<DomainWallFermion,LatticeFermion> HermIndefOp(Ddwf);
  CoarsenedMatrix<vSpinColourVector,vTComplex,nbasis> LDOp(*Coarse5d);
  LDOp.CoarsenOperator(FGrid,HermIndefOp,Aggregates);

  std::cout << "**************************************************"<< std::endl;
  std::cout << "Testing some coarse space solvers  " <<std::endl;
  std::cout << "**************************************************"<< std::endl;
  CoarseVector c_src (Coarse5d);
  CoarseVector c_res (Coarse5d);
  gaussian(CRNG,c_src);
  c_res=zero;

  std::cout << "**************************************************"<< std::endl;
  std::cout << "Solving posdef-CG on coarse space "<< std::endl;
  std::cout << "**************************************************"<< std::endl;
  MdagMLinearOperator<CoarseOperator,CoarseVector> PosdefLdop(LDOp);
  ConjugateGradient<CoarseVector> CG(1.0e-6,10000);
  CG(PosdefLdop,c_src,c_res);

  std::cout << "**************************************************"<< std::endl;
  std::cout << "Solving indef-MCR on coarse space "<< std::endl;
  std::cout << "**************************************************"<< std::endl;
  HermitianLinearOperator<CoarseOperator,CoarseVector> HermIndefLdop(LDOp);
  ConjugateResidual<CoarseVector> MCR(1.0e-6,10000);
  //MCR(HermIndefLdop,c_src,c_res);

  std::cout << "**************************************************"<< std::endl;
  std::cout << "Building deflation preconditioner "<< std::endl;
  std::cout << "**************************************************"<< std::endl;

  MultiGridPreconditioner <vSpinColourVector,vTComplex,nbasis> Precon(Aggregates, LDOp,HermIndefOp);

  std::cout << "**************************************************"<< std::endl;
  std::cout << "Building a one level PGCR "<< std::endl;
  std::cout << "**************************************************"<< std::endl;
  TrivialPrecon<LatticeFermion> simple;
  PrecGeneralisedConjugateResidual<LatticeFermion> GCR(1.0e-6,10000,simple,8,64);
  GCR(HermIndefOp,src,result);

  std::cout << "**************************************************"<< std::endl;
  std::cout << "Building a two level PGCR "<< std::endl;
  std::cout << "**************************************************"<< std::endl;
  PrecGeneralisedConjugateResidual<LatticeFermion> PGCR(1.0e-6,10000,Precon,8,64);
  PGCR(HermIndefOp,src,result);

  std::cout << "**************************************************"<< std::endl;
  std::cout << "Done "<< std::endl;
  std::cout << "**************************************************"<< std::endl;
  Grid_finalize();
}
