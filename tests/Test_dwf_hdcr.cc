#include <Grid.h>
#include <algorithms/iterative/PrecGeneralisedConjugateResidual.h>
//#include <algorithms/iterative/PrecConjugateResidual.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

RealD InverseApproximation(RealD x){
  return 1.0/x;
}

template<class Fobj,class CComplex,int nbasis, class Matrix>
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
  Matrix         & _Matrix;
  FineOperator   & _FineOperator;

  // Constructor
  MultiGridPreconditioner(Aggregates &Agg, CoarseOperator &Coarse, FineOperator &Fine,Matrix &FineMatrix) 
    : _Aggregates(Agg),
      _CoarseOperator(Coarse),
      _FineOperator(Fine),
      _Matrix(FineMatrix)
  {
  }

  void PowerMethod(const FineField &in) {

    FineField p1(in._grid);
    FineField p2(in._grid);

    MdagMLinearOperator<Matrix,FineField>   fMdagMOp(_Matrix);

    p1=in;
    RealD absp2;
    for(int i=0;i<20;i++){
      RealD absp1=std::sqrt(norm2(p1));
      fMdagMOp.HermOp(p1,p2);// this is the G5 herm bit      
      //      _FineOperator.Op(p1,p2);// this is the G5 herm bit      
      RealD absp2=std::sqrt(norm2(p2));
      if(i%10==9)
	std::cout<<GridLogMessage << "Power method on mdagm "<<i<<" " << absp2/absp1<<std::endl;
      p1=p2*(1.0/std::sqrt(absp2));
    }
  }

#if 0
  void operator()(const FineField &in, FineField & out) {

    FineField Min(in._grid);
    FineField tmp(in._grid);

    CoarseVector Csrc(_CoarseOperator.Grid());
    CoarseVector Ctmp(_CoarseOperator.Grid());
    CoarseVector Csol(_CoarseOperator.Grid());

    // Monitor completeness of low mode space
    _Aggregates.ProjectToSubspace  (Csrc,in);
    _Aggregates.PromoteFromSubspace(Csrc,out);
    std::cout<<GridLogMessage<<"Completeness: "<<std::sqrt(norm2(out)/norm2(in))<<std::endl;

    // Build some solvers
    ConjugateGradient<FineField>    fCG(1.0e-3,1000);
    ConjugateGradient<CoarseVector>  CG(1.0e-8,100000);

    ////////////////////////////////////////////////////////////////////////
    // ADEF2: [PTM+Q] in = [1 - Q A] M in + Q in = Min + Q [ in -A Min]
    ////////////////////////////////////////////////////////////////////////

    // Smoothing step, followed by coarse grid correction
    MdagMLinearOperator<Matrix,FineField> MdagMOp(_Matrix);

    Min=in;
    std::cout<<GridLogMessage<< " Preconditioner in  " << norm2(in)<<std::endl; 
    _FineOperator.AdjOp(Min,tmp);
    std::cout<<GridLogMessage<< " Preconditioner tmp  " << norm2(in)<<std::endl; 

    fCG(MdagMOp,tmp,out);

    _FineOperator.Op(out,tmp);

    std::cout<<GridLogMessage<< " Preconditioner in  " << norm2(in)<<std::endl; 
    std::cout<<GridLogMessage<< " Preconditioner out " << norm2(out)<<std::endl; 
    std::cout<<GridLogMessage<< " Preconditioner Aout" << norm2(tmp)<<std::endl; 

    tmp = tmp - in;
    
    std::cout<<GridLogMessage<<"preconditioner thinks residual is "<<std::sqrt(norm2(tmp)/norm2(in))<<std::endl;

    /*
    //    _FineOperator.Op(Min,out);
    //    out = in -out; // out = in - A Min
    out = in;

        MdagMLinearOperator<CoarseOperator,CoarseVector> MdagMOp(_CoarseOperator);
    HermitianLinearOperator<CoarseOperator,CoarseVector> HermOp(_CoarseOperator);
    Csol=zero;
    _Aggregates.ProjectToSubspace  (Csrc,out);
    HermOp.AdjOp(Csrc,Ctmp);// Normal equations
    CG(MdagMOp  ,Ctmp,Csol);
    _Aggregates.PromoteFromSubspace(Csol,out);

    out = Min + out;;
    */

  }
#endif

    ////////////////////////////////////////////////////////////////////////
    // ADEF2: [PTM+Q] in = [1 - Q A] M in + Q in = Min + Q [ in -A Min]
    // ADEF1: [MP+Q ] in =M [1 - A Q] in + Q in  
    ////////////////////////////////////////////////////////////////////////
#if 0
  void operator()(const FineField &in, FineField & out) {

    CoarseVector Csrc(_CoarseOperator.Grid());
    CoarseVector Ctmp(_CoarseOperator.Grid());
    CoarseVector Csol(_CoarseOperator.Grid());

    ConjugateGradient<CoarseVector>  CG(1.0e-10,100000);
    ConjugateGradient<FineField>    fCG(3.0e-2,1000);

    HermitianLinearOperator<CoarseOperator,CoarseVector>  HermOp(_CoarseOperator);
    MdagMLinearOperator<CoarseOperator,CoarseVector>     MdagMOp(_CoarseOperator);
    MdagMLinearOperator<Matrix,FineField>               fMdagMOp(_Matrix);

    FineField tmp(in._grid);
    FineField res(in._grid);
    FineField Min(in._grid);

    // Monitor completeness of low mode space
    _Aggregates.ProjectToSubspace  (Csrc,in);
    _Aggregates.PromoteFromSubspace(Csrc,out);
    std::cout<<GridLogMessage<<"Coarse Grid Preconditioner\nCompleteness in: "<<std::sqrt(norm2(out)/norm2(in))<<std::endl;

    // [PTM+Q] in = [1 - Q A] M in + Q in = Min + Q [ in -A Min]
    _FineOperator.Op(in,tmp);// this is the G5 herm bit
    fCG(fMdagMOp,tmp,Min);    // solves  MdagM = g5 M g5M

    // Monitor completeness of low mode space
    _Aggregates.ProjectToSubspace  (Csrc,Min);
    _Aggregates.PromoteFromSubspace(Csrc,out);
    std::cout<<GridLogMessage<<"Completeness Min: "<<std::sqrt(norm2(out)/norm2(Min))<<std::endl;

    _FineOperator.Op(Min,tmp);
    tmp = in - tmp;   // in - A Min

    Csol=zero;
    _Aggregates.ProjectToSubspace  (Csrc,tmp);
    HermOp.AdjOp(Csrc,Ctmp);// Normal equations
    CG(MdagMOp,Ctmp,Csol);

    HermOp.Op(Csol,Ctmp);
    Ctmp=Ctmp-Csrc;
    std::cout<<GridLogMessage<<"coarse space true residual "<<std::sqrt(norm2(Ctmp)/norm2(Csrc))<<std::endl;
    _Aggregates.PromoteFromSubspace(Csol,out);

    _FineOperator.Op(out,res);
    res=res-tmp;
    std::cout<<GridLogMessage<<"promoted sol residual "<<std::sqrt(norm2(res)/norm2(tmp))<<std::endl;
    _Aggregates.ProjectToSubspace  (Csrc,res);
    std::cout<<GridLogMessage<<"coarse space proj of residual "<<norm2(Csrc)<<std::endl;

    
    out = out+Min; // additive coarse space correction
    //    out = Min; // no additive coarse space correction

    _FineOperator.Op(out,tmp);
    tmp=tmp-in;         // tmp is new residual

    std::cout<<GridLogMessage<< " Preconditioner in  " << norm2(in)<<std::endl; 
    std::cout<<GridLogMessage<< " Preconditioner out " << norm2(out)<<std::endl; 
    std::cout<<GridLogMessage<<"preconditioner thinks residual is "<<std::sqrt(norm2(tmp)/norm2(in))<<std::endl;

  }
#endif
  // ADEF1: [MP+Q ] in =M [1 - A Q] in + Q in  
#if 0
  void operator()(const FineField &in, FineField & out) {

    CoarseVector Csrc(_CoarseOperator.Grid());
    CoarseVector Ctmp(_CoarseOperator.Grid());
    CoarseVector Csol(_CoarseOperator.Grid()); Csol=zero;

    ConjugateGradient<CoarseVector>  CG(1.0e-10,100000);
    ConjugateGradient<FineField>    fCG(3.0e-2,1000);

    HermitianLinearOperator<CoarseOperator,CoarseVector>  HermOp(_CoarseOperator);
    MdagMLinearOperator<CoarseOperator,CoarseVector>     MdagMOp(_CoarseOperator);
    ShiftedMdagMLinearOperator<Matrix,FineField>        fMdagMOp(_Matrix,0.1);

    FineField tmp(in._grid);
    FineField res(in._grid);
    FineField Qin(in._grid);

    // Monitor completeness of low mode space
    //    _Aggregates.ProjectToSubspace  (Csrc,in);
    //    _Aggregates.PromoteFromSubspace(Csrc,out);
    //    std::cout<<GridLogMessage<<"Coarse Grid Preconditioner\nCompleteness in: "<<std::sqrt(norm2(out)/norm2(in))<<std::endl;
    
    _Aggregates.ProjectToSubspace  (Csrc,in);
    HermOp.AdjOp(Csrc,Ctmp);// Normal equations
    CG(MdagMOp,Ctmp,Csol);
    _Aggregates.PromoteFromSubspace(Csol,Qin);

    //    Qin=0;
    _FineOperator.Op(Qin,tmp);// A Q in
    tmp = in - tmp;            // in - A Q in

    _FineOperator.Op(tmp,res);// this is the G5 herm bit
    fCG(fMdagMOp,res,out);    // solves  MdagM = g5 M g5M

    out = out + Qin;

    _FineOperator.Op(out,tmp);
    tmp=tmp-in;         // tmp is new residual

    std::cout<<GridLogMessage<<"preconditioner thinks residual is "<<std::sqrt(norm2(tmp)/norm2(in))<<std::endl;

  }
#endif

  void SmootherTest (const FineField & in){
    
    FineField vec1(in._grid);
    FineField vec2(in._grid);
    RealD lo[3] = { 0.5, 1.0, 2.0};

    //    MdagMLinearOperator<Matrix,FineField>        fMdagMOp(_Matrix);
    ShiftedMdagMLinearOperator<Matrix,FineField> fMdagMOp(_Matrix,0.5);

    RealD Ni,r;

    Ni = norm2(in);

    for(int ilo=0;ilo<3;ilo++){
      for(int ord=5;ord<50;ord*=2){

	_FineOperator.AdjOp(in,vec1);

	Chebyshev<FineField> Cheby  (lo[ilo],70.0,ord,InverseApproximation);
	Cheby(fMdagMOp,vec1,vec2);    // solves  MdagM = g5 M g5M

	_FineOperator.Op(vec2,vec1);// this is the G5 herm bit
	vec1  = in - vec1;   // tmp  = in - A Min
	r=norm2(vec1);
	std::cout<<GridLogMessage << "Smoother resid "<<std::sqrt(r/Ni)<<std::endl;

      }
    }
  }

  void operator()(const FineField &in, FineField & out) {

    CoarseVector Csrc(_CoarseOperator.Grid());
    CoarseVector Ctmp(_CoarseOperator.Grid());
    CoarseVector Csol(_CoarseOperator.Grid()); Csol=zero;

    ConjugateGradient<CoarseVector>  CG(1.0e-3,100000);
    //    ConjugateGradient<FineField>    fCG(3.0e-2,1000);

    HermitianLinearOperator<CoarseOperator,CoarseVector>  HermOp(_CoarseOperator);
    MdagMLinearOperator<CoarseOperator,CoarseVector>     MdagMOp(_CoarseOperator);
    //    MdagMLinearOperator<Matrix,FineField>        fMdagMOp(_Matrix);
    ShiftedMdagMLinearOperator<Matrix,FineField> fMdagMOp(_Matrix,0.0);

    FineField vec1(in._grid);
    FineField vec2(in._grid);

    //    Chebyshev<FineField> Cheby    (0.5,70.0,30,InverseApproximation);
    //    Chebyshev<FineField> ChebyAccu(0.5,70.0,30,InverseApproximation);
    Chebyshev<FineField> Cheby    (2.0,70.0,10,InverseApproximation);
    Chebyshev<FineField> ChebyAccu(2.0,70.0,10,InverseApproximation);
    Cheby.JacksonSmooth();
    ChebyAccu.JacksonSmooth();

    _Aggregates.ProjectToSubspace  (Csrc,in);
    _Aggregates.PromoteFromSubspace(Csrc,out);
    std::cout<<GridLogMessage<<"Completeness: "<<std::sqrt(norm2(out)/norm2(in))<<std::endl;
    
    //    ofstream fout("smoother");
    //    Cheby.csv(fout);

    // V11 multigrid.
    // Use a fixed chebyshev and hope hermiticity helps.

    // To make a working smoother for indefinite operator
    // must multiply by "Mdag" (ouch loses all low mode content)
    // and apply to poly approx of (mdagm)^-1.
    // so that we end up with an odd polynomial.

    RealD Ni = norm2(in);

    _FineOperator.AdjOp(in,vec1);// this is the G5 herm bit
    ChebyAccu(fMdagMOp,vec1,out);    // solves  MdagM = g5 M g5M

    std::cout<<GridLogMessage << "Smoother norm "<<norm2(out)<<std::endl;

    // Update with residual for out
    _FineOperator.Op(out,vec1);// this is the G5 herm bit
    vec1  = in - vec1;   // tmp  = in - A Min

    RealD r = norm2(vec1);

    std::cout<<GridLogMessage << "Smoother resid "<<std::sqrt(r/Ni)<< " " << r << " " << Ni <<std::endl;
    
    _Aggregates.ProjectToSubspace  (Csrc,vec1);
    HermOp.AdjOp(Csrc,Ctmp);// Normal equations
    CG(MdagMOp,Ctmp,Csol);
    _Aggregates.PromoteFromSubspace(Csol,vec1); // Ass^{-1} [in - A Min]_s
                                             // Q = Q[in - A Min]  
    out = out+vec1;

    // Three preconditioner smoothing -- hermitian if C3 = C1
    // Recompute error
    _FineOperator.Op(out,vec1);// this is the G5 herm bit
    vec1  = in - vec1;   // tmp  = in - A Min
    r=norm2(vec1);

    std::cout<<GridLogMessage << "Coarse resid "<<std::sqrt(r/Ni)<<std::endl;

    // Reapply smoother
    _FineOperator.Op(vec1,vec2);  // this is the G5 herm bit
    ChebyAccu(fMdagMOp,vec2,vec1);    // solves  MdagM = g5 M g5M

    out =out+vec1;
    _FineOperator.Op(out,vec1);// this is the G5 herm bit
    vec1  = in - vec1;   // tmp  = in - A Min
    r=norm2(vec1);
    std::cout<<GridLogMessage << "Smoother resid "<<std::sqrt(r/Ni)<<std::endl;

  }

};

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  const int Ls=8;

  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);

  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  ///////////////////////////////////////////////////
  // Construct a coarsened grid; utility for this?
  ///////////////////////////////////////////////////
  const int block=2;
  std::vector<int> clatt = GridDefaultLatt();
  for(int d=0;d<clatt.size();d++){
    clatt[d] = clatt[d]/block;
  }
  GridCartesian *Coarse4d =  SpaceTimeGrid::makeFourDimGrid(clatt, GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());;
  GridCartesian *Coarse5d =  SpaceTimeGrid::makeFiveDimGrid(1,Coarse4d);

  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});
  std::vector<int> cseeds({5,6,7,8});
  GridParallelRNG          RNG5(FGrid);   RNG5.SeedFixedIntegers(seeds5);
  GridParallelRNG          RNG4(UGrid);   RNG4.SeedFixedIntegers(seeds4);
  GridParallelRNG          CRNG(Coarse5d);CRNG.SeedFixedIntegers(cseeds);

  Gamma g5(Gamma::Gamma5);

  LatticeFermion    src(FGrid); gaussian(RNG5,src);// src=src+g5*src;
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

  RealD mass=0.01;
  RealD M5=1.8;

  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  std::cout<<GridLogMessage << "Building g5R5 hermitian DWF operator" <<std::endl;
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  DomainWallFermionR Ddwf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);

  const int nbasis = 32;
  //  const int nbasis = 4;

  typedef Aggregation<vSpinColourVector,vTComplex,nbasis>              Subspace;
  typedef CoarsenedMatrix<vSpinColourVector,vTComplex,nbasis>          CoarseOperator;
  typedef CoarseOperator::CoarseVector                                 CoarseVector;

  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  std::cout<<GridLogMessage << "Calling Aggregation class to build subspace" <<std::endl;
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  MdagMLinearOperator<DomainWallFermionR,LatticeFermion> HermDefOp(Ddwf);
  Subspace Aggregates(Coarse5d,FGrid);
  //  Aggregates.CreateSubspace(RNG5,HermDefOp,nbasis);
  assert ( (nbasis & 0x1)==0);
  int nb=nbasis/2;
  std::cout<<GridLogMessage << " nbasis/2 = "<<nb<<std::endl;
  Aggregates.CreateSubspace(RNG5,HermDefOp,nb);
  for(int n=0;n<nb;n++){
    G5R5(Aggregates.subspace[n+nb],Aggregates.subspace[n]);
    std::cout<<GridLogMessage<<n<<" subspace "<<norm2(Aggregates.subspace[n+nb])<<" "<<norm2(Aggregates.subspace[n]) <<std::endl;
  }
  for(int n=0;n<nbasis;n++){
    std::cout<<GridLogMessage << "vec["<<n<<"] = "<<norm2(Aggregates.subspace[n])  <<std::endl;
  }

//  for(int i=0;i<nbasis;i++){
//    result =     Aggregates.subspace[i];
//    Aggregates.subspace[i]=result+g5*result;
//  }
  result=zero;
  
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  std::cout<<GridLogMessage << "Building coarse representation of Indef operator" <<std::endl;
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  Gamma5R5HermitianLinearOperator<DomainWallFermionR,LatticeFermion> HermIndefOp(Ddwf);
  CoarsenedMatrix<vSpinColourVector,vTComplex,nbasis> LDOp(*Coarse5d);
  LDOp.CoarsenOperator(FGrid,HermIndefOp,Aggregates);

  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  std::cout<<GridLogMessage << "Testing some coarse space solvers  " <<std::endl;
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  CoarseVector c_src (Coarse5d);
  CoarseVector c_res (Coarse5d);
  gaussian(CRNG,c_src);
  c_res=zero;

  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  std::cout<<GridLogMessage << "Solving posdef-CG on coarse space "<< std::endl;
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  MdagMLinearOperator<CoarseOperator,CoarseVector> PosdefLdop(LDOp);
  ConjugateGradient<CoarseVector> CG(1.0e-6,100000);
  CG(PosdefLdop,c_src,c_res);

  //  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  //  std::cout<<GridLogMessage << "Solving indef-MCR on coarse space "<< std::endl;
  //  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  //  HermitianLinearOperator<CoarseOperator,CoarseVector> HermIndefLdop(LDOp);
  //  ConjugateResidual<CoarseVector> MCR(1.0e-6,100000);
  //MCR(HermIndefLdop,c_src,c_res);

  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  std::cout<<GridLogMessage << "Building deflation preconditioner "<< std::endl;
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;

  MultiGridPreconditioner <vSpinColourVector,vTComplex,nbasis,DomainWallFermionR> Precon(Aggregates, LDOp,HermIndefOp,Ddwf);
  TrivialPrecon<LatticeFermion> simple;

  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  std::cout<<GridLogMessage << "Testing smoother efficacy"<< std::endl;
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  Precon.SmootherTest(src);

  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  std::cout<<GridLogMessage << "Unprec CG "<< std::endl;
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  //  TrivialPrecon<LatticeFermion> simple;
  //  ConjugateGradient<LatticeFermion> fCG(1.0e-8,100000);
  //  fCG(HermDefOp,src,result);
  //  exit(0);

  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  std::cout<<GridLogMessage << "Testing GCR on indef matrix "<< std::endl;
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  //  PrecGeneralisedConjugateResidual<LatticeFermion> UPGCR(1.0e-8,100000,simple,8,128);
  //  UPGCR(HermIndefOp,src,result);

  
  /// Get themax eval
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  std::cout<<GridLogMessage <<" Applying power method to find spectral range      "<<std::endl;
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  Precon.PowerMethod(src);

  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  std::cout<<GridLogMessage << "Building a two level PGCR "<< std::endl;
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  PrecGeneralisedConjugateResidual<LatticeFermion> PGCR(1.0e-8,100000,Precon,8,128);
  std::cout<<GridLogMessage<<"checking norm src "<<norm2(src)<<std::endl;
  PGCR(HermIndefOp,src,result);

  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  std::cout<<GridLogMessage << "Red Black Prec CG "<< std::endl;
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  SchurDiagMooeeOperator<DomainWallFermionR,LatticeFermion> HermOpEO(Ddwf);
  ConjugateGradient<LatticeFermion> pCG(1.0e-8,10000);

  LatticeFermion    src_o(FrbGrid);
  LatticeFermion result_o(FrbGrid);
  pickCheckerboard(Odd,src_o,src);
  result_o=zero;

  pCG(HermOpEO,src_o,result_o);

  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  std::cout<<GridLogMessage << "Done "<< std::endl;
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  Grid_finalize();
}
