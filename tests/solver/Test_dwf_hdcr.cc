/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_dwf_hdcr.cc

    Copyright (C) 2015

Author: Antonin Portelli <antonin.portelli@me.com>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
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
#include <Grid/algorithms/iterative/PrecGeneralisedConjugateResidual.h>
//#include <algorithms/iterative/PrecConjugateResidual.h>

using namespace std;
using namespace Grid;

class myclass: Serializable {
public:

  GRID_SERIALIZABLE_CLASS_MEMBERS(myclass,
			  int, domaindecompose,
			  int, domainsize,
			  int, order,
			  int, Ls,
			  double, mq,
			  double, lo,
			  double, hi,
			  int, steps);

  myclass(){};

};

RealD InverseApproximation(RealD x){
  return 1.0/x;
}

template<class Fobj,class CComplex,int nbasis, class Matrix, class Guesser>
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
  Matrix         & _FineMatrix;
  FineOperator   & _FineOperator;
  Matrix         & _SmootherMatrix;
  FineOperator   & _SmootherOperator;
  Guesser        & _Guess;

  double cheby_hi;
  double cheby_lo;
  int    cheby_ord;

  myclass _params;

  // Constructor
  MultiGridPreconditioner(Aggregates &Agg, CoarseOperator &Coarse, 
			  FineOperator &Fine,Matrix &FineMatrix,
			  FineOperator &Smooth,Matrix &SmootherMatrix,
			  Guesser &Guess_,
			  myclass params_)
    : _Aggregates(Agg),
      _CoarseOperator(Coarse),
      _FineOperator(Fine),
      _FineMatrix(FineMatrix),
      _SmootherOperator(Smooth),
      _SmootherMatrix(SmootherMatrix),
      _Guess(Guess_),
      _params(params_)
  {
  }

  void PowerMethod(const FineField &in) {

    FineField p1(in.Grid());
    FineField p2(in.Grid());

    MdagMLinearOperator<Matrix,FineField>   fMdagMOp(_FineMatrix);

    p1=in;
    for(int i=0;i<50;i++){
      RealD absp1=std::sqrt(norm2(p1));
      fMdagMOp.HermOp(p1,p2);// this is the G5 herm bit      
      //      _FineOperator.Op(p1,p2);// this is the G5 herm bit      
      RealD absp2=std::sqrt(norm2(p2));
      if(i%10==9)
	std::cout<<GridLogMessage << "Power method on mdagm "<<i<<" " << absp2/absp1<<std::endl;
      p1=p2*(1.0/std::sqrt(absp2));
    }
  }

  void operator()(const FineField &in, FineField & out ) {
    operatorCheby(in,out);
    //operatorADEF2(in,out);
  }

    ////////////////////////////////////////////////////////////////////////
    // ADEF2: [PTM+Q] in = [1 - Q A] M in + Q in = Min + Q [ in -A Min]
    // ADEF1: [MP+Q ] in =M [1 - A Q] in + Q in  
    ////////////////////////////////////////////////////////////////////////
#if 1
  void operatorADEF2(const FineField &in, FineField & out) {

    CoarseVector Csrc(_CoarseOperator.Grid());
    CoarseVector Ctmp(_CoarseOperator.Grid());
    CoarseVector Csol(_CoarseOperator.Grid());

    ConjugateGradient<CoarseVector>  CG(1.0e-3,100,false);
    ConjugateGradient<FineField>    fCG(1.0e-3,10,false);

    HermitianLinearOperator<CoarseOperator,CoarseVector>  HermOp(_CoarseOperator);
    MdagMLinearOperator<CoarseOperator,CoarseVector>     MdagMOp(_CoarseOperator);
    MdagMLinearOperator<Matrix,FineField>               fMdagMOp(_FineMatrix);

    FineField tmp(in.Grid());
    FineField res(in.Grid());
    FineField Min(in.Grid());

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

    _Aggregates.ProjectToSubspace  (Csrc,tmp);
    HermOp.AdjOp(Csrc,Ctmp);// Normal equations
    _Guess(Ctmp,Csol);
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
#if 1
  void operatorADEF1(const FineField &in, FineField & out) {

    CoarseVector Csrc(_CoarseOperator.Grid());
    CoarseVector Ctmp(_CoarseOperator.Grid());
    CoarseVector Csol(_CoarseOperator.Grid()); Csol=Zero();

    ConjugateGradient<CoarseVector>  CG(1.0e-10,100000);
    ConjugateGradient<FineField>    fCG(1.0e-3,1000);

    HermitianLinearOperator<CoarseOperator,CoarseVector>  HermOp(_CoarseOperator);
    MdagMLinearOperator<CoarseOperator,CoarseVector>     MdagMOp(_CoarseOperator);
    ShiftedMdagMLinearOperator<Matrix,FineField>        fMdagMOp(_FineMatrix,0.1);

    FineField tmp(in.Grid());
    FineField res(in.Grid());
    FineField Qin(in.Grid());

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
    
    FineField vec1(in.Grid());
    FineField vec2(in.Grid());
    RealD lo[3] = { 0.5, 1.0, 2.0};

    //    MdagMLinearOperator<Matrix,FineField>        fMdagMOp(_FineMatrix);
    ShiftedMdagMLinearOperator<Matrix,FineField> fMdagMOp(_SmootherMatrix,0.0);

    RealD Ni,r;

    Ni = norm2(in);

    for(int ilo=0;ilo<3;ilo++){
      for(int ord=5;ord<50;ord*=2){

	std::cout << " lo "<<lo[ilo]<<" order "<<ord<<std::endl;

	_SmootherOperator.AdjOp(in,vec1);

	Chebyshev<FineField> Cheby  (lo[ilo],70.0,ord,InverseApproximation);
	Cheby(fMdagMOp,vec1,vec2);    // solves  MdagM = g5 M g5M

	_FineOperator.Op(vec2,vec1);// this is the G5 herm bit
	vec1  = in - vec1;   // tmp  = in - A Min
	r=norm2(vec1);
	std::cout<<GridLogMessage << "Smoother resid "<<std::sqrt(r/Ni)<<std::endl;

      }
    }
  }

  void operatorCheby(const FineField &in, FineField & out) {

    CoarseVector Csrc(_CoarseOperator.Grid());
    CoarseVector Ctmp(_CoarseOperator.Grid());
    CoarseVector Ctmp1(_CoarseOperator.Grid());
    CoarseVector Csol(_CoarseOperator.Grid()); 
    
    ConjugateGradient<CoarseVector>  CG(5.0e-2,100000);

    HermitianLinearOperator<CoarseOperator,CoarseVector>  HermOp(_CoarseOperator);
    MdagMLinearOperator<CoarseOperator,CoarseVector>     MdagMOp(_CoarseOperator);
    //    MdagMLinearOperator<Matrix,FineField>        fMdagMOp(_FineMatrix);
    ShiftedMdagMLinearOperator<Matrix,FineField> fMdagMOp(_SmootherMatrix,0.0);

    FineField vec1(in.Grid());
    FineField vec2(in.Grid());

    Chebyshev<FineField> Cheby    (_params.lo,_params.hi,_params.order,InverseApproximation);
    Chebyshev<FineField> ChebyAccu(_params.lo,_params.hi,_params.order,InverseApproximation);

    //    _Aggregates.ProjectToSubspace  (Csrc,in);
    //    _Aggregates.PromoteFromSubspace(Csrc,out);
    //    std::cout<<GridLogMessage<<"Completeness: "<<std::sqrt(norm2(out)/norm2(in))<<std::endl;
    
    //    ofstream fout("smoother");
    //    Cheby.csv(fout);

    // V11 multigrid.
    // Use a fixed chebyshev and hope hermiticity helps.

    // To make a working smoother for indefinite operator
    // must multiply by "Mdag" (ouch loses all low mode content)
    // and apply to poly approx of (mdagm)^-1.
    // so that we end up with an odd polynomial.

    RealD Ni = norm2(in);

    std::cout<<GridLogMessage << "Smoother calling Cheby" <<std::endl;
    _SmootherOperator.AdjOp(in,vec1);// this is the G5 herm bit
    ChebyAccu(fMdagMOp,vec1,out);    // solves  MdagM = g5 M g5M
    std::cout<<GridLogMessage << "Smoother called Cheby" <<std::endl;

    // Update with residual for out
    _FineOperator.Op(out,vec1);// this is the G5 herm bit
    vec1  = in - vec1;   // tmp  = in - A Min

    RealD r = norm2(vec1);

    std::cout<<GridLogMessage << "Smoother resid "<<std::sqrt(r/Ni)<< " " << r << " " << Ni <<std::endl;
    
    std::cout<<GridLogMessage << "ProjectToSubspace" <<std::endl;
    _Aggregates.ProjectToSubspace  (Csrc,vec1);
    std::cout<<GridLogMessage << "ProjectToSubspaceDone" <<std::endl;
    
    HermOp.AdjOp(Csrc,Ctmp1);// Normal equations

    _Guess(Ctmp1,Csol);
    CG(MdagMOp,Ctmp1,Csol);

    //////////////////////////////
    // Recompute true residual
    //////////////////////////////
    MdagMOp.HermOp(Csol,Ctmp);
    Ctmp = Ctmp1 - Ctmp;      // r=Csrc - M^dagM sol // This is already computed inside CG
    HermOp.AdjOp(Ctmp,Ctmp1);// Normal equations
    _Guess(Ctmp1,Ctmp);      // sol = sol' + MdagM^-1 (Csrc' - MdagM sol')
    Csol = Csol + Ctmp;

    std::cout<<GridLogMessage << "PromoteFromSubspace" <<std::endl;
    _Aggregates.PromoteFromSubspace(Csol,vec1); // Ass^{-1} [in - A Min]_s
                                                // Q = Q[in - A Min]  
    std::cout<<GridLogMessage << "PromoteFromSubspaceDone" <<std::endl;
    out = out+vec1;

    // Three preconditioner smoothing -- hermitian if C3 = C1
    // Recompute error
    _FineOperator.Op(out,vec1);// this is the G5 herm bit
    std::cout<<GridLogMessage << "FineOp" <<std::endl;
    vec1  = in - vec1;   // tmp  = in - A Min
    r=norm2(vec1);

    std::cout<<GridLogMessage << "Coarse resid "<<std::sqrt(r/Ni)<<std::endl;

    // Reapply smoother
    std::cout<<GridLogMessage << "Smoother calling Cheby" <<std::endl;
    _SmootherOperator.Op(vec1,vec2);  // this is the G5 herm bit
    ChebyAccu(fMdagMOp,vec2,vec1);    // solves  MdagM = g5 M g5M
    std::cout<<GridLogMessage << "Smoother called Cheby" <<std::endl;

    out =out+vec1;
    vec1  = in - vec1;   // tmp  = in - A Min
    r=norm2(vec1);
    std::cout<<GridLogMessage << "Smoother resid "<<std::sqrt(r/Ni)<<std::endl;

  }

};

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  myclass params;
  myclass cparams;

  XmlReader RD("params.xml");
  read(RD,"params",params);
  std::cout<<"Params: Order "<<params.order<<"["<<params.lo<<","<<params.hi<<"]"<< " steps "<<params.steps<<std::endl;

  const int Ls=params.Ls;

  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);

  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  ///////////////////////////////////////////////////
  // Construct a coarsened grid; utility for this?
  ///////////////////////////////////////////////////
  std::vector<int> block ({2,2,2,2});
  const int nbasis= 32;
  auto clatt = GridDefaultLatt();
  for(int d=0;d<clatt.size();d++){
    clatt[d] = clatt[d]/block[d];
  }


  GridCartesian *Coarse4d =  SpaceTimeGrid::makeFourDimGrid(clatt, GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());;
  GridCartesian *Coarse5d =  SpaceTimeGrid::makeFiveDimGrid(1,Coarse4d);

  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});
  std::vector<int> cseeds({5,6,7,8});
  GridParallelRNG          RNG5(FGrid);   RNG5.SeedFixedIntegers(seeds5);
  GridParallelRNG          RNG4(UGrid);   RNG4.SeedFixedIntegers(seeds4);
  GridParallelRNG          CRNG(Coarse5d);CRNG.SeedFixedIntegers(cseeds);

  Gamma g5(Gamma::Algebra::Gamma5);

  LatticeFermion    src(FGrid); gaussian(RNG5,src);// src=src+g5*src;
  LatticeFermion result(FGrid); result=Zero();
  LatticeFermion    ref(FGrid); ref=Zero();
  LatticeFermion    tmp(FGrid);
  LatticeFermion    err(FGrid);
  LatticeGaugeField Umu(UGrid); 
  LatticeGaugeField UmuDD(UGrid); 
  LatticeColourMatrix U(UGrid);
  LatticeColourMatrix zz(UGrid);

  FieldMetaData header;
  std::string file("./ckpoint_lat.4000");
  NerscIO::readConfiguration(Umu,header,file);


  if ( params.domaindecompose ) { 
    Lattice<iScalar<vInteger> > coor(UGrid);
    zz=Zero();
    for(int mu=0;mu<Nd;mu++){
      LatticeCoordinate(coor,mu);
      U = PeekIndex<LorentzIndex>(Umu,mu);
      U = where(mod(coor,params.domainsize)==(Integer)0,zz,U);
      PokeIndex<LorentzIndex>(UmuDD,U,mu);
    }
  } else { 
    UmuDD = Umu;
  }
  //  SU3::ColdConfiguration(RNG4,Umu);
  //  SU3::TepidConfiguration(RNG4,Umu);
  //  SU3::HotConfiguration(RNG4,Umu);
  //  Umu=Zero();

  RealD mass=params.mq;
  RealD M5=1.8;

  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  std::cout<<GridLogMessage << "Building g5R5 hermitian DWF operator" <<std::endl;
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  DomainWallFermionR Ddwf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);
  DomainWallFermionR DdwfDD(UmuDD,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);

  typedef Aggregation<vSpinColourVector,vTComplex,nbasis>              Subspace;
  typedef CoarsenedMatrix<vSpinColourVector,vTComplex,nbasis>          CoarseOperator;
  typedef CoarseOperator::CoarseVector                                 CoarseVector;
  typedef CoarseOperator::siteVector siteVector;
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  std::cout<<GridLogMessage << "Calling Aggregation class to build subspace" <<std::endl;
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  MdagMLinearOperator<DomainWallFermionR,LatticeFermion> HermDefOp(Ddwf);

  Subspace Aggregates(Coarse5d,FGrid,0);

  assert ( (nbasis & 0x1)==0);

  {
    int nb=nbasis/2;
    std::cout<<GridLogMessage << " nbasis/2 = "<<nb<<std::endl;

   Aggregates.CreateSubspaceChebyshev(RNG5,HermDefOp,nb,60.0,0.02,500,110);
    for(int n=0;n<nb;n++){
      G5R5(Aggregates.subspace[n+nb],Aggregates.subspace[n]);
    }
  }

  
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  std::cout<<GridLogMessage << "Building coarse representation of Indef operator" <<std::endl;
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  Gamma5R5HermitianLinearOperator<DomainWallFermionR,LatticeFermion> HermIndefOp(Ddwf);
  Gamma5R5HermitianLinearOperator<DomainWallFermionR,LatticeFermion> HermIndefOpDD(DdwfDD);
  CoarsenedMatrix<vSpinColourVector,vTComplex,nbasis> LDOp(*Coarse5d,1); // Hermitian matrix
  LDOp.CoarsenOperator(FGrid,HermIndefOp,Aggregates);
  exit(0);

  CoarseVector c_src (Coarse5d);
  CoarseVector c_res (Coarse5d);
  gaussian(CRNG,c_src);
  result=Zero();
  c_res=Zero();

  //////////////////////////////////////////////////
  // Deflate the course space. Recursive multigrid?
  //////////////////////////////////////////////////

  typedef CoarsenedMatrix<vSpinColourVector,vTComplex,nbasis> Level1Op;
  typedef CoarsenedMatrix<siteVector,iScalar<vTComplex>,nbasis> Level2Op;

  auto cclatt = clatt;
  for(int d=0;d<clatt.size();d++){
    cclatt[d] = clatt[d]/block[d];
  }
  GridCartesian *CoarseCoarse4d =  SpaceTimeGrid::makeFourDimGrid(cclatt, GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());;
  GridCartesian *CoarseCoarse5d =  SpaceTimeGrid::makeFiveDimGrid(1,CoarseCoarse4d);

  typedef Aggregation<siteVector,iScalar<vTComplex>,nbasis>                   CoarseSubspace;
  CoarseSubspace CoarseAggregates(CoarseCoarse5d,Coarse5d,0);

  double c_first = 0.2;
  double c_div   = 1.2;
  std::vector<double> c_lo(nbasis/2);
  c_lo[0] = c_first;
  for(int b=1;b<nbasis/2;b++) {
    c_lo[b] = c_lo[b-1]/c_div;
  }
  std::vector<int> c_ord(nbasis/2,200);
  c_ord[0]=500;

#define RECURSIVE_MULTIGRID
#ifdef RECURSIVE_MULTIGRID
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  std::cout<<GridLogMessage << "Build deflation space in coarse operator "<< std::endl;
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;

  MdagMLinearOperator<CoarseOperator,CoarseVector> PosdefLdop(LDOp);
  //  CoarseAggregates.CreateSubspaceChebyshev(CRNG,PosdefLdop,nbasis,14.0,c_lo,c_ord);
  //  CoarseAggregates.CreateSubspaceRandom(CRNG);

  //  Level2Op L2Op(*CoarseCoarse5d,1); // Hermitian matrix
  //  HermitianLinearOperator<Level1Op,CoarseVector> L1LinOp(LDOp);
  //  L2Op.CoarsenOperator(Coarse5d,L1LinOp,CoarseAggregates);
#endif


  //  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  //  std::cout<<GridLogMessage << "Unprec CG "<< std::endl;
  //  std::cout<<GridLogMessage << "**************************************************"<< std::endl;

  //  TrivialPrecon<LatticeFermion> simple;
  //  ConjugateGradient<LatticeFermion> fCG(1.0e-8,100000);
  //  fCG(HermDefOp,src,result);

    std::cout<<GridLogMessage << "**************************************************"<< std::endl;
    std::cout<<GridLogMessage << "Red Black Prec CG "<< std::endl;
    std::cout<<GridLogMessage << "**************************************************"<< std::endl;
    LatticeFermion    src_o(FrbGrid);
    LatticeFermion result_o(FrbGrid);
    pickCheckerboard(Odd,src_o,src);
    result_o=Zero();
    SchurDiagMooeeOperator<DomainWallFermionR,LatticeFermion> HermOpEO(Ddwf);
    ConjugateGradient<LatticeFermion> pCG(1.0e-8,10000);
    //    pCG(HermOpEO,src_o,result_o);
  
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  std::cout<<GridLogMessage << " Running coarse grid Lanczos "<< std::endl;
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  MdagMLinearOperator<Level1Op,CoarseVector> IRLHermOp(LDOp);
  Chebyshev<CoarseVector> IRLCheby(0.005,16.0,51);
  //  IRLCheby.InitLowPass(0.01,18.0,51);
  FunctionHermOp<CoarseVector> IRLOpCheby(IRLCheby,IRLHermOp);
     PlainHermOp<CoarseVector> IRLOp    (IRLHermOp);

     int Nstop=24;
     int Nk=24;
  int Nm=48;
  ImplicitlyRestartedLanczos<CoarseVector> IRL(IRLOpCheby,IRLOp,Nstop,Nk,Nm,1.0e-3,20);
  int Nconv;
  std::vector<RealD>          eval(Nm);
  std::vector<CoarseVector>   evec(Nm,Coarse5d);
  IRL.calc(eval,evec,c_src,Nconv);


  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  std::cout<<GridLogMessage << "coarse grid CG "<< std::endl;
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;

  //  ConjugateGradient<CoarseVector> CG(3.0e-3,100000);
  //  CG(PosdefLdop,c_src,c_res);

  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  std::cout<<GridLogMessage << "coarse grid Deflated CG with "<< eval.size() << " evecs" << std::endl;
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  
  c_res=Zero();
  DeflatedGuesser<CoarseVector> DeflCoarseGuesser(evec,eval);
  DeflCoarseGuesser(c_src,c_res);
  //  CG(PosdefLdop,c_src,c_res);

  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  std::cout<<GridLogMessage <<" Applying Fine power method to find spectral range      "<<std::endl;
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;

  ZeroGuesser<CoarseVector> CoarseZeroGuesser;

  MultiGridPreconditioner <vSpinColourVector,vTComplex,nbasis,DomainWallFermionR,
			   ZeroGuesser<CoarseVector> >
    Precon  (Aggregates, LDOp,
	     HermIndefOp,Ddwf,
	     HermIndefOp,Ddwf,
	     CoarseZeroGuesser,
	     params);

  //  Precon.PowerMethod(src);
  
  /*  
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  std::cout<<GridLogMessage <<" Applying Coarse power method to find spectral range      "<<std::endl;
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  
  cparams = params;
  cparams.hi   = 20.0;
  cparams.lo   =  0.2;
  cparams.order=  20;

  MultiGridPreconditioner <siteVector,iScalar<vTComplex>,nbasis,Level1Op,ZeroGuesser<CoarseVector> > 
  CoarsePrecon (CoarseAggregates, 
		L2Op,
		L1LinOp,LDOp,
		L1LinOp,LDOp,
		CoarseZeroGuesser,
		cparams);
  
  CoarsePrecon.PowerMethod(c_src);
  */

  /*
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  std::cout<<GridLogMessage << "Building a two level PGCR "<< std::endl;
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  PrecGeneralisedConjugateResidual<LatticeFermion> PGCR(1.0e-8,100000,Precon,8,8);
  std::cout<<GridLogMessage<<"checking norm src "<<norm2(src)<<std::endl;
  result=Zero();
  PGCR(HermIndefOp,src,result);
  */
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  std::cout<<GridLogMessage << "Building a two level deflated PGCR "<< std::endl;
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;

  MultiGridPreconditioner <vSpinColourVector,vTComplex,nbasis,DomainWallFermionR, DeflatedGuesser<CoarseVector> >
    DeflatedPrecon  (Aggregates, LDOp,
		     HermIndefOp,Ddwf,
		     HermIndefOp,Ddwf,
		     DeflCoarseGuesser,
		     params);

  PrecGeneralisedConjugateResidual<LatticeFermion> deflPGCR(1.0e-8,100000,DeflatedPrecon,16,16);

  std::cout<<GridLogMessage<<"checking norm src "<<norm2(src)<<std::endl;
  result=Zero();
  deflPGCR(HermIndefOp,src,result);


  /*
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  std::cout<<GridLogMessage << "Building deflation preconditioner "<< std::endl;
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;

  PrecGeneralisedConjugateResidual<CoarseVector> CPGCR(1.0e-3,10000,CoarsePrecon,8,8);
  std::cout<<GridLogMessage<<"checking norm src "<<norm2(c_src)<<std::endl;
  c_res=Zero();
  CPGCR(L1LinOp,c_src,c_res);
  */

  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  std::cout<<GridLogMessage << "Done "<< std::endl;
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  Grid_finalize();
}
