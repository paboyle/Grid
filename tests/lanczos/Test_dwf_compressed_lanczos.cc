/*
  Authors: Christoph Lehner
  Date: 2017

  Multigrid Lanczos



  TODO:

  High priority:
  - Explore filtering of starting vector again, should really work:  If cheby has 4 for low mode region and 1 for high mode, applying 15 iterations has 1e9 suppression
    of high modes, which should create the desired invariant subspace already?  Missing something here???  Maybe dynamic range dangerous, i.e., could also kill interesting
    eigenrange if not careful.

    Better: Use all Cheby up to order N in order to approximate a step function; try this!  Problem: width of step function.  Can kill eigenspace > 1e-3 and have < 1e-5 equal
            to 1

  Low priority:
  - Given that I seem to need many restarts and high degree poly to create the base and this takes about 1 day, seriously consider a simple method to create a basis
    (ortho krylov low poly); and then fix up lowest say 200 eigenvalues by 1 run with high-degree poly (600 could be enough)
*/
#include <Grid/Grid.h>
#include <Grid/algorithms/iterative/ImplicitlyRestartedLanczos.h>
/////////////////////////////////////////////////////////////////////////////
// The following are now decoupled from the Lanczos and deal with grids.
// Safe to replace functionality
/////////////////////////////////////////////////////////////////////////////
#include "BlockedGrid.h"
#include "FieldBasisVector.h"
#include "BlockProjector.h"
#include "FieldVectorIO.h"
#include "Params.h"

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

bool read_evals(GridBase* _grid, char* fn, std::vector<RealD>& evals) {

  FILE* f = 0;
  uint32_t status = 0;
  if (_grid->IsBoss()) {
    f = fopen(fn,"rt");
    status = f ? 1 : 0;
  }
  _grid->GlobalSum(status);

  if (!status)
    return false;

  uint32_t N;
  if (f)
    assert(fscanf(f,"%d\n",&N)==1);
  else
    N = 0;
  _grid->GlobalSum(N);

  std::cout << "Reading " << N << " eigenvalues" << std::endl;

  evals.resize(N);

  for (int i=0;i<N;i++) {
    if (f)
      assert(fscanf(f,"%lf",&evals[i])==1);
    else
      evals[i] = 0;
  }

  _grid->GlobalSumVector(&evals[0],evals.size());

  if (f)
    fclose(f);
  return true;
}

void write_evals(char* fn, std::vector<RealD>& evals) {
  FILE* f = fopen(fn,"wt");
  assert(f);

  int N = (int)evals.size();
  fprintf(f,"%d\n",N);

  for (int i=0;i<N;i++) {
    fprintf(f,"%.15E\n",evals[i]);
  }

  fclose(f);
}

void write_history(char* fn, std::vector<RealD>& hist) {
  FILE* f = fopen(fn,"wt");
  assert(f);

  int N = (int)hist.size();
  for (int i=0;i<N;i++) {
    fprintf(f,"%d %.15E\n",i,hist[i]);
  }

  fclose(f);
}


template<typename Field>
class CheckpointedLinearFunction : public LinearFunction<Field> {
public:
  LinearFunction<Field>& _op;
  std::string _dir;
  int _max_apply;
  int _apply, _apply_actual;
  GridBase* _grid;
  FILE* _f;

  CheckpointedLinearFunction(GridBase* grid, LinearFunction<Field>& op, const char* dir,int max_apply) : _op(op), _dir(dir), _grid(grid), _f(0),
													 _max_apply(max_apply), _apply(0), _apply_actual(0) {

    FieldVectorIO::conditionalMkDir(dir);

    char fn[4096];
    sprintf(fn,"%s/ckpt_op.%4.4d",_dir.c_str(),_grid->ThisRank());
    printf("CheckpointLinearFunction:: file %s\n",fn);
    _f = fopen(fn,"r+b");
    if (!_f)
      _f = fopen(fn,"w+b");
    assert(_f);
    fseek(_f,0,SEEK_CUR);

  }

  ~CheckpointedLinearFunction() {
    if (_f) {
      fclose(_f);
      _f = 0;
    }
  }

  bool load_ckpt(const Field& in, Field& out) {

    off_t cur = ftello(_f);
    fseeko(_f,0,SEEK_END);
    if (cur == ftello(_f))
      return false;
    fseeko(_f,cur,SEEK_SET);

    size_t sz = sizeof(out._odata[0]) * out._odata.size();

    GridStopWatch gsw;
    gsw.Start();
    uint32_t crc_exp;
    assert(fread(&crc_exp,4,1,_f)==1);
    assert(fread(&out._odata[0],sz,1,_f)==1);
    assert(FieldVectorIO::crc32_threaded((unsigned char*)&out._odata[0],sz,0x0)==crc_exp);
    gsw.Stop();

    printf("CheckpointLinearFunction:: reading %lld\n",(long long)sz);
    std::cout << GridLogMessage << "Loading " << ((RealD)sz/1024./1024./1024.) << " GB in " << gsw.Elapsed() << std::endl;
    return true;
  }

  void save_ckpt(const Field& in, Field& out) {

    fseek(_f,0,SEEK_CUR); // switch to write

    size_t sz = sizeof(out._odata[0]) * out._odata.size();

    GridStopWatch gsw;
    gsw.Start();
    uint32_t crc = FieldVectorIO::crc32_threaded((unsigned char*)&out._odata[0],sz,0x0);
    assert(fwrite(&crc,4,1,_f)==1);
    assert(fwrite(&out._odata[0],sz,1,_f)==1);
    fflush(_f); // try this on the GPFS to suppress OPA usage for disk during dslash; this is not needed at Lustre/JLAB
    gsw.Stop();

    printf("CheckpointLinearFunction:: writing %lld\n",(long long)sz);
    std::cout << GridLogMessage << "Saving " << ((RealD)sz/1024./1024./1024.) << " GB in " << gsw.Elapsed() << std::endl;
  }

  void operator()(const Field& in, Field& out) {

    _apply++;

    if (load_ckpt(in,out))
      return;

    _op(in,out);
    
    save_ckpt(in,out);

    if (_apply_actual++ >= _max_apply) {
      std::cout << GridLogMessage << "Maximum application of operator reached, checkpoint and finish in future job" << std::endl;
      if (_f) { fclose(_f); _f=0; }
      in._grid->Barrier();
      Grid_finalize();
      exit(3);
    }
  }
};

template<typename CoarseField,typename Field>
class ProjectedFunctionHermOp : public LinearFunction<CoarseField> {
public:
  OperatorFunction<Field>   & _poly;
  LinearOperatorBase<Field> &_Linop;
  BlockProjector<Field>& _pr;

  ProjectedFunctionHermOp(BlockProjector<Field>& pr,OperatorFunction<Field> & poly,LinearOperatorBase<Field>& linop) : _poly(poly), _Linop(linop), _pr(pr) {
  }

  void operator()(const CoarseField& in, CoarseField& out) {
    assert(_pr._bgrid._o_blocks == in._grid->oSites());

    Field fin(_pr._bgrid._grid);
    Field fout(_pr._bgrid._grid);

    GridStopWatch gsw1,gsw2,gsw3;
    // fill fin
    gsw1.Start();
    _pr.coarseToFine(in,fin);
    gsw1.Stop();

    // apply poly
    gsw2.Start();
    _poly(_Linop,fin,fout);
    gsw2.Stop();

    // fill out
    gsw3.Start();
    _pr.fineToCoarse(fout,out);
    gsw3.Stop();

    auto eps = innerProduct(in,out);
    std::cout << GridLogMessage << "Operator timing details: c2f = " << gsw1.Elapsed() << " poly = " << gsw2.Elapsed() << " f2c = " << gsw3.Elapsed() << 
      "   Complimentary Hermiticity check: " << eps.imag() / std::abs(eps) << std::endl;

  }
};

template<typename CoarseField,typename Field>
class ProjectedHermOp : public LinearFunction<CoarseField> {
public:
  LinearOperatorBase<Field> &_Linop;
  BlockProjector<Field>& _pr;

  ProjectedHermOp(BlockProjector<Field>& pr,LinearOperatorBase<Field>& linop) : _Linop(linop), _pr(pr) {
  }

  void operator()(const CoarseField& in, CoarseField& out) {
    assert(_pr._bgrid._o_blocks == in._grid->oSites());
    Field fin(_pr._bgrid._grid);
    Field fout(_pr._bgrid._grid);
    _pr.coarseToFine(in,fin);
    _Linop.HermOp(fin,fout);
    _pr.fineToCoarse(fout,out);

  }
};

template<typename vtype, int N > using CoarseSiteFieldGeneral = iScalar< iVector<vtype, N> >;
template<int N> using CoarseSiteFieldD = CoarseSiteFieldGeneral< vComplexD, N >;
template<int N> using CoarseSiteFieldF = CoarseSiteFieldGeneral< vComplexF, N >;
template<int N> using CoarseSiteField  = CoarseSiteFieldGeneral< vComplex,  N >;
template<int N> using CoarseLatticeFermion  = Lattice< CoarseSiteField<N> >;
template<int N> using CoarseLatticeFermionD = Lattice< CoarseSiteFieldD<N> >;

template<typename Field,int Nstop1>
void CoarseGridLanczos(BlockProjector<Field>& pr,RealD alpha2,RealD beta,int Npoly2,
		       int Nstop2,int Nk2,int Nm2,RealD resid2,RealD betastp2,int MaxIt,int MinRes2,
		       LinearOperatorBase<Field>& HermOp, std::vector<RealD>& eval1, bool cg_test_enabled, 
		       int cg_test_maxiter,int nsingle,int SkipTest2, int MaxApply2,bool smoothed_eval_enabled,
		       int smoothed_eval_inner,int smoothed_eval_outer,int smoothed_eval_begin,
		       int smoothed_eval_end,RealD smoothed_eval_inner_resid) {

  BlockedGrid<Field>& bgrid = pr._bgrid;
  BasisFieldVector<Field>& basis = pr._evec;


  std::vector<int> coarseFourDimLatt;
  for (int i=0;i<4;i++)
    coarseFourDimLatt.push_back(bgrid._nb[1+i] * bgrid._grid->_processors[1+i]);
  assert(bgrid._grid->_processors[0] == 1);

  std::cout << GridLogMessage << "CoarseGrid = " << coarseFourDimLatt << " with basis = " << Nstop1 << std::endl;
  GridCartesian         * UCoarseGrid   = SpaceTimeGrid::makeFourDimGrid(coarseFourDimLatt, GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridCartesian         * FCoarseGrid   = SpaceTimeGrid::makeFiveDimGrid(bgrid._nb[0],UCoarseGrid);

  Chebyshev<Field> Cheb2(alpha2,beta,Npoly2);
  CoarseLatticeFermion<Nstop1> src_coarse(FCoarseGrid);

  // Second round of Lanczos in blocked space
  std::vector<RealD>         eval2(Nm2);
  std::vector<RealD>         eval3(Nm2);
  BasisFieldVector<CoarseLatticeFermion<Nstop1> > coef(Nm2,FCoarseGrid);

  ProjectedFunctionHermOp<CoarseLatticeFermion<Nstop1>,LatticeFermion> Op2plain(pr,Cheb2,HermOp);
  CheckpointedLinearFunction<CoarseLatticeFermion<Nstop1> > Op2ckpt(src_coarse._grid,Op2plain,"checkpoint",MaxApply2);
  LinearFunction< CoarseLatticeFermion<Nstop1> >* Op2;
  if (MaxApply2) {
    Op2 = &Op2ckpt;
  } else {
    Op2 = &Op2plain;
  }
  ProjectedHermOp<CoarseLatticeFermion<Nstop1>,LatticeFermion> Op2nopoly(pr,HermOp);
  ImplicitlyRestartedLanczos<CoarseLatticeFermion<Nstop1> > IRL2(*Op2,*Op2,Nstop2,Nk2,Nm2,resid2,MaxIt,betastp2,MinRes2);


  src_coarse = 1.0;
  
  // Precision test
  {
    Field tmp(bgrid._grid);
    CoarseLatticeFermion<Nstop1> tmp2(FCoarseGrid);
    CoarseLatticeFermion<Nstop1> tmp3(FCoarseGrid);
    tmp2 = 1.0;
    tmp3 = 1.0;

    pr.coarseToFine(tmp2,tmp);
    pr.fineToCoarse(tmp,tmp2);

    tmp2 -= tmp3;
    std::cout << GridLogMessage << "Precision Test c->f->c: " << norm2(tmp2) / norm2(tmp3) << std::endl;

    //bgrid._grid->Barrier();
    //return;
  }

  int Nconv;
  if (!FieldVectorIO::read_compressed_vectors("lanczos.output",pr,coef) ||
      !read_evals(UCoarseGrid,(char *)"lanczos.output/eigen-values.txt",eval3) ||
      !read_evals(UCoarseGrid,(char *)"lanczos.output/eigen-values.txt.linear",eval1) ||
      !read_evals(UCoarseGrid,(char *)"lanczos.output/eigen-values.txt.poly",eval2)
      ) {
    

    IRL2.calc(eval2,coef._v,src_coarse,Nconv,true);

    coef.resize(Nstop2);
    eval2.resize(Nstop2);
    eval3.resize(Nstop2);

    std::vector<Field> step3_cache;

    // reconstruct eigenvalues of original operator
    for (int i=0;i<Nstop2;i++){
      RealD eval2_linear;

      if (i<Nstop1) {
	eval2_linear = eval1[i];
      } else {
	eval2_linear = eval2[i-1];
      }

      RealD eval2_poly = eval2[i];
      RealD eval_reconstruct = Cheb2.approxInv(eval2_poly,eval2_linear,100,1e-10);
      std::cout << i << " Reconstructed eval = " << eval_reconstruct << " from quess " << eval2_linear << std::endl;
      eval2[i] = eval_reconstruct;
    }
    
    // as demonstrated in CG test below, best result from mixed determination
    for (int i=0;i<Nstop2;i++)
      eval3[i] = (i < Nstop1) ? eval1[i] : eval2[i];
    
    for(int i=0;i<Nstop2;i++){
      std::cout << i<<" / "<< Nstop2<< " eigenvalue "<< eval3[i] <<std::endl;
    };
    
    // write
    mkdir("lanczos.output",ACCESSPERMS);
    FieldVectorIO::write_compressed_vectors("lanczos.output",pr,coef,nsingle);
    if (bgrid._grid->IsBoss()) {
      write_evals((char *)"lanczos.output/eigen-values.txt",eval3);
      write_evals((char *)"lanczos.output/eigen-values.txt.linear",eval1);
      write_evals((char *)"lanczos.output/eigen-values.txt.poly",eval2);
    }

  }

  // fix up eigenvalues
  if (!read_evals(UCoarseGrid,(char *)"lanczos.output/eigen-values.txt.smoothed",eval3) && smoothed_eval_enabled) {

    ConjugateGradient<LatticeFermion> CG(smoothed_eval_inner_resid, smoothed_eval_inner, false);

    LatticeFermion v_i(basis[0]._grid);
    auto tmp = v_i;
    auto tmp2 = v_i;

    for (int i=smoothed_eval_begin;i<smoothed_eval_end;i++) {

      GridStopWatch gsw;

      gsw.Start();

      pr.coarseToFine(coef[i],v_i);
      v_i.checkerboard = Odd;
      
      for (int j=0;j<smoothed_eval_outer;j++) {
	tmp=zero;
	//pr.deflate(coef,eval3,Nstop2,v_i,tmp);
	CG(HermOp, v_i, tmp);

	v_i = 1.0 / ::sqrt( norm2(tmp) ) * tmp;
      }

      tmp = v_i;

      HermOp.HermOp(tmp,tmp2);

      RealD ev = innerProduct(tmp,tmp2).real();

      gsw.Stop();

      std::cout << GridLogMessage << "Smoothed eigenvalue " << i << " from " << eval3[i] << " to " << ev << " in " << gsw.Elapsed() << std::endl;
      //	" with effective smoother precision " << (CG.ResHistory.back() / CG.ResHistory.front() ) << std::endl;
      //      CG.ResHistory.clear();

      eval3[i] = ev;
    }

    if (bgrid._grid->IsBoss()) {
      write_evals((char *)"lanczos.output/eigen-values.txt.smoothed",eval3);
      write_evals((char *)"lanczos.output/eigen-values.txt",eval3); // also reset this to the best ones we have available
    }
  }

  // do CG test with and without deflation
  if (cg_test_enabled) {
    ConjugateGradient<LatticeFermion> CG(1.0e-8, cg_test_maxiter, false);
    LatticeFermion src_orig(bgrid._grid);
    src_orig.checkerboard = Odd;
    src_orig = 1.0;
    src_orig = src_orig * (1.0 / ::sqrt(norm2(src_orig)) );
    auto result = src_orig; 

    // undeflated solve
    std::cout << GridLogMessage << " Undeflated solve "<<std::endl;
    result = zero;
    CG(HermOp, src_orig, result);
    //    if (UCoarseGrid->IsBoss())
    //      write_history("cg_test.undefl",CG.ResHistory);
    //    CG.ResHistory.clear();

    // deflated solve with all eigenvectors
    std::cout << GridLogMessage << " Deflated solve with all evectors"<<std::endl;
    result = zero;
    pr.deflate(coef,eval2,Nstop2,src_orig,result);
    CG(HermOp, src_orig, result);
    //    if (UCoarseGrid->IsBoss())
    //      write_history("cg_test.defl_all",CG.ResHistory);
    //    CG.ResHistory.clear();

    // deflated solve with non-blocked eigenvectors
    std::cout << GridLogMessage << " Deflated solve with non-blocked evectors"<<std::endl;
    result = zero;
    pr.deflate(coef,eval1,Nstop1,src_orig,result);
    CG(HermOp, src_orig, result);
    //    if (UCoarseGrid->IsBoss())
    //      write_history("cg_test.defl_full",CG.ResHistory);
    //    CG.ResHistory.clear();

    // deflated solve with all eigenvectors and original eigenvalues from proj
    std::cout << GridLogMessage << " Deflated solve with all eigenvectors and original eigenvalues from proj"<<std::endl;
    result = zero;
    pr.deflate(coef,eval3,Nstop2,src_orig,result);
    CG(HermOp, src_orig, result);
    //    if (UCoarseGrid->IsBoss())
    //      write_history("cg_test.defl_all_ev3",CG.ResHistory);
    //    CG.ResHistory.clear();

  }
  
}


template<typename Field>
void quick_krylov_basis(BasisFieldVector<Field>& evec,Field& src,LinearFunction<Field>& Op,int Nstop) {
  Field tmp = src;
  Field tmp2 = tmp;

  for (int i=0;i<Nstop;i++) {
    GridStopWatch gsw;
    gsw.Start();
    Op(tmp,tmp2);
    gsw.Stop();
    evec.orthogonalize(tmp2,i);

    RealD nn = norm2(tmp2);
    nn = Grid::sqrt(nn);
    tmp2 = tmp2 * (1.0/nn);

    evec[i] = tmp2;
    tmp = tmp2;
    std::cout << GridLogMessage << "Quick_krylov_basis: " << i << "/" << Nstop << " timing of operator=" << gsw.Elapsed() << std::endl;
  }

}



int main (int argc, char ** argv) {

  Grid_init(&argc,&argv);

  const int MaxIt = 10000;

  int Ls;
  RealD mass;
  RealD M5;
  std::vector < std::complex<double>  > omega;
  
  RealD alpha1, alpha2, beta;
  int Npoly1, Npoly2;
  int Nstop1, Nstop2;
  int Nk1, Nk2;
  int Np1, Np2;
  int MinRes1, MinRes2;
  int SkipTest2, MaxApply2;
  bool checkpoint_basis;
  bool cg_test_enabled;
  bool exit_after_basis_calculation;
  bool simple_krylov_basis;
  int cg_test_maxiter;
  int nsingle; // store in single precision, the rest in FP16
  int max_cheb_time_ms;
  bool smoothed_eval_enabled;
  int smoothed_eval_inner;
  int smoothed_eval_outer;
  int smoothed_eval_begin;
  int smoothed_eval_end;
  RealD smoothed_eval_inner_resid;

  // vector representation
  std::vector<int> block_size; // 5d block size

  RealD resid1, resid2, betastp1, betastp2, basis_norm_threshold;

  std::string config;
  
  Params jp("params.txt");
  PADD(jp,Npoly1); PADD(jp,Npoly2);
  PADD(jp,max_cheb_time_ms);
  PADD(jp,Nstop1); PADD(jp,Nstop2); PADD(jp,MaxApply2);
  PADD(jp,Nk1); PADD(jp,Nk2); PADD(jp,betastp1); PADD(jp,betastp2);
  PADD(jp,Np1); PADD(jp,Np2); basis_norm_threshold = 1e-5; //PADD(jp,basis_norm_threshold);
  PADD(jp,block_size); PADD(jp,smoothed_eval_enabled); PADD(jp,smoothed_eval_inner);
  PADD(jp,resid1); PADD(jp,resid2); PADD(jp,smoothed_eval_outer);
  PADD(jp,alpha1); PADD(jp,alpha2); PADD(jp,smoothed_eval_begin);
  PADD(jp,MinRes1); PADD(jp,MinRes2); PADD(jp,smoothed_eval_end);
  PADD(jp,beta); PADD(jp,mass); PADD(jp,smoothed_eval_inner_resid);
  PADD(jp,omega); PADD(jp,config); 
  PADD(jp,M5); PADD(jp,cg_test_enabled);
  PADD(jp,cg_test_maxiter); PADD(jp,checkpoint_basis);
  PADD(jp,nsingle); PADD(jp,exit_after_basis_calculation);
  PADD(jp,simple_krylov_basis); PADD(jp,SkipTest2);

  Ls = (int)omega.size();

  // Grids
  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridCartesian         * UGridHP = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplexD::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridRedBlackCartesian * UrbGridHP = SpaceTimeGrid::makeFourDimRedBlackGrid(UGridHP);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridCartesian         * FGridHP   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGridHP);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGridHP = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGridHP);

  // Gauge field
  LatticeGaugeField Umu(UGrid);
  FieldMetaData header;
  NerscIO::readConfiguration(Umu,header,config);
  std::cout << GridLogMessage << "Lattice dimensions: " << GridDefaultLatt()
            << "   Ls: " << Ls << std::endl;

  // ZMobius EO Operator
  ZMobiusFermionR Ddwf(Umu, *FGrid, *FrbGrid, *UGrid, *UrbGrid, mass, M5, omega,1.,0.);
  SchurDiagTwoOperator<ZMobiusFermionR,LatticeFermion> HermOp(Ddwf);

  // Eigenvector storage
  const int Nm1 = Np1 + Nk1;
  const int Nm2 = Np2 + Nk2; // maximum number of vectors we need to keep
  std::cout << GridLogMessage << "Keep " << Nm1 << " full vectors" << std::endl;
  std::cout << GridLogMessage << "Keep " << Nm2 << " total vectors" << std::endl;
  assert(Nm2 >= Nm1);
  BasisFieldVector<LatticeFermion> evec(Nm1,FrbGrid); // start off with keeping full vectors

  // First and second cheby
  Chebyshev<LatticeFermion> Cheb1(alpha1,beta,Npoly1);
  FunctionHermOp<LatticeFermion> Op1(Cheb1,HermOp);
  PlainHermOp<LatticeFermion> Op1test(HermOp);

  // Eigenvalue storage
  std::vector<RealD>          eval1(evec.size());

  // Construct source vector
  LatticeFermion    src(FrbGrid);
  {
    src=1.0;
    src.checkerboard = Odd;

    // normalize
    RealD nn = norm2(src);
    nn = Grid::sqrt(nn);
    src = src * (1.0/nn);
  }

  // Do a benchmark and a quick exit if performance is too little (ugly but needed due to performance fluctuations)
  if (max_cheb_time_ms) {
    // one round of warmup
    auto tmp = src;
    GridStopWatch gsw1,gsw2;
    gsw1.Start();
    Cheb1(HermOp,src,tmp);
    gsw1.Stop();
    Ddwf.ZeroCounters();
    gsw2.Start();
    Cheb1(HermOp,src,tmp);
    gsw2.Stop();
    Ddwf.Report();
    std::cout << GridLogMessage << "Performance check; warmup = " << gsw1.Elapsed() << "  test = " << gsw2.Elapsed() << std::endl;
    int ms = (int)(gsw2.useconds()/1e3);
    if (ms > max_cheb_time_ms) {
      std::cout << GridLogMessage << "Performance too poor: " << ms << " ms, cutoff = " << max_cheb_time_ms << " ms" << std::endl;
      Grid_finalize();
      return 2;
    }

  }

  // First round of Lanczos to get low mode basis
  ImplicitlyRestartedLanczos<LatticeFermion> IRL1(Op1,Op1test,Nstop1,Nk1,Nm1,resid1,MaxIt,betastp1,MinRes1);
  int Nconv;

  char tag[1024];
  if (!FieldVectorIO::read_argonne(evec,(char *)"checkpoint") || !read_evals(UGrid,(char *)"checkpoint/eigen-values.txt",eval1)) {

    if (simple_krylov_basis) {
      quick_krylov_basis(evec,src,Op1,Nstop1);
    } else {
      IRL1.calc(eval1,evec._v,src,Nconv,false);
    }
    evec.resize(Nstop1); // and throw away superfluous
    eval1.resize(Nstop1);
    if (checkpoint_basis)
      FieldVectorIO::write_argonne(evec,(char *)"checkpoint");
    if (UGrid->IsBoss() && checkpoint_basis)
      write_evals((char *)"checkpoint/eigen-values.txt",eval1);

    Ddwf.Report();

    if (exit_after_basis_calculation) {
      Grid_finalize();
      return 0;
    }
  }

  // now test eigenvectors
  if (!simple_krylov_basis) {
    for (int i=0;i<Nstop1;i++){
      auto B = evec[i];
      auto tmp = B;
      auto v = B;
      
      {
	HermOp.HermOp(B,v);
	
	RealD vnum = real(innerProduct(B,v)); // HermOp.
	RealD vden = norm2(B);
	RealD vv0 = norm2(v);
	RealD eval2 = vnum/vden;
	v -= eval2*B;
	RealD vv = norm2(v);
	
	std::cout << i << " OP eval = " << eval2 << " (" << eval1[i] << ") "
		  << "res2 = " << vv << " norm2 = " << norm2(B) << std::endl;
      }
    }
  }

  // do second step only if needed
  if (Nstop1 <= Nstop2) {
    
    // Now setup blocking
    assert(evec.size() == Nstop1);
    BlockedGrid<LatticeFermion> bgrid(FrbGrid, block_size);
    BlockProjector<LatticeFermion> pr(evec,bgrid);
    pr.createOrthonormalBasis(basis_norm_threshold);
    pr.createOrthonormalBasis(basis_norm_threshold); // another round due to precision issues created by local coherence

    constexpr int common_basis_sizes[] = { 60, 250, 400 };
    constexpr int n_common_basis_sizes = sizeof(common_basis_sizes) / sizeof(common_basis_sizes[0]);
    switch (Nstop1) {
#define BASIS(n) case common_basis_sizes[n]:\
      CoarseGridLanczos<LatticeFermion,common_basis_sizes[n]>\
	(pr,alpha2,beta,Npoly2,Nstop2,Nk2,Nm2,resid2,betastp2,MaxIt,MinRes2,HermOp,eval1, \
	 cg_test_enabled,cg_test_maxiter,nsingle,SkipTest2, \
	 MaxApply2,smoothed_eval_enabled,smoothed_eval_inner,smoothed_eval_outer, \
	 smoothed_eval_begin,smoothed_eval_end,smoothed_eval_inner_resid); break;
      BASIS(0);
      BASIS(1);
      BASIS(2);
    default:
      std::cout << GridLogMessage << "Basis size " << Nstop1 << " must be added at compile-time" << std::endl;
      std::cout << GridLogMessage << "Currently available sizes: " << std::endl;
      for (int i=0;i<n_common_basis_sizes;i++) {
	std::cout << GridLogMessage << "  " << common_basis_sizes[i] << std::endl;
      }
    }

  }
    
  Grid_finalize();
}

