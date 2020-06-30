/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_dwf_hdcr.cc

    Copyright (C) 2015

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
#include <Grid/Grid.h>
#include <Grid/algorithms/iterative/PrecGeneralisedConjugateResidual.h>
//#include <Grid/algorithms/iterative/PrecGeneralisedConjugateResidualNonHermitian.h>
#include <Grid/algorithms/iterative/BiCGSTAB.h>

using namespace std;
using namespace Grid;

// TODO
//
// Coarse Grid axpby_ssp_pminus // Inherit from spProj5pm
// Coarse Grid axpby_ssp_pplus

template<class Field>
class CayleyBase
{
public:
  int Ls;
  //    protected:
  RealD mass;
  // Save arguments to SetCoefficientsInternal
  Vector<Coeff_t> _gamma;
  RealD                _zolo_hi;
  RealD                _b;
  RealD                _c;

  // Cayley form Moebius (tanh and zolotarev)
  Vector<Coeff_t> omega;
  Vector<Coeff_t> bs;    // S dependent coeffs
  Vector<Coeff_t> cs;
  Vector<Coeff_t> as;
  // For preconditioning Cayley form
  Vector<Coeff_t> bee;
  Vector<Coeff_t> cee;
  Vector<Coeff_t> aee;
  Vector<Coeff_t> beo;
  Vector<Coeff_t> ceo;
  Vector<Coeff_t> aeo;
  // LDU factorisation of the eeoo matrix
  Vector<Coeff_t> lee;
  Vector<Coeff_t> leem;
  Vector<Coeff_t> uee;
  Vector<Coeff_t> ueem;
  Vector<Coeff_t> dee;
public:
  /////////////////////////////////////////////////////////
  // Replicates functionality
  // Use a common base class approach
  /////////////////////////////////////////////////////////
  // Tanh
  void SetCoefficientsTanh(Approx::zolotarev_data *zdata,RealD b,RealD c)
  {
    Vector<Coeff_t> gamma(this->Ls);
    for(int s=0;s<this->Ls;s++) gamma[s] = zdata->gamma[s];
    SetCoefficientsInternal(1.0,gamma,b,c);
  }
  //Zolo
  void SetCoefficientsZolotarev(RealD zolo_hi,Approx::zolotarev_data *zdata,RealD b,RealD c)
  {
    Vector<Coeff_t> gamma(this->Ls);
    for(int s=0;s<this->Ls;s++) gamma[s] = zdata->gamma[s];
    SetCoefficientsInternal(zolo_hi,gamma,b,c);
  }
  //Zolo
  void SetCoefficientsInternal(RealD zolo_hi,Vector<Coeff_t> & gamma,RealD b,RealD c)
  {
    int Ls=this->Ls;

    ///////////////////////////////////////////////////////////
    // The Cayley coeffs (unprec)
    ///////////////////////////////////////////////////////////
    assert(gamma.size()==Ls);

    omega.resize(Ls);
    bs.resize(Ls);
    cs.resize(Ls);
    as.resize(Ls);
    
    double bpc = b+c;
    double bmc = b-c;
    _b = b;
    _c = c;
    _gamma  = gamma; // Save the parameters so we can change mass later.
    _zolo_hi= zolo_hi;
    for(int i=0; i < Ls; i++){
      as[i] = 1.0;
      omega[i] = _gamma[i]*_zolo_hi; //NB reciprocal relative to Chroma NEF code
      assert(omega[i]!=Coeff_t(0.0));
      bs[i] = 0.5*(bpc/omega[i] + bmc);
      cs[i] = 0.5*(bpc/omega[i] - bmc);
    }

    ////////////////////////////////////////////////////////
    // Constants for the preconditioned matrix Cayley form
    ////////////////////////////////////////////////////////
    bee.resize(Ls);
    cee.resize(Ls);
    beo.resize(Ls);
    ceo.resize(Ls);
    
    for(int i=0;i<Ls;i++){
      bee[i]=as[i]*(bs[i]*(4.0-this->M5) +1.0);     
      assert(bee[i]!=Coeff_t(0.0));
      cee[i]=as[i]*(1.0-cs[i]*(4.0-this->M5));
      beo[i]=as[i]*bs[i];
      ceo[i]=-as[i]*cs[i];
    }
    aee.resize(Ls);
    aeo.resize(Ls);
    for(int i=0;i<Ls;i++){
      aee[i]=cee[i];
      aeo[i]=ceo[i];
    }
    
    //////////////////////////////////////////
    // LDU decomposition of eeoo
    //////////////////////////////////////////
    dee.resize(Ls);
    lee.resize(Ls);
    leem.resize(Ls);
    uee.resize(Ls);
    ueem.resize(Ls);
  
    for(int i=0;i<Ls;i++){
      
      dee[i] = bee[i];
      
      if ( i < Ls-1 ) {
	
	assert(bee[i]!=Coeff_t(0.0));
	assert(bee[0]!=Coeff_t(0.0));
      
	lee[i] =-cee[i+1]/bee[i]; // sub-diag entry on the ith column
      
	leem[i]=mass*cee[Ls-1]/bee[0];
	for(int j=0;j<i;j++) {
	  assert(bee[j+1]!=Coeff_t(0.0));
	  leem[i]*= aee[j]/bee[j+1];
	}
      
	uee[i] =-aee[i]/bee[i];   // up-diag entry on the ith row
	
	ueem[i]=mass;
	for(int j=1;j<=i;j++) ueem[i]*= cee[j]/bee[j];
	ueem[i]*= aee[0]/bee[0];
      
      } else { 
	lee[i] =0.0;
	leem[i]=0.0;
	uee[i] =0.0;
	ueem[i]=0.0;
      }
    }
    
    { 
      Coeff_t delta_d=mass*cee[Ls-1];
      for(int j=0;j<Ls-1;j++) {
	assert(bee[j] != Coeff_t(0.0));
	delta_d *= cee[j]/bee[j];
      }
      dee[Ls-1] += delta_d;
    }  
  };

  //////////////////////////////
  // M and Mdag
  //////////////////////////////
  void M    (const Field &psi, Field &chi)
  {
    Field Din(psi.Grid());
    Meooe5D(psi,Din);
    DW(Din,chi);
    axpby(chi,1.0,1.0,chi,psi); 
    M5D(psi,chi);
  }
  void Mdag (const Field &psi, Field &chi)
  {
    Field Din(psi.Grid());
    DWDag(psi,Din); 
    MeooeDag5D(Din,chi);
    M5Ddag(psi,chi);
    axpby (chi,1.0,1.0,chi,psi); 
  }
  /////////////////////////////////
  // P and Pdag - might be needed
  /////////////////////////////////
  void P(const Field &psi, Field &chi)
  {
    int Ls= this->Ls;
    chi=Zero();
    for(int s=0;s<Ls;s++){
      axpby_ssp_pminus(chi,1.0,chi,1.0,psi,s,s);
      axpby_ssp_pplus (chi,1.0,chi,1.0,psi,s,(s+1)%Ls);
    }
  }
  void Pdag(const Field &psi, Field &chi)
  {
    int Ls= this->Ls;
    chi=Zero();
    for(int s=0;s<Ls;s++){
      axpby_ssp_pminus(chi,1.0,chi,1.0,psi,s,s);
      axpby_ssp_pplus (chi,1.0,chi,1.0,psi,s,(s-1+Ls)%Ls);
    }
  }
  ////////////////////////////////////////////////////////
  // Depends: Dw, M5D,  M5Ddag, Meooe5D, MeooeDag5D,
  ////////////////////////////////////////////////////////
  void M5D   (const Field &psi, Field &chi)
  {
    int Ls=this->Ls;
    Vector<Coeff_t> diag (Ls,1.0);
    Vector<Coeff_t> upper(Ls,-1.0); upper[Ls-1]=mass;
    Vector<Coeff_t> lower(Ls,-1.0); lower[0]   =mass;
    M5D(psi,chi,chi,lower,diag,upper);
  }
  void M5Ddag (const Field &psi, Field &chi)
  {
    int Ls=this->Ls;
    Vector<Coeff_t> diag(Ls,1.0);
    Vector<Coeff_t> upper(Ls,-1.0);
    Vector<Coeff_t> lower(Ls,-1.0);
    upper[Ls-1]=-mass*upper[Ls-1];
    lower[0]   =-mass*lower[0];
    M5Ddag(psi,chi,chi,lower,diag,upper);
  }
  void Meooe5D    (const Field &psi, Field &Din)
  {
    int Ls=this->Ls;
    Vector<Coeff_t> diag = bs;
    Vector<Coeff_t> upper= cs;
    Vector<Coeff_t> lower= cs; 
    upper[Ls-1]=-mass*upper[Ls-1];
    lower[0]   =-mass*lower[0];
    M5D(psi,psi,Din,lower,diag,upper);
  }
  void MeooeDag5D    (const Field &psi, Field &Din)
  {
    int Ls=this->Ls;
    Vector<Coeff_t> diag =bs;
    Vector<Coeff_t> upper=cs;
    Vector<Coeff_t> lower=cs; 
    
    for (int s=0;s<Ls;s++){
      if ( s== 0 ) {
	upper[s] = cs[s+1];
	lower[s] =-mass*cs[Ls-1];
      } else if ( s==(Ls-1) ) { 
	upper[s] =-mass*cs[0];
	lower[s] = cs[s-1];
      } else { 
	upper[s] = cs[s+1];
	lower[s] = cs[s-1];
      }
      upper[s] = conjugate(upper[s]);
      lower[s] = conjugate(lower[s]);
      diag[s]  = conjugate(diag[s]);
    }
    M5Ddag(psi,psi,Din,lower,diag,upper);
  }
  
  virtual void M5D(const Field &psi_i,
		   const Field &phi_i, 
		   Field &chi_i,
		   Vector<Coeff_t> &lower,
		   Vector<Coeff_t> &diag,
		   Vector<Coeff_t> &upper);
  virtual void M5Ddag(const Field &psi_i,
		      const Field &phi_i, 
		      Field &chi_i,
		      Vector<Coeff_t> &lower,
		      Vector<Coeff_t> &diag,
		      Vector<Coeff_t> &upper);

};

template<class Fobj,class CComplex,int nbasis>
class CoarseCayleyFermion 
  : public CayleyBase< Lattice<siteVector> >,
    public SparseMatrixBase<Lattice<iVector<CComplex,nbasis > > > 
{
public:
  typedef iVector<CComplex,nbasis >           siteVector;
  typedef Lattice<CComplex >                  CoarseComplexField;
  typedef Lattice<siteVector>                 CoarseVector;
  typedef Lattice<iMatrix<CComplex,nbasis > > CoarseMatrix;
  typedef iMatrix<CComplex,nbasis >  Cobj;
  typedef Lattice< CComplex >   CoarseScalar; // used for inner products on fine field
  typedef Lattice<Fobj >        FineField;

  // Similar to the CoarseOperator but add 5D support.
  Geometry  geom;
  GridBase *Coarse5D;
  GridBase *Coarse4D;
  CartesianStencil<siteVector,siteVector,int> Stencil; 
  CoarsenedMatrix<Fobj,CComplex,nbasis> &Dw;

  GridBase * Grid(void)         { return Coarse5D; };   // this is all the linalg routines need to know

  CoarseCayleyFermion(GridCartesian &CoarseGrid4,
		      GridCartesian &CoarseGrid5
		      CoarsenedMatrix<Fobj,CComplex,nbasis> &_Dw ):
    Coarse4D(&CoarseGrid4),
    Coarse5D(&CoarseGrid5),
    Dw(_Dw),
    geom(CoarseGrid4._ndimension),
    Stencil( &CoarseGrid4,geom.npoint,Even,geom.directions,geom.displacements,0),
    A(2*Nd+1,&CoarseGrid4)
  { 
    Ls=Coarse5D->_fdimensions[0];
    RealD eps = 1.0;
    Approx::zolotarev_data *zdata = Approx::higham(eps,Ls);// eps is ignored for higham
    SetCoefficientsTanh(zdata,1.0,0.0);
    Approx::zolotarev_free(zdata);
  };

public:
  
  ////////////////////////////////////////////////
  // This is specific to Coarse Grid Cayley
  ////////////////////////////////////////////////
  void DW (const CoarseVector &in, CoarseVector &out)
  {
    conformable(Coarse5D,in.Grid());
    conformable(in.Grid(),out.Grid());

    SimpleCompressor<siteVector> compressor;

    Stencil.HaloExchange(in,compressor);
    autoView( in_v , in, AcceleratorRead);
    autoView( out_v , out, AcceleratorWrite);
    typedef LatticeView<Cobj> Aview;
      
    Vector<Aview> AcceleratorViewContainer;
  
    for(int p=0;p<geom.npoint;p++) AcceleratorViewContainer.push_back(Dw.A[p].View(AcceleratorRead));
    Aview *Aview_p = & AcceleratorViewContainer[0];

    const int Nsimd = CComplex::Nsimd();
    typedef decltype(coalescedRead(in_v[0])) calcVector;
    typedef decltype(coalescedRead(in_v[0](0))) calcComplex;

    int osites=Grid()->oSites();
    
    // Ls loop for2D
    int Ls=this->Ls;
    accelerator_for2d(sF, oSites*Ls, b, nbasis, Nsimd, {

      int sU = sF/Ls;
      int  s = sF%Ls;
      calcComplex res = Zero();
      calcVector nbr;
      int ptype;
      StencilEntry *SE;

      for(int point=0;point<geom.npoint;point++){

	SE=Stencil.GetEntry(ptype,point,sF);
	  
	if(SE->_is_local) { 
	  nbr = coalescedReadPermute(in_v[SE->_offset],ptype,SE->_permute);
	} else {
	  nbr = coalescedRead(Stencil.CommBuf()[SE->_offset]);
	}
	acceleratorSynchronise();

	for(int bb=0;bb<nbasis;bb++) {
	  res = res + coalescedRead(Aview_p[point][sU](b,bb))*nbr(bb);
	}
      }
      coalescedWrite(out_v[sF](b),res);
      });

    for(int p=0;p<geom.npoint;p++) AcceleratorViewContainer[p].ViewClose();
  };

  void DWDag (const CoarseVector &in, CoarseVector &out)
  {
    // Inefficient G5 hermitian use
    CoarseVector tmp(Grid());
    G5C(tmp, in); //There has to be a better way
    M(tmp, out);
    G5C(out, out);
  };

  void M5D(const FermionField &psi_i,
	   const FermionField &phi_i, 
	   FermionField &chi_i,
	   Vector<Coeff_t> &lower,
	   Vector<Coeff_t> &diag,
	   Vector<Coeff_t> &upper)
  {
    
    chi_i.Checkerboard()=psi_i.Checkerboard();
    GridBase *grid=psi_i.Grid();
    autoView(psi , psi_i,AcceleratorRead);
    autoView(phi , phi_i,AcceleratorRead);
    autoView(chi , chi_i,AcceleratorWrite);
    assert(phi.Checkerboard() == psi.Checkerboard());

    auto pdiag = &diag[0];
    auto pupper = &upper[0];
    auto plower = &lower[0];

    int Ls =this->Ls;
    
    // 10 = 3 complex mult + 2 complex add
    // Flops = 10.0*(Nc*Ns) *Ls*vol (/2 for red black counting)
    uint64_t nloop = grid->oSites()/Ls;
    accelerator_for(sss,nloop,Simd::Nsimd(),{
	uint64_t ss= sss*Ls;
	typedef decltype(coalescedRead(psi[0])) spinor;
	spinor tmp1, tmp2;
	for(int s=0;s<Ls;s++){
	  uint64_t idx_u = ss+((s+1)%Ls);
	  uint64_t idx_l = ss+((s+Ls-1)%Ls);
	  spProj5m(tmp1,psi(idx_u)); // Need routines for this
	  spProj5p(tmp2,psi(idx_l));
	  coalescedWrite(chi[ss+s],pdiag[s]*phi(ss+s)+pupper[s]*tmp1+plower[s]*tmp2);
	}
      });
  }
  void M5Ddag(const FermionField &psi_i,
	      const FermionField &phi_i, 
	      FermionField &chi_i,
	      Vector<Coeff_t> &lower,
	      Vector<Coeff_t> &diag,
	      Vector<Coeff_t> &upper)
  {
    chi_i.Checkerboard()=psi_i.Checkerboard();
    GridBase *grid=psi_i.Grid();
    autoView(psi , psi_i,AcceleratorRead);
    autoView(phi , phi_i,AcceleratorRead);
    autoView(chi , chi_i,AcceleratorWrite);
    assert(phi.Checkerboard() == psi.Checkerboard());
    
    auto pdiag = &diag[0];
    auto pupper = &upper[0];
    auto plower = &lower[0];
    
    int Ls=this->Ls;
    
    uint64_t nloop = grid->oSites()/Ls;
    accelerator_for(sss,nloop,Simd::Nsimd(),{
	uint64_t ss=sss*Ls;
	typedef decltype(coalescedRead(psi[0])) spinor;
	spinor tmp1,tmp2;
	for(int s=0;s<Ls;s++){
	  uint64_t idx_u = ss+((s+1)%Ls);
	  uint64_t idx_l = ss+((s+Ls-1)%Ls);
	  spProj5p(tmp1,psi(idx_u));
	  spProj5m(tmp2,psi(idx_l));
	  coalescedWrite(chi[ss+s],pdiag[s]*phi(ss+s)+pupper[s]*tmp1+plower[s]*tmp2);
	}
      });
  }


};

template<class Field> class SolverWrapper : public LinearFunction<Field> {
private:
  LinearOperatorBase<Field> & _Matrix;
  OperatorFunction<Field> & _Solver;
  LinearFunction<Field>   & _Guess;
public:

  /////////////////////////////////////////////////////
  // Wrap the usual normal equations trick
  /////////////////////////////////////////////////////
  SolverWrapper(LinearOperatorBase<Field> &Matrix,
	      OperatorFunction<Field> &Solver,
	      LinearFunction<Field> &Guess) 
   :  _Matrix(Matrix), _Solver(Solver), _Guess(Guess) {}; 

  void operator() (const Field &in, Field &out){
 
    _Guess(in,out);
    _Solver(_Matrix,in,out);  // Mdag M out = Mdag in

  }     
};


// Must use a non-hermitian solver
template<class Matrix,class Field>
class PVdagMLinearOperator : public LinearOperatorBase<Field> {
  Matrix &_Mat;
  Matrix &_PV;
public:
  PVdagMLinearOperator(Matrix &Mat,Matrix &PV): _Mat(Mat),_PV(PV){};

  void OpDiag (const Field &in, Field &out) {
    assert(0);
  }
  void OpDir  (const Field &in, Field &out,int dir,int disp) {
    assert(0);
  }
  void OpDirAll  (const Field &in, std::vector<Field> &out){
    assert(0);
  };
  void Op     (const Field &in, Field &out){
    Field tmp(in.Grid());
    _Mat.M(in,tmp);
    _PV.Mdag(tmp,out);
  }
  void AdjOp     (const Field &in, Field &out){
    Field tmp(in.Grid());
    _PV.M(tmp,out);
    _Mat.Mdag(in,tmp);
  }
  void HermOpAndNorm(const Field &in, Field &out,RealD &n1,RealD &n2){
    assert(0);
  }
  void HermOp(const Field &in, Field &out){
    assert(0);
  }
};


RealD InverseApproximation(RealD x){
  return 1.0/x;
}

template<class Field,class Matrix> class ChebyshevSmoother : public LinearFunction<Field>
{
public:
  typedef LinearOperatorBase<Field>                            FineOperator;
  Matrix         & _SmootherMatrix;
  FineOperator   & _SmootherOperator;
  
  Chebyshev<Field> Cheby;

  ChebyshevSmoother(RealD _lo,RealD _hi,int _ord, FineOperator &SmootherOperator,Matrix &SmootherMatrix) :
    _SmootherOperator(SmootherOperator),
    _SmootherMatrix(SmootherMatrix),
    Cheby(_lo,_hi,_ord,InverseApproximation)
  {};

  void operator() (const Field &in, Field &out) 
  {
    Field tmp(in.Grid());
    MdagMLinearOperator<Matrix,Field>   MdagMOp(_SmootherMatrix); 
    _SmootherOperator.AdjOp(in,tmp);
    Cheby(MdagMOp,tmp,out);         
  }
};

template<class Field,class Matrix> class MirsSmoother : public LinearFunction<Field>
{
public:
  typedef LinearOperatorBase<Field>                            FineOperator;
  Matrix         & SmootherMatrix;
  FineOperator   & SmootherOperator;
  RealD tol;
  RealD shift;
  int   maxit;

  MirsSmoother(RealD _shift,RealD _tol,int _maxit,FineOperator &_SmootherOperator,Matrix &_SmootherMatrix) :
    shift(_shift),tol(_tol),maxit(_maxit),
    SmootherOperator(_SmootherOperator),
    SmootherMatrix(_SmootherMatrix)
  {};

  void operator() (const Field &in, Field &out) 
  {
    ZeroGuesser<Field> Guess;
    ConjugateGradient<Field>  CG(tol,maxit,false);
 
    Field src(in.Grid());

    ShiftedMdagMLinearOperator<SparseMatrixBase<Field>,Field> MdagMOp(SmootherMatrix,shift);
    SmootherOperator.AdjOp(in,src);
    Guess(src,out);
    CG(MdagMOp,src,out); 
  }
};

#define GridLogLevel std::cout << GridLogMessage <<std::string(level,'\t')<< " Level "<<level <<" "

template<class Fobj,class CComplex,int nbasis, class CoarseSolver>
class HDCRPreconditioner : public LinearFunction< Lattice<Fobj> > {
public:

  typedef Aggregation<Fobj,CComplex,nbasis> Aggregates;
  typedef CoarsenedMatrix<Fobj,CComplex,nbasis> CoarseOperator;
  typedef typename Aggregation<Fobj,CComplex,nbasis>::CoarseVector CoarseVector;
  typedef typename Aggregation<Fobj,CComplex,nbasis>::CoarseMatrix CoarseMatrix;
  typedef typename Aggregation<Fobj,CComplex,nbasis>::FineField    FineField;
  typedef LinearOperatorBase<FineField>                            FineOperator;
  typedef LinearFunction    <FineField>                            FineSmoother;

  Aggregates     & _Aggregates;
  FineOperator   & _FineOperator;
  FineSmoother   & _Smoother;
  CoarseSolver   & _CoarseSolve;

  int    level;  void Level(int lv) {level = lv; };


  HDCRPreconditioner(Aggregates &Agg,
		     FineOperator &Fine,
		     FineSmoother &Smoother,
		     CoarseSolver &CoarseSolve_)
    : _Aggregates(Agg),
      _FineOperator(Fine),
      _Smoother(Smoother),
      _CoarseSolve(CoarseSolve_),
      level(1)  {  }

  virtual void operator()(const FineField &in, FineField & out) 
  {
    auto CoarseGrid = _Aggregates.CoarseGrid;
    CoarseVector Csrc(CoarseGrid);
    CoarseVector Csol(CoarseGrid);
    FineField vec1(in.Grid());
    FineField vec2(in.Grid());

    double t;
    // Fine Smoother
    t=-usecond();
    _Smoother(in,out);
    t+=usecond();
    GridLogLevel << "Smoother took "<< t/1000.0<< "ms" <<std::endl;

    // Update the residual
    _FineOperator.Op(out,vec1);  sub(vec1, in ,vec1);   

    // Fine to Coarse 
    t=-usecond();
    _Aggregates.ProjectToSubspace  (Csrc,vec1);
    t+=usecond();
    GridLogLevel << "Project to coarse took "<< t/1000.0<< "ms" <<std::endl;

    // Coarse correction
    t=-usecond();
    _CoarseSolve(Csrc,Csol);
    t+=usecond();
    GridLogLevel << "Coarse solve took "<< t/1000.0<< "ms" <<std::endl;

    // Coarse to Fine
    t=-usecond();
    _Aggregates.PromoteFromSubspace(Csol,vec1); 
    add(out,out,vec1);
    t+=usecond();
    GridLogLevel << "Promote to this level took "<< t/1000.0<< "ms" <<std::endl;

    // Residual
    _FineOperator.Op(out,vec1);  sub(vec1 ,in , vec1);  

    // Fine Smoother
    t=-usecond();
    _Smoother(vec1,vec2);
    t+=usecond();
    GridLogLevel << "Smoother took "<< t/1000.0<< "ms" <<std::endl;

    add( out,out,vec2);
  }
};

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  const int Ls=16;

  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);

  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  ///////////////////////////////////////////////////
  // Construct a coarsened grid; utility for this?
  ///////////////////////////////////////////////////
  std::vector<int> block ({2,2,2,2});
  const int nbasis= 8;

  auto clatt = GridDefaultLatt();
  for(int d=0;d<clatt.size();d++){
    clatt[d] = clatt[d]/block[d];
  }

  GridCartesian *Coarse4d =  SpaceTimeGrid::makeFourDimGrid(clatt, GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());;
  GridCartesian *Coarse5d =  SpaceTimeGrid::makeFiveDimGrid(Ls,Coarse4d);

  std::vector<int> seeds({1,2,3,4});
  GridParallelRNG          RNG5(FGrid);   RNG5.SeedFixedIntegers(seeds);
  GridParallelRNG          RNG4(UGrid);   RNG4.SeedFixedIntegers(seeds);
  GridParallelRNG          CRNG(Coarse5d);CRNG.SeedFixedIntegers(seeds);

  LatticeGaugeField Umu(UGrid); 
  FieldMetaData header;
  std::string file("./ckpoint_lat.4000");
  NerscIO::readConfiguration(Umu,header,file);

  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  std::cout<<GridLogMessage << "Building g5R5 hermitian DWF operator" <<std::endl;
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  RealD mass=0.001;
  RealD M5=1.8;
  WilsonFermionR    Dw(Umu,*UGrid,*UrbGrid,-M5);
  DomainWallFermionR Ddwf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);
  DomainWallFermionR Dpv (Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,1.0,M5);

  typedef Aggregation<vSpinColourVector,vTComplex,nbasis>              Subspace;
  typedef CoarsenedMatrix<vSpinColourVector,vTComplex,nbasis>          CoarseOperator;
  typedef CoarseOperator::CoarseVector                                 CoarseVector;
  typedef CoarseOperator::siteVector siteVector;

  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  std::cout<<GridLogMessage << "Calling Aggregation class to build subspace" <<std::endl;
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  MdagMLinearOperator<WilsonFermionR,LatticeFermion> SubspaceOp(Dw);

  Subspace Aggregates4D(Coarse4d,UGrid,0);
  Subspace Aggregates5D(Coarse5d,FGrid,0);

  assert ( (nbasis & 0x1)==0);
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  std::cout<<GridLogMessage << " 4D subspace build                                " <<std::endl;
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  int nb=nbasis/2;
  Gamma g5(Gamma::Algebra::Gamma5);
  Aggregates4D.CreateSubspaceChebyshev(RNG4,SubspaceOp,nb,60.0,0.02,500,100,100,0.0);
  for(int n=0;n<nb;n++){
    Aggregates4D.subspace[n+nb]= Aggregates4D.subspace[n] - g5 * Aggregates4D.subspace[n];
    Aggregates4D.subspace[n]   = Aggregates4D.subspace[n] + g5 * Aggregates4D.subspace[n];
  }

  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  std::cout<<GridLogMessage << " Coarsen the operator                          " <<std::endl;
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  typedef CoarsenedMatrix<vSpinColourVector,vTComplex,nbasis>    Level1Op;

  NonHermitianLinearOperator<DomainWallFermionR,LatticeFermion>  LinOpDwf(Ddwf);

  Level1Op c_Dw  (*Coarse4d,0);   
  
  std::cout<<GridLogMessage << " Callinig Coarsen the operator                          " <<std::endl;
  LDOp.CoarsenOperator(FGrid,LinOpDwf,Aggregates4D);

  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  std::cout<<GridLogMessage << " Solve                                            " <<std::endl;
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  
  LatticeFermion    src(FGrid); gaussian(RNG5,src);
  LatticeFermion result(FGrid); 

  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  std::cout<<GridLogMessage << "Done "<< std::endl;
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  Grid_finalize();
  
}
