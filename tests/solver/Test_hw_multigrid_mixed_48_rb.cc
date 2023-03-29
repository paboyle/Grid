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
#include <Grid/algorithms/iterative/PrecGeneralisedConjugateResidualNonHermitian.h>
#include <Grid/algorithms/iterative/BiCGSTAB.h>

using namespace std;
using namespace Grid;

// TODO
//
// Coarse Grid axpby_ssp_pminus // Inherit from spProj5pm
// Coarse Grid axpby_ssp_pplus

template<class Field,class Coeff_t>
class CayleyBase : public SparseMatrixBase<Field> 
{
public:
  int Ls;
  //    protected:
  RealD mass;
  RealD M5;
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
  CayleyBase(RealD _M5, RealD _mass, int _Ls, RealD b_, RealD c_) :
    M5(_M5),
    mass(_mass),
    Ls(_Ls),
    _b(b_),
    _c(c_)
  {
    RealD eps = 1.0;
    Approx::zolotarev_data *zdata = Approx::higham(eps,this->Ls);// eps is ignored for higham
    this->SetCoefficientsTanh(zdata,1.0,0.0);
    Approx::zolotarev_free(zdata);
  }
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
  virtual  void Mdiag    (const Field &in, Field &out) {assert(0);}
  virtual  void Mdir     (const Field &in, Field &out,int dir, int disp){assert(0);};
  virtual  void MdirAll  (const Field &in, std::vector<Field> &out){assert(0);};
  virtual  void DW       (const Field &psi, Field &chi)=0;
  virtual  void DWDag    (const Field &psi, Field &chi)=0;

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

  void M5D(const Field &psi_i,
	   const Field &phi_i, 
	   Field &chi_i,
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

    const int Nsimd = Field::vector_type::Nsimd();
    accelerator_for(sss,nloop,Nsimd,{
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
  void M5Ddag(const Field &psi_i,
	      const Field &phi_i, 
	      Field &chi_i,
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
    const int Nsimd = Field::vector_type::Nsimd();
    accelerator_for(sss,nloop,Nsimd,{
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

template<class Fobj,class CComplex,int nbasis>
class CoarseCayleyFermion  : public CayleyBase< Lattice<iVector<CComplex,nbasis > > , ComplexD >
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
		      GridCartesian &CoarseGrid5,
		      CoarsenedMatrix<Fobj,CComplex,nbasis> &_Dw,
		      RealD M5, RealD mass, int Ls, RealD b, RealD c) :
    CayleyBase<CoarseVector,ComplexD>(M5,mass,Ls,b,c),
    Coarse4D(&CoarseGrid4),
    Coarse5D(&CoarseGrid5),
    Dw(_Dw),
    geom(CoarseGrid5._ndimension),
    Stencil( &CoarseGrid5,geom.npoint,Even,geom.directions,geom.displacements,0)
  { 
  };

public:
  void Project( CoarseVector &C )
  {
    const int Nsimd = CComplex::Nsimd();
    autoView(Cv,C, AcceleratorWrite);
    int Ls = this->Ls;
    for(int s=0;s<Ls;s++){
      accelerator_for(sU, Coarse4D->oSites(), Nsimd, {
	  int sF= sU*Ls+s;
	  auto tmp = coalescedRead(Cv[sF]);
	  coalescedWrite(Cv[sF],tmp);
      });
    }
  }
  ////////////////////////////////////////////////
  // This is specific to Coarse Grid Cayley
  ////////////////////////////////////////////////
  virtual  void Mdiag    (const CoarseVector &in, CoarseVector &out)
  {
    std::vector<CoarseVector> allout(9,in.Grid());
    this->MdirAll(in,allout);
    out = allout[8];
  }
  virtual  void Mdir     (const CoarseVector &in, CoarseVector &out,int dir, int disp)
  {
    assert(0);
  }
  virtual  void MdirAll  (const CoarseVector &in, std::vector<CoarseVector> &out)
  {
    conformable(Coarse5D,in.Grid());

    SimpleCompressor<siteVector> compressor;

    Stencil.HaloExchange(in,compressor);
    typedef LatticeView<Cobj> Aview;
      
    const int Nsimd = CComplex::Nsimd();
    
    // Ls loop for2D
    int Ls=this->Ls;

    siteVector *CBp=Stencil.CommBuf();			

    //int ptype;
    //    int nb2=nbasis/2;
    
    autoView(in_v ,   in, AcceleratorRead);
    autoView(st, Stencil, AcceleratorRead);
    for(int point=0;point<geom.npoint;point++){
      
      autoView(out_v,  out[point], AcceleratorWrite);
      autoView(Aview,Dw.A[point],AcceleratorRead);

      accelerator_for2d(sF, Coarse5D->oSites(), b, nbasis, Nsimd, {

	  typedef decltype(coalescedRead(in_v[0])) calcVector;
	  typedef decltype(coalescedRead(in_v[0](0))) calcComplex;
	  int sU = sF/Ls;
	  //	  int  s = sF%Ls;

	  calcComplex res = Zero();
	  calcVector  nbr;
	  int ptype;
	    
	  StencilEntry *SE=st.GetEntry(ptype,point,sF);
	  
	  if(SE->_is_local) { 
	    nbr = coalescedReadPermute(in_v[SE->_offset],ptype,SE->_permute);
	  } else {
	    nbr = coalescedRead(CBp[SE->_offset]);
	  }
	  acceleratorSynchronise();

	  for(int bb=0;bb<nbasis;bb++) {
	    res = res + coalescedRead(Aview[sU](b,bb))*nbr(bb);
	  }
	  
	  coalescedWrite(out_v[sF](b),res);
      });
    }      
  }
  void DW (const CoarseVector &in, CoarseVector &out)
  {
    conformable(Coarse5D,in.Grid());
    conformable(in.Grid(),out.Grid());

    SimpleCompressor<siteVector> compressor;

    Stencil.HaloExchange(in,compressor);
    typedef LatticeView<Cobj> Aview;
      
    const int Nsimd = CComplex::Nsimd();
    
    // Ls loop for2D
    int Ls=this->Ls;

    Vector<Aview> AcceleratorViewContainer;
    for(int p=0;p<geom.npoint;p++) AcceleratorViewContainer.push_back(Dw.A[p].View(AcceleratorRead));
    Aview *Aview_p = & AcceleratorViewContainer[0];
    autoView(in_v ,   in, AcceleratorRead);
    autoView(out_v,  out, AcceleratorWrite);
    autoView(st, Stencil, AcceleratorRead);
    siteVector *CBp=Stencil.CommBuf();			

    //    int ptype;
    //    int nb2=nbasis/2;
    accelerator_for2d(sF, Coarse5D->oSites(), b, nbasis, Nsimd, {

      typedef decltype(coalescedRead(in_v[0])) calcVector;
      typedef decltype(coalescedRead(in_v[0](0))) calcComplex;
      int sU = sF/Ls;
      //      int  s = sF%Ls;

      calcComplex res = Zero();

      {
      	calcVector  nbr;
	int ptype;

	for(int point=0;point<geom.npoint;point++){

	  StencilEntry *SE=st.GetEntry(ptype,point,sF);
	  
	  if(SE->_is_local) { 
	    nbr = coalescedReadPermute(in_v[SE->_offset],ptype,SE->_permute);
	  } else {
	    nbr = coalescedRead(CBp[SE->_offset]);
	  }
	  acceleratorSynchronise();

	  for(int bb=0;bb<nbasis;bb++) {
	    res = res + coalescedRead(Aview_p[point][sU](b,bb))*nbr(bb);
	  }	  
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
    DW(tmp, out);
    G5C(out, out);
  };

  typedef Aggregation<Fobj,CComplex,nbasis> Aggregates;

  void PromoteFromSubspace(Aggregates &_Aggregates,CoarseVector &C,FineField &F) 
  {
    auto FineGrid4 = _Aggregates.FineGrid;
    FineField F4(FineGrid4);
    CoarseVector C4(Coarse4D);
    for(int s=0;s<this->Ls;s++){
      ExtractSlice(C4,C,s,0);
      _Aggregates.PromoteFromSubspace(C4,F4); 
      InsertSlice(F4,F,s,0);
    }      
  }
  void ProjectToSubspace(Aggregates &_Aggregates,CoarseVector &C,FineField &F) 
  {
    auto FineGrid4 = _Aggregates.FineGrid;
    FineField F4(FineGrid4);
    CoarseVector C4(Coarse4D);
    for(int s=0;s<this->Ls;s++){
      ExtractSlice(F4,F,s,0);
      _Aggregates.ProjectToSubspace  (C4,F4);
      InsertSlice(C4,C,s,0);
    }
    Project(C);
  }
  template<class Ddwf>
  void Test(Aggregates &_Aggregates,GridBase *FineGrid, Ddwf &_Ddwf)
  {
    typedef Lattice<Fobj> FineField;
    CoarseVector Cin(Coarse5D);
    CoarseVector Cout(Coarse5D);
    CoarseVector CFout(Coarse5D);

    FineField Fin(FineGrid);
    FineField Fout(FineGrid);


    std::vector<int> seeds({1,2,3,4,5});
    GridParallelRNG RNG(Coarse5D);  RNG.SeedFixedIntegers(seeds);

    gaussian(RNG,Cin);
    PromoteFromSubspace(_Aggregates,Cin,Fin);
    ProjectToSubspace(_Aggregates,Cin,Fin);

    std::cout << GridLogMessage<< "************  "<<std::endl;
    std::cout << GridLogMessage<< " Testing M  "<<std::endl;
    std::cout << GridLogMessage<< "************  "<<std::endl;
    // Coarse operator
    this->M(Cin,Cout);
    this->Project(Cout);
    std::cout << GridLogMessage<< " Cout  "<<norm2(Cout)<<std::endl;

    // Fine projected operator
    PromoteFromSubspace(_Aggregates,Cin,Fin);
    _Ddwf.M(Fin,Fout);
    ProjectToSubspace(_Aggregates,CFout,Fout);
    std::cout << GridLogMessage<< " CFout "<<norm2(CFout)<<std::endl;
    CFout = CFout-Cout;
    std::cout << GridLogMessage<< " diff  "<<norm2(CFout)<<std::endl;

    std::cout << GridLogMessage<< "************  "<<std::endl;
    std::cout << GridLogMessage<< " Testing Mdag  "<<std::endl;
    std::cout << GridLogMessage<< "************  "<<std::endl;
    // Coarse operator
    this->Mdag(Cin,Cout);
    this->Project(Cout);
    std::cout << GridLogMessage<< " Cout  "<<norm2(Cout)<<std::endl;

    // Fine operator
    _Ddwf.Mdag(Fin,Fout);
    ProjectToSubspace(_Aggregates,CFout,Fout);
    std::cout << GridLogMessage<< " CFout "<<norm2(CFout)<<std::endl;
    CFout = CFout-Cout;
    std::cout << GridLogMessage<< " diff  "<<norm2(CFout)<<std::endl;
 
  }
  virtual std::vector<int> Directions(void)   { return geom.directions;};
  virtual std::vector<int> Displacements(void){ return geom.displacements;};
};

template<class Field> class SchurSolverWrapper : public LinearFunction<Field> {
private:
  CheckerBoardedSparseMatrixBase<Field> & _Matrix;
  SchurRedBlackBase<Field> & _Solver;
public:
  using LinearFunction<Field>::operator();
  /////////////////////////////////////////////////////
  // Wrap the usual normal equations trick
  /////////////////////////////////////////////////////
  SchurSolverWrapper(CheckerBoardedSparseMatrixBase<Field> &Matrix,
		SchurRedBlackBase<Field> &Solver)
   :  _Matrix(Matrix), _Solver(Solver) {}; 

  void operator() (const Field &in, Field &out){
 
    _Solver(_Matrix,in,out);  // Mdag M out = Mdag in

  }     
};

template<class Field> class SolverWrapper : public LinearFunction<Field> {
private:
  LinearOperatorBase<Field> & _Matrix;
  OperatorFunction<Field> & _Solver;
  LinearFunction<Field>   & _Guess;
public:
  using LinearFunction<Field>::operator();

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

  virtual std::vector<int> Directions(void)   { return _Mat.Directions();};
  virtual std::vector<int> Displacements(void){ return _Mat.Displacements();};

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
  using LinearFunction<Field>::operator();
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
template<class Fobj,class CComplex,int nbasis, class CoarseSolver>
class MGPreconditioner : public LinearFunction< Lattice<Fobj> > {
public:
  using LinearFunction<Lattice<Fobj> >::operator();

  typedef Aggregation<Fobj,CComplex,nbasis> Aggregates;
  typedef typename Aggregation<Fobj,CComplex,nbasis>::CoarseVector CoarseVector;
  typedef typename Aggregation<Fobj,CComplex,nbasis>::CoarseMatrix CoarseMatrix;
  typedef typename Aggregation<Fobj,CComplex,nbasis>::FineField    FineField;
  typedef LinearOperatorBase<FineField>                            FineOperator;
  typedef LinearFunction    <FineField>                            FineSmoother;
  typedef CoarseCayleyFermion<Fobj,CComplex,nbasis> CoarseOperator;
  //  typedef SparseMatrixBase<CoarseVector> CoarseOperator;

  Aggregates     & _Aggregates;
  FineOperator   & _FineOperator;
  FineSmoother   & _PreSmoother;
  FineSmoother   & _PostSmoother;
  CoarseOperator & _CoarseOperator;
  CoarseSolver   & _CoarseSolve;

  int    level;  void Level(int lv) {level = lv; };

  MGPreconditioner(Aggregates &Agg,
		   FineOperator &Fine,
		   FineSmoother &PreSmoother,
		   FineSmoother &PostSmoother,
		   CoarseOperator &CoarseOperator_,
		   CoarseSolver &CoarseSolve_)
    : _Aggregates(Agg),
      _FineOperator(Fine),
      _PreSmoother(PreSmoother),
      _PostSmoother(PostSmoother),
      _CoarseOperator(CoarseOperator_),
      _CoarseSolve(CoarseSolve_),
      level(1)  {  }

  virtual void operator()(const FineField &in, FineField & out) 
  {
    auto CoarseGrid = _CoarseOperator.Grid();
    CoarseVector Csrc(CoarseGrid);
    CoarseVector Csol(CoarseGrid);
    FineField vec1(in.Grid());
    FineField vec2(in.Grid());

    std::cout<<GridLogMessage << "Calling PreSmoother " <<std::endl;

    //    std::cout<<GridLogMessage << "Calling PreSmoother input residual "<<norm2(in) <<std::endl;
    double t;
    // Fine Smoother
    t=-usecond();
    _PreSmoother(in,out);
    t+=usecond();

    std::cout<<GridLogMessage << "PreSmoother took "<< t/1000.0<< "ms" <<std::endl;

    // Update the residual
    _FineOperator.Op(out,vec1);  sub(vec1, in ,vec1);   
    //    std::cout<<GridLogMessage <<"Residual-1 now " <<norm2(vec1)<<std::endl;

    // Fine to Coarse 
    t=-usecond();
    _CoarseOperator.ProjectToSubspace(_Aggregates,Csrc,vec1);
    //    _Aggregates.ProjectToSubspace  (Csrc,vec1);
    t+=usecond();
    std::cout<<GridLogMessage << "Project to coarse took "<< t/1000.0<< "ms" <<std::endl;

    // Coarse correction
    t=-usecond();
    _CoarseSolve(Csrc,Csol);
    //Csol=Zero();
    t+=usecond();
    std::cout<<GridLogMessage << "Coarse solve took "<< t/1000.0<< "ms" <<std::endl;

    // Coarse to Fine
    t=-usecond();  
    _CoarseOperator.PromoteFromSubspace(_Aggregates,Csol,vec1);
    //    _Aggregates.PromoteFromSubspace(Csol,vec1); 
    add(out,out,vec1);
    t+=usecond();
    std::cout<<GridLogMessage << "Promote to this level took "<< t/1000.0<< "ms" <<std::endl;

    // Residual
    _FineOperator.Op(out,vec1);  sub(vec1 ,in , vec1);  
    //    std::cout<<GridLogMessage <<"Residual-2 now " <<norm2(vec1)<<std::endl;

    // Fine Smoother
    t=-usecond();
    _PostSmoother(vec1,vec2);
    t+=usecond();
    std::cout<<GridLogMessage << "PostSmoother took "<< t/1000.0<< "ms" <<std::endl;

    add( out,out,vec2);
    std::cout<<GridLogMessage << "Done " <<std::endl;
  }
};

template<class Fobj,class CComplex,int nbasis, class CoarseSolver>
class HDCRPreconditioner : public LinearFunction< Lattice<Fobj> > {
public:
  using LinearFunction<Lattice<Fobj> >::operator();
  
  typedef Aggregation<Fobj,CComplex,nbasis> Aggregates;
  typedef typename Aggregation<Fobj,CComplex,nbasis>::CoarseVector CoarseVector;
  typedef typename Aggregation<Fobj,CComplex,nbasis>::CoarseMatrix CoarseMatrix;
  typedef typename Aggregation<Fobj,CComplex,nbasis>::FineField    FineField;
  typedef LinearOperatorBase<FineField>                            FineOperator;
  typedef LinearFunction    <FineField>                            FineSmoother;
  //typedef CoarseCayleyFermion<Fobj,CComplex,nbasis> CoarseOperator;
  typedef SparseMatrixBase<CoarseVector> CoarseOperator;

  Aggregates     & _Aggregates;
  FineOperator   & _FineOperator;
  FineSmoother   & _PreSmoother;
  FineSmoother   & _PostSmoother;
  CoarseOperator & _CoarseOperator;
  CoarseSolver   & _CoarseSolve;

  int    level;  void Level(int lv) {level = lv; };

  HDCRPreconditioner(Aggregates &Agg,
		   FineOperator &Fine,
		   FineSmoother &PreSmoother,
		   FineSmoother &PostSmoother,
		   CoarseOperator &CoarseOperator_,
		   CoarseSolver &CoarseSolve_)
    : _Aggregates(Agg),
      _FineOperator(Fine),
      _PreSmoother(PreSmoother),
      _PostSmoother(PostSmoother),
      _CoarseOperator(CoarseOperator_),
      _CoarseSolve(CoarseSolve_),
      level(1)  {  }

  virtual void operator()(const FineField &in, FineField & out) 
  {
    auto CoarseGrid = _CoarseOperator.Grid();
    CoarseVector Csrc(CoarseGrid);
    CoarseVector g5Csrc(CoarseGrid);
    CoarseVector Csol(CoarseGrid);
    FineField vec1(in.Grid());
    FineField vec2(in.Grid());

    std::cout<<GridLogMessage<<"\t\t\t" << "Calling PreSmoother " <<std::endl;

    double t;
    // Fine Smoother
    t=-usecond();
    _PreSmoother(in,out);
    t+=usecond();

    std::cout<<GridLogMessage<<"\t\t\t" << "PreSmoother took "<< t/1000.0<< "ms" <<std::endl;
    // Update the residual
    _FineOperator.Op(out,vec1);  sub(vec1, in ,vec1);   

    // Fine to Coarse 
    // Based on a coarsening of G5R5 D
    // Solves Ddwf out = in
    // Coarse operator is g5R5 Ddwf solves g5R5 Ddwf out = in  
    t=-usecond();
    G5R5(vec2,vec1);
    _Aggregates.ProjectToSubspace  (Csrc,vec2);
    t+=usecond();
    std::cout<<GridLogMessage<<"\t\t\t" << "Project to coarse took "<< t/1000.0<< "ms" <<std::endl;

    // Coarse correction
    t=-usecond();
    _CoarseSolve(Csrc,Csol);
    t+=usecond();
    std::cout<<GridLogMessage<<"\t\t\t" << "Coarse solve took "<< t/1000.0<< "ms" <<std::endl;

    // Coarse to Fine
    t=-usecond();  
    _Aggregates.PromoteFromSubspace(Csol,vec1); 
    add(out,out,vec1);
    t+=usecond();
    std::cout<<GridLogMessage<<"\t\t\t" << "Promote to this level took "<< t/1000.0<< "ms" <<std::endl;

    // Residual
    _FineOperator.Op(out,vec1);  sub(vec1 ,in , vec1);  

    // Fine Smoother
    t=-usecond();
    _PostSmoother(vec1,vec2);
    t+=usecond();
    std::cout<<GridLogMessage<<"\t\t\t" << "PostSmoother took "<< t/1000.0<< "ms" <<std::endl;

    add( out,out,vec2);
  }
};

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  const int Ls=24;

  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);

  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  ///////////////////////////////////////////////////
  // Construct a coarsened grid; utility for this?
  ///////////////////////////////////////////////////
  std::vector<int> block ({2,2,2,2}); // 4,2,2,2 gets worse
  std::vector<int> blockc ({1,1,1,1});
  const int nbasis= 24;
  const int nbasisc= 40; // decrease, not improvement

  auto clatt = GridDefaultLatt();
  for(int d=0;d<clatt.size();d++){
    clatt[d] = clatt[d]/block[d];
  }
  auto cclatt = clatt;
  for(int d=0;d<clatt.size();d++){
    cclatt[d] = clatt[d]/blockc[d];
  }

  GridCartesian *Coarse4d =  SpaceTimeGrid::makeFourDimGrid(clatt, GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());;
  GridCartesian *Coarse5d =  SpaceTimeGrid::makeFiveDimGrid(Ls,Coarse4d);
  //  GridRedBlackCartesian * Coarse4dRB = SpaceTimeGrid::makeFourDimRedBlackGrid(Coarse4d);
  //  GridRedBlackCartesian * Coarse5dRB = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,Coarse4d);

  GridCartesian *CoarseCoarse4d =  SpaceTimeGrid::makeFourDimGrid(cclatt, GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());;
  GridCartesian *CoarseCoarse5d =  SpaceTimeGrid::makeFiveDimGrid(1,CoarseCoarse4d);
  //  GridRedBlackCartesian * CoarseCoarse4dRB = SpaceTimeGrid::makeFourDimRedBlackGrid(CoarseCoarse4d);
  GridRedBlackCartesian * CoarseCoarse5dRB = SpaceTimeGrid::makeFiveDimRedBlackGrid(1,CoarseCoarse4d);

  std::vector<int> seeds({1,2,3,4});
  GridParallelRNG          RNG5(FGrid);   RNG5.SeedFixedIntegers(seeds);
  GridParallelRNG          RNG4(UGrid);   RNG4.SeedFixedIntegers(seeds);
  GridParallelRNG          CRNG(Coarse4d);CRNG.SeedFixedIntegers(seeds);

  LatticeGaugeField Umu(UGrid); 
#if 0
  SU3::TepidConfiguration(RNG4,Umu);
  RealD M5=1.0;
#else
  std::string file("./ckpoint_lat.1000");
  FieldMetaData header;
  NerscIO::readConfiguration(Umu,header,file);
  RealD M5=1.8;
#endif

  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  std::cout<<GridLogMessage << "Building g5R5 hermitian DWF operator" <<std::endl;
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;

  RealD mass=0.00078;

  WilsonFermionR     Dw(Umu,*UGrid,*UrbGrid,-M5);
  DomainWallFermionR Ddwf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);
  DomainWallFermionR Dpv (Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,1.0,M5);

  typedef Aggregation<vSpinColourVector,vTComplex,nbasis>              Subspace;
  typedef CoarsenedMatrix<vSpinColourVector,vTComplex,nbasis>          CoarseOperator;
  typedef CoarseOperator::CoarseVector                                 CoarseVector;
  typedef CoarseOperator::siteVector siteVector;

  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  std::cout<<GridLogMessage << "Calling Aggregation class to build subspace" <<std::endl;
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  // How to find criticall mass?
  // WilsonFermionR     Dw_null(Umu,*UGrid,*UrbGrid,-0.75); //   600 iters
  // WilsonFermionR     Dw_null(Umu,*UGrid,*UrbGrid,-0.80); //   800 iters
  //  WilsonFermionR     Dw_null(Umu,*UGrid,*UrbGrid,-0.82); // 1023 iters
  //  WilsonFermionR     Dw_null(Umu,*UGrid,*UrbGrid,-0.85); // 1428 iters
  //  WilsonFermionR     Dw_null(Umu,*UGrid,*UrbGrid,-0.87); //  1900 iters
  //  WilsonFermionR     Dw_null(Umu,*UGrid,*UrbGrid,-0.90); // 3900   iters
  //  WilsonFermionR     Dw_null(Umu,*UGrid,*UrbGrid,-0.92); // 6200   iters
  //  WilsonFermionR     Dw_null(Umu,*UGrid,*UrbGrid,-0.94);  // 8882 iters
  WilsonFermionR     Dw_null(Umu,*UGrid,*UrbGrid,-0.95);  // 9170  iters
  //  WilsonFermionR     Dw_null(Umu,*UGrid,*UrbGrid,-0.96);  // 8882   iters
  //  WilsonFermionR     Dw_null(Umu,*UGrid,*UrbGrid,-0.97);  // 8406  iters
  //  WilsonFermionR     Dw_null(Umu,*UGrid,*UrbGrid,-0.99); // 6900   iters
  //  WilsonFermionR     Dw_null(Umu,*UGrid,*UrbGrid,-1.01); // 6397   iters
  //  WilsonFermionR     Dw_null(Umu,*UGrid,*UrbGrid,-1.00); // 5900   iters
  MdagMLinearOperator<WilsonFermionR,LatticeFermion> MdagM_Dw(Dw_null);

  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  std::cout<<GridLogMessage << "Testing Wilson criticality " <<std::endl;
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  /*
  ConjugateGradient<LatticeFermion>          WilsonCG(1.0e-10,40000);
  LatticeFermion w_src(UGrid); w_src=1.0;
  LatticeFermion w_res(UGrid);
  WilsonCG(MdagM_Dw,w_src,w_res);
  exit(0);
  */
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  std::cout<<GridLogMessage << " 4D subspace build                                " <<std::endl;
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  Subspace Aggregates4D(Coarse4d,UGrid,0);
  assert ( (nbasis & 0x1)==0);
  int nb=nbasis/2;
  Gamma g5(Gamma::Algebra::Gamma5);

  Aggregates4D.CreateSubspaceChebyshev(RNG4,MdagM_Dw,nb,60.0,3.0,300,150,150,0.0); // now at 26 iter
  //  Aggregates4D.CreateSubspaceChebyshev(RNG4,MdagM_Dw,nb,60.0,3.0,300,150,150,0.0); // now at 26 iter
  //  Aggregates4D.CreateSubspaceChebyshev(RNG4,MdagM_Dw,nb,60.0,3.0,500,150,150,0.0); // now at 26 iter
  //  Aggregates4D.CreateSubspaceChebyshev(RNG4,MdagM_Dw,nb,60.0,4.0,500,150,150,0.0); // now at 26 iter
  //  Aggregates4D.CreateSubspaceChebyshev(RNG4,MdagM_Dw,nb,60.0,4.0,500,150,150,0.0); //35
  //  Aggregates4D.CreateSubspaceChebyshev(RNG4,MdagM_Dw,nb,60.0,4.0,600,100,250,0.0); //39
  //  Aggregates4D.CreateSubspaceChebyshev(RNG4,MdagM_Dw,nb,60.0,4.0,600,250,100,0.0); //39
  //  Aggregates4D.CreateSubspaceChebyshev(RNG4,MdagM_Dw,nb,60.0,4.0,600,100,100,0.0); //36
  //  Aggregates4D.CreateSubspaceChebyshev(RNG4,MdagM_Dw,nb,60.0,4.0,600,250,100,0.0); //39
  //   Aggregates4D.CreateSubspaceChebyshev(RNG4,MdagM_Dw,nb,60.0,4.0,600,250,50,0.0);// 39
  //  Aggregates4D.CreateSubspaceChebyshev(RNG4,MdagM_Dw,nb,60.0,4.0,600,250,250,0.0);// 35
  //  Aggregates4D.CreateSubspaceChebyshev(RNG4,MdagM_Dw,nb,60.0,4.0,250,250,250,0.0);//  38 iter
  //  Aggregates4D.CreateSubspaceChebyshev(RNG4,MdagM_Dw,nb,50.0,4.0,250,250,100,0.0);//  40 iter
  //  Aggregates4D.CreateSubspaceChebyshev(RNG4,MdagM_Dw,nb,50.0,4.0,500,100,100,0.0);// 38 iter
  //  Aggregates4D.CreateSubspaceChebyshev(RNG4,MdagM_Dw,nb,60.0,4.0,500,100,100,0.0);// 36 iter
  //  Aggregates4D.CreateSubspaceChebyshev(RNG4,MdagM_Dw,nb,60.0,3.0,500,100,100,0.0);// 37 iter
  //  Aggregates4D.CreateSubspaceChebyshev(RNG4,MdagM_Dw,nb,60.0,1.0,1000,400,400,0.0);// 38 iter
  //  Aggregates4D.CreateSubspaceChebyshev(RNG4,MdagM_Dw,nb,60.0,0.5,1000,400,400,0.0);// 39 iter HDCR smooth 14
  //  Aggregates4D.CreateSubspaceChebyshev(RNG4,MdagM_Dw,nb,60.0,0.1,1000,400,400,0.0);// 41
  for(int n=0;n<nb;n++){
    Aggregates4D.subspace[nbasis-1-n]= Aggregates4D.subspace[n] - g5 * Aggregates4D.subspace[n];
    Aggregates4D.subspace[n]         = Aggregates4D.subspace[n] + g5 * Aggregates4D.subspace[n];
  }
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  std::cout<<GridLogMessage << " Coarsen the Dw operator                          " <<std::endl;
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  typedef CoarsenedMatrix<vSpinColourVector,vTComplex,nbasis>    Level1Op4;
  typedef CoarseCayleyFermion<vSpinColourVector,vTComplex,nbasis> Level1Op5;
  Level1Op4 c_Dw    (*Coarse4d,0);
  NonHermitianLinearOperator<WilsonFermionR,LatticeFermion>  LinOpDw(Dw);
  c_Dw.CoarsenOperator(UGrid,LinOpDw,Aggregates4D); // contains the M5 from Dw(-M5)
  //  c_Dw.Test(Aggregates4D,UGrid,LinOpDw);

  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  std::cout<<GridLogMessage << " Build coarse DWF operator                          " <<std::endl;
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  Level1Op5 c_Dwf  (*Coarse4d,*Coarse5d,c_Dw,M5, mass, Ls, 1.0,0.0);
  //  c_Dwf.Test(Aggregates4D,FGrid,Ddwf);

  MdagMLinearOperator<Level1Op5,CoarseVector> MdagM_cDwf(c_Dwf);

  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  std::cout<<GridLogMessage << "Build 5D coarse deflation space" << std::endl;
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  int nbc=nbasisc/2;
  typedef CoarsenedMatrix<siteVector,iScalar<vTComplex>,nbasisc>     Level2Op;
  typedef Aggregation<siteVector,iScalar<vTComplex>,nbasisc> CoarseSubspace;
  CoarseSubspace CoarseAggregates(CoarseCoarse5d,Coarse5d,0);

  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  std::cout<<GridLogMessage << "Build Chebyshev space in coarse operator "<< std::endl;
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  //  CoarseAggregates.CreateSubspaceChebyshev(CRNG,MdagM_cDwf,nbc,40.0,0.01,300,150,100,0.0);
  CoarseAggregates.CreateSubspaceChebyshev(CRNG,MdagM_cDwf,nbc,40.0,0.005,500,150,100,0.0);

  {
    std::cout<<GridLogMessage << "**************************************************"<< std::endl;
    std::cout<<GridLogMessage << "Applying G5R5 projection of coarse operator "<< std::endl;
    std::cout<<GridLogMessage << "**************************************************"<< std::endl;
    CoarseVector A(Coarse5d), B(Coarse5d);
    for(int n=0;n<nbc;n++){
      G5R5(B,CoarseAggregates.subspace[n]);
      A = CoarseAggregates.subspace[n];
      CoarseAggregates.subspace[n]    = A+B; // 1+G5R5 // eigen value of G5R5 is +1
      CoarseAggregates.subspace[n+nbc]= A-B; // 1-G5R5 // eigen value of G5R5 is -1
    }
  }

  Gamma5R5HermitianLinearOperator<Level1Op5,CoarseVector> L1Hdwf(c_Dwf);
  Level2Op cc_Dwf  (*CoarseCoarse5d,*CoarseCoarse5dRB,1); // say it is hermitian
  cc_Dwf.CoarsenOperator(Coarse5d,L1Hdwf,CoarseAggregates);
  //  cc_Dwf.Test(CoarseAggregates,Coarse5d,L1Hdwf);

  typedef Level2Op::CoarseVector CoarseCoarseVector;

  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  std::cout<<GridLogMessage << " Testing fine and coarse solvers " <<std::endl;
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;

  RealD tol=1.0e-8;
  int MaxIt = 10000;

  CoarseVector c_src(Coarse5d); c_src=1.0;
  CoarseVector c_res(Coarse5d);

  LatticeFermion f_src(FGrid); f_src=1.0;
  LatticeFermion f_res(FGrid);

  LatticeFermion f_src_e(FrbGrid); f_src_e=1.0;
  LatticeFermion f_res_e(FrbGrid);

  CoarseCoarseVector cc_src(CoarseCoarse5d); cc_src=1.0;

  ConjugateGradient<CoarseVector>            CoarseCG(tol,MaxIt);
  ConjugateGradient<LatticeFermion>          FineCG(tol,MaxIt);
  
  NonHermitianLinearOperator<DomainWallFermionR,LatticeFermion> FineM(Ddwf);
  MdagMLinearOperator<DomainWallFermionR,LatticeFermion>    FineMdagM(Ddwf);     //  M^\dag M

  NonHermitianLinearOperator<Level1Op5,CoarseVector> CoarseM(c_Dwf);
  MdagMLinearOperator<Level1Op5,CoarseVector> CoarseMdagM(c_Dwf);

  NonHermitianLinearOperator<Level2Op,CoarseCoarseVector> CoarseCoarseM(cc_Dwf);
  MdagMLinearOperator<Level2Op,CoarseCoarseVector> CoarseCoarseMdagM(cc_Dwf);


  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  std::cout<<GridLogMessage << " Fine    Hw PowerMethod           "<< std::endl;
  LatticeFermion w_src(UGrid); 
  w_src=1.0;
  PowerMethod<LatticeFermion>       PM;   PM(MdagM_Dw,w_src);
  std::cout<<GridLogMessage << " Coarse       PowerMethod           "<< std::endl;
  c_src=1.0;
  PowerMethod<CoarseVector>        cPM;  cPM(CoarseMdagM,c_src);

  cc_src=1.0;
  PowerMethod<CoarseCoarseVector>        ccPM;  ccPM(CoarseCoarseMdagM,cc_src);

  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  std::cout<<GridLogMessage << " Running CoarseCoarse grid Lanczos "<< std::endl;
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;

  // 37s, 26 iter
  int cNk=128;
  int cNm=256;
  int cNstop=128;
  RealD IRL_lo=0.01;
  RealD IRL_hi=16.0;
  int IRL_ord=201;

  /*
  //  int cNk=100; -- slower 27 iters
  int cNk=128; //-- 26 iters, but slower
  int cNm=192;
  int cNstop=128;
  RealD IRL_lo=0.005;
  RealD IRL_hi=10.0;
  int IRL_ord=101;
  */

  MdagMLinearOperator<Level2Op,CoarseCoarseVector> IRLHermOpL2(cc_Dwf);
  Chebyshev<CoarseCoarseVector> IRLChebyL2(IRL_lo,IRL_hi,IRL_ord);
  FunctionHermOp<CoarseCoarseVector> IRLOpChebyL2(IRLChebyL2,IRLHermOpL2);
  PlainHermOp<CoarseCoarseVector> IRLOpL2    (IRLHermOpL2);
  ImplicitlyRestartedLanczos<CoarseCoarseVector> IRLL2(IRLOpChebyL2,IRLOpL2,cNstop,cNk,cNm,1.0e-3,20);

  cNm=0;
  std::vector<RealD>          eval2(cNm);
  std::vector<CoarseCoarseVector>   evec2(cNm,CoarseCoarse5d);
  cc_src=1.0;
  //  int cNconv;
  //  IRLL2.calc(eval2,evec2,cc_src,cNconv);
  
  std::vector<RealD> tols ({0.005,0.001});
  std::vector<RealD> c_los  ({0.1,0.05});
  std::vector<RealD> c_his  ({22.0});
  std::vector<RealD> f_los  ({0.5,0.2});
  std::vector<RealD> f_his  ({60.0});
  std::vector<int> ws ({2,3});
  std::vector<int> c_ords ({32,24});
  std::vector<int> f_ords ({20,16});

  for(auto w : ws ) {
  for(auto tol : tols ) {
  for(auto f_ord : f_ords ) {
  for(auto c_ord : c_ords ) {
  for(auto c_lo : c_los ) {
  for(auto c_hi : c_his ) {
  for(auto f_lo : f_los ) {
  for(auto f_hi : f_his ) {
    //  ZeroGuesser<CoarseVector> CoarseZeroGuesser;
    //  ZeroGuesser<CoarseCoarseVector>       CoarseCoarseZeroGuesser;
  ConjugateGradient<CoarseCoarseVector>  CoarseCoarseCG(tol,10000);
  //  ZeroGuesser<CoarseCoarseVector> CoarseCoarseGuesser;
  SchurRedBlackDiagMooeeSolve<CoarseCoarseVector> CoarseCoarseRBCG(CoarseCoarseCG);
  SchurSolverWrapper<CoarseCoarseVector> CoarseCoarseSolver(cc_Dwf,CoarseCoarseRBCG);

  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  std::cout<<GridLogMessage << "Building 3 level hdcr                             "<< std::endl;
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  //  NormalEquations<CoarseCoarseVector>   CoarseCoarseCGNE(cc_Dwf,CoarseCoarseCG,CoarseCoarseZeroGuesser);
  {
typedef HDCRPreconditioner<siteVector,iScalar<vTComplex>,nbasisc,LinearFunction<CoarseCoarseVector> > CoarseMG;
  typedef MGPreconditioner<vSpinColourVector,  vTComplex,nbasis, LinearFunction<CoarseVector> >     ThreeLevelMG;

  // MultiGrid preconditioner acting on the coarse space <-> coarsecoarse space
  //  ChebyshevSmoother<CoarseVector,  Level1Op5 >       CoarseSmoother1(0.5,22.0,c_ord,CoarseM,c_Dwf); // 37s, 26 iter
  //  ChebyshevSmoother<CoarseVector,  Level1Op5 >       CoarseSmoother2(0.5,22.0,c_ord,CoarseM,c_Dwf);
  ChebyshevSmoother<CoarseVector,  Level1Op5 >       CoarseSmoother(c_lo,c_hi,c_ord,CoarseM,c_Dwf); // 37s, 26 iter

  //  ChebyshevSmoother<CoarseVector,  Level1Op5 >       CoarseSmoother1(0.5,22.0,7,CoarseM,c_Dwf); // 38s, 26 iter
  //  ChebyshevSmoother<CoarseVector,  Level1Op5 >       CoarseSmoother2(0.5,22.0,7,CoarseM,c_Dwf);
  //  ChebyshevSmoother<CoarseVector,  Level1Op5 >       CoarseSmoother1(0.4,22.0,7,CoarseM,c_Dwf); // 41s, 27 iter
  //  ChebyshevSmoother<CoarseVector,  Level1Op5 >       CoarseSmoother2(0.4,22.0,7,CoarseM,c_Dwf);
  //  ChebyshevSmoother<CoarseVector,  Level1Op5 >       CoarseSmoother1(0.6,22.0,6,CoarseM,c_Dwf); // 26 iter
  //  ChebyshevSmoother<CoarseVector,  Level1Op5 >       CoarseSmoother2(0.6,22.0,6,CoarseM,c_Dwf);
  //  ChebyshevSmoother<CoarseVector,  Level1Op5 >       CoarseSmoother1(0.5,22.0,5,CoarseM,c_Dwf); // 33 iter, 55s
  //  ChebyshevSmoother<CoarseVector,  Level1Op5 >       CoarseSmoother2(0.5,22.0,5,CoarseM,c_Dwf);


  CoarseMG Level2Precon (CoarseAggregates,
			 CoarseM,
			 CoarseSmoother,
			 CoarseSmoother,
			 cc_Dwf,
			 CoarseCoarseSolver);
  Level2Precon.Level(2);

  //PrecGeneralisedConjugateResidualNonHermitian<CoarseVector>  L2PGCR(0.5, 100, CoarseM,Level2Precon,16,16); // 26 iter, 37s
  //  PrecGeneralisedConjugateResidualNonHermitian<CoarseVector>  L2PGCR(0.0, 1, CoarseM,Level2Precon,2,2);  // 296 s, 50 iter
  //  PrecGeneralisedConjugateResidualNonHermitian<CoarseVector>  L2PGCR(0.0, 1, CoarseM,Level2Precon,2,2);  // 250 s, 37 iter
  PrecGeneralisedConjugateResidualNonHermitian<CoarseVector>  L2PGCR(0.0, 1, CoarseM,Level2Precon,2,2); 

  //PrecGeneralisedConjugateResidualNonHermitian<CoarseVector>  L2PGCR(1.0, 100, CoarseM,Level2Precon,16,16); // 35 iter, 45s
  //PrecGeneralisedConjugateResidualNonHermitian<CoarseVector>  L2PGCR(0.6, 100, CoarseM,Level2Precon,16,16); // 26,38 (diifferene is measurement noise)
  //PrecGeneralisedConjugateResidualNonHermitian<CoarseVector>  L2PGCR(0.2, 100, CoarseM,Level2Precon,16,16); // 26 iter, 47s
  L2PGCR.Level(2);

  // Wrap the 2nd level solver in a MultiGrid preconditioner acting on the fine space

  //  ChebyshevSmoother<LatticeFermion,DomainWallFermionR> FineSmoother1(0.5,60.0,14,FineM,Ddwf); // 26 iter, 39s
  //  ChebyshevSmoother<LatticeFermion,DomainWallFermionR> FineSmoother2(0.5,60.0,14,FineM,Ddwf);

  //  ChebyshevSmoother<LatticeFermion,DomainWallFermionR> FineSmoother1(0.5,60.0,12,FineM,Ddwf); // 25 iter, 38s
  //  ChebyshevSmoother<LatticeFermion,DomainWallFermionR> FineSmoother2(0.5,60.0,16,FineM,Ddwf);

  //  ChebyshevSmoother<LatticeFermion,DomainWallFermionR> FineSmoother1(0.5,60.0,12,FineM,Ddwf); // 23 iter, 39s
  //  ChebyshevSmoother<LatticeFermion,DomainWallFermionR> FineSmoother2(0.5,60.0,20,FineM,Ddwf);

  //  ChebyshevSmoother<LatticeFermion,DomainWallFermionR> FineSmoother1(0.5,60.0,10,FineM,Ddwf);24 iter, 44s
  //  ChebyshevSmoother<LatticeFermion,DomainWallFermionR> FineSmoother2(0.5,60.0,24,FineM,Ddwf);

  //  ChebyshevSmoother<LatticeFermion,DomainWallFermionR> FineSmoother1(0.5,60.0,12,FineM,Ddwf); // odd convergence tail at 10^-9 ish
  //  ChebyshevSmoother<LatticeFermion,DomainWallFermionR> FineSmoother2(0.1,60.0,24,FineM,Ddwf); // 33 iter, waas O(10-9 by 26)

  //  ChebyshevSmoother<LatticeFermion,DomainWallFermionR> FineSmoother1(0.5,60.0,12,FineM,Ddwf); // 25 iter, 39s
  //  ChebyshevSmoother<LatticeFermion,DomainWallFermionR> FineSmoother2(0.5,60.0,18,FineM,Ddwf); //

  ChebyshevSmoother<LatticeFermion,DomainWallFermionR> FineSmoother(f_lo,f_hi,f_ord,FineM,Ddwf); 

  //  ChebyshevSmoother<LatticeFermion,DomainWallFermionR> FineSmoother1(0.5,60.0,11,FineM,Ddwf); // 33 iter, 49s
  //  ChebyshevSmoother<LatticeFermion,DomainWallFermionR> FineSmoother2(0.5,60.0,11,FineM,Ddwf);
  //  ChebyshevSmoother<LatticeFermion,DomainWallFermionR> FineSmoother1(0.5,60.0,12,FineM,Ddwf); // 26 iter, 37s
  //  ChebyshevSmoother<LatticeFermion,DomainWallFermionR> FineSmoother2(0.5,60.0,12,FineM,Ddwf);
  //  ChebyshevSmoother<LatticeFermion,DomainWallFermionR> FineSmoother1(0.4,60.0,12,FineM,Ddwf); //  iter 26 no change in final residual
  //  ChebyshevSmoother<LatticeFermion,DomainWallFermionR> FineSmoother2(0.4,60.0,12,FineM,Ddwf);
  //  ChebyshevSmoother<LatticeFermion,DomainWallFermionR> FineSmoother1(0.3,60.0,12,FineM,Ddwf); // 27 iter 39s.
  //  ChebyshevSmoother<LatticeFermion,DomainWallFermionR> FineSmoother2(0.3,60.0,12,FineM,Ddwf);
  //  ChebyshevSmoother<LatticeFermion,DomainWallFermionR> FineSmoother1(0.3,60.0,13,FineM,Ddwf); // 26 iter, but slower
  //  ChebyshevSmoother<LatticeFermion,DomainWallFermionR> FineSmoother2(0.3,60.0,13,FineM,Ddwf);
  //  ChebyshevSmoother<LatticeFermion,DomainWallFermionR> FineSmoother1(1.0,60.0,12,FineM,Ddwf); // 34 iter, slower
  //  ChebyshevSmoother<LatticeFermion,DomainWallFermionR> FineSmoother2(1.0,60.0,12,FineM,Ddwf);

  ThreeLevelMG ThreeLevelPrecon(Aggregates4D,
				FineM,
				FineSmoother,
				FineSmoother,
				c_Dwf,
				L2PGCR);
  ThreeLevelPrecon.Level(1);

  PrecGeneralisedConjugateResidualNonHermitian<LatticeFermion> L1PGCR(1.0e-8,1000,FineM,ThreeLevelPrecon,16,16);
  L1PGCR.Level(1);

  f_res=Zero();
  L1PGCR(f_src,f_res);
  }
  }}}}  
  }}}
  }
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  std::cout<<GridLogMessage << "Done "<< std::endl;
  std::cout<<GridLogMessage << "**************************************************"<< std::endl;
  Grid_finalize();
  
}
