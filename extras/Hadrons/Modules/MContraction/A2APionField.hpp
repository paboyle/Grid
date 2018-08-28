#ifndef Hadrons_MContraction_A2APionField_hpp_
#define Hadrons_MContraction_A2APionField_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>
#include <Grid/Hadrons/AllToAllVectors.hpp>

#include <unsupported/Eigen/CXX11/Tensor>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         A2APionField                                       *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

typedef std::pair<Gamma::Algebra, Gamma::Algebra> GammaPair;


class A2APionFieldPar : Serializable
{
  public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2APionFieldPar,
				    int, cacheBlock,
				    int, schurBlock,
				    int, Nmom,
                                    std::string, A2A_i,
                                    std::string, A2A_j,
                                    std::string, output);
};

template <typename FImpl>
class TA2APionField : public Module<A2APionFieldPar>
{
  public:
    FERM_TYPE_ALIASES(FImpl, );
    SOLVER_TYPE_ALIASES(FImpl, );

    typedef A2AModesSchurDiagTwo<typename FImpl::FermionField, FMat, Solver> A2ABase;

  int d_unroll ;

  public:
    // constructor
    TA2APionField(const std::string name);
    // destructor
    virtual ~TA2APionField(void){};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
  
    virtual void DeltaFeq2(int dt_min,int dt_max,
			   Eigen::Tensor<ComplexD,2> &dF2_fig8,
			   Eigen::Tensor<ComplexD,2> &dF2_trtr,
			   Eigen::Tensor<ComplexD,3> &WW_sd, 
			   const LatticeFermion *vs,
			   const LatticeFermion *vd,
			   int orthogdim);

    virtual void DeltaFeq2_alt(int dt_min,int dt_max,
			       Eigen::Tensor<ComplexD,2> &dF2_fig8,
			       Eigen::Tensor<ComplexD,2> &dF2_trtr,
			       Eigen::Tensor<ComplexD,1> &den0,
			       Eigen::Tensor<ComplexD,1> &den1,
			       Eigen::Tensor<ComplexD,3> &WW_sd, 
			       const LatticeFermion *vs,
			       const LatticeFermion *vd,
			       int orthogdim);

    ///////////////////////////////////////
    // Arithmetic help. Move to Grid??
    ///////////////////////////////////////
    virtual void PionFieldXX(Eigen::Tensor<ComplexD,3> &mat, 
			     const LatticeFermion *wi,
			     const LatticeFermion *vj,
			     int orthogdim,
			     int g5);      
    ///////////////////////////
    // Simple wrappers
    ///////////////////////////
    virtual void PionFieldWVmom(Eigen::Tensor<ComplexD,4> &mat, 
				const LatticeFermion *wi,
				const LatticeFermion *vj,
				const std::vector<LatticeComplex > &mom,
				int orthogdim);      

    virtual void PionFieldWV(Eigen::Tensor<ComplexD,3> &mat, 
			     const LatticeFermion *wi,
			     const LatticeFermion *vj,
			     int orthogdim);      

    virtual void PionFieldVV(Eigen::Tensor<ComplexD,3> &mat, 
			     const LatticeFermion *vi,
			     const LatticeFermion *vj,
			     int orthogdim);      

    virtual void PionFieldWW(Eigen::Tensor<ComplexD,3> &mat, 
			     const LatticeFermion *wi,
			     const LatticeFermion *wj,
			     int orthogdim);      

};

MODULE_REGISTER(A2APionField, ARG(TA2APionField<FIMPL>), MContraction);
MODULE_REGISTER(ZA2APionField, ARG(TA2APionField<ZFIMPL>), MContraction);

/******************************************************************************
*                  TA2APionField implementation                             *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TA2APionField<FImpl>::TA2APionField(const std::string name)
    : Module<A2APionFieldPar>(name)
{
}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TA2APionField<FImpl>::getInput(void)
{
  std::vector<std::string> in;
  in.push_back(par().A2A_i + "_class");
  in.push_back(par().A2A_i + "_w_high_4d");
  in.push_back(par().A2A_i + "_v_high_4d");
  in.push_back(par().A2A_j + "_class");
  in.push_back(par().A2A_j + "_w_high_4d");
  in.push_back(par().A2A_j + "_v_high_4d");
  
  return in;
}

template <typename FImpl>
std::vector<std::string> TA2APionField<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {};

    return out;
}


// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2APionField<FImpl>::setup(void)
{
  d_unroll = 32; // empirical default. Can be overridden if desired

  // Four D fields
  envTmp(std::vector<FermionField>, "wi", 1, par().schurBlock, FermionField(env().getGrid(1)));
  envTmp(std::vector<FermionField>, "vi", 1, par().schurBlock, FermionField(env().getGrid(1)));
  envTmp(std::vector<FermionField>, "wj", 1, par().schurBlock, FermionField(env().getGrid(1)));
  envTmp(std::vector<FermionField>, "vj", 1, par().schurBlock, FermionField(env().getGrid(1)));

  // 5D tmp
  int Ls_i = env().getObjectLs(par().A2A_i + "_class");
  envTmpLat(FermionField, "tmp_5d", Ls_i);
  
  int Ls_j= env().getObjectLs(par().A2A_j + "_class");
  assert ( Ls_i == Ls_j ); 
}


//////////////////////////////////////////////////////////////////////////////////
// Cache blocked arithmetic routine
// Could move to Grid ???
//////////////////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2APionField<FImpl>::PionFieldWVmom(Eigen::Tensor<ComplexD,4> &mat, 
					  const LatticeFermion *wi,
					  const LatticeFermion *vj,
					  const std::vector<LatticeComplex > &mom,
					  int orthogdim) 
{
  double t0,t1,t2,t3;  t0=t1=t2=t3= 0.0;

  typedef typename FImpl::SiteSpinor vobj;

  typedef typename vobj::scalar_object sobj;
  typedef typename vobj::scalar_type scalar_type;
  typedef typename vobj::vector_type vector_type;

  typedef iSpinMatrix<vector_type> SpinMatrix_v;
  typedef iSpinMatrix<scalar_type> SpinMatrix_s;
  
  int Lblock = mat.dimension(2); 
  int Rblock = mat.dimension(3);

  GridBase *grid = wi[0]._grid;
  
  const int    nd = grid->_ndimension;
  const int Nsimd = grid->Nsimd();

  int Nt     = grid->GlobalDimensions()[orthogdim];
  int Nmom   = mom.size();

  int fd=grid->_fdimensions[orthogdim];
  int ld=grid->_ldimensions[orthogdim];
  int rd=grid->_rdimensions[orthogdim];

  // will locally sum vectors first
  // sum across these down to scalars
  // splitting the SIMD
  int MFrvol = rd*Lblock*Rblock*Nmom;
  int MFlvol = ld*Lblock*Rblock*Nmom;

  Vector<vector_type > lvSum(MFrvol);
  parallel_for (int r = 0; r < MFrvol; r++){
    lvSum[r] = zero;
  }

  Vector<scalar_type > lsSum(MFlvol);             
  parallel_for (int r = 0; r < MFlvol; r++){
    lsSum[r]=scalar_type(0.0);
  }

  int e1=    grid->_slice_nblock[orthogdim];
  int e2=    grid->_slice_block [orthogdim];
  int stride=grid->_slice_stride[orthogdim];

  t0-=usecond();
  parallel_for(int r=0;r<rd;r++){

    int so=r*grid->_ostride[orthogdim]; // base offset for start of plane 

    for(int n=0;n<e1;n++){
      for(int b=0;b<e2;b++){

	int ss= so+n*stride+b;

	for(int i=0;i<Lblock;i++){

	  auto w = conjugate(wi[i]._odata[ss]);

	  for(int j=0;j<Rblock;j++){

	    auto v = vj[j]._odata[ss];

	    auto vv = w()(0)(0) * v()(0)(0)// Gamma5 Dirac basis explicitly written out
	      +       w()(0)(1) * v()(0)(1)
	      +       w()(0)(2) * v()(0)(2)
	      +       w()(1)(0) * v()(1)(0)
	      +       w()(1)(1) * v()(1)(1)
	      +       w()(1)(2) * v()(1)(2)
	      -       w()(2)(0) * v()(2)(0)
	      -       w()(2)(1) * v()(2)(1)
	      -       w()(2)(2) * v()(2)(2)
	      -       w()(3)(0) * v()(3)(0)
	      -       w()(3)(1) * v()(3)(1)
	      -       w()(3)(2) * v()(3)(2);

	    
	    // After getting the sitewise product do the mom phase loop
	    int base = Nmom*i+Nmom*Lblock*j+Nmom*Lblock*Rblock*r;
	    for ( int m=0;m<Nmom;m++){
	      int idx = m+base;
	      auto phase = mom[m]._odata[ss];
	      mac(&lvSum[idx],&vv,&phase()()());
	    }
	  }
	}
      }
    }
  }
  t0+=usecond();


  // Sum across simd lanes in the plane, breaking out orthog dir.
  t1-=usecond();
  parallel_for(int rt=0;rt<rd;rt++){

    std::vector<int> icoor(nd);
    iScalar<vector_type> temp; 
    std::vector<iScalar<scalar_type> > extracted(Nsimd);               

      //    std::vector<scalar_type> extracted(Nsimd);               

    for(int i=0;i<Lblock;i++){
    for(int j=0;j<Rblock;j++){
    for(int m=0;m<Nmom;m++){

      int ij_rdx = m+Nmom*i+Nmom*Lblock*j+Nmom*Lblock*Rblock*rt;

      temp._internal = lvSum[ij_rdx];
      extract(temp,extracted);

      for(int idx=0;idx<Nsimd;idx++){

	grid->iCoorFromIindex(icoor,idx);

	int ldx    = rt+icoor[orthogdim]*rd;

	int ij_ldx = m+Nmom*i+Nmom*Lblock*j+Nmom*Lblock*Rblock*ldx;

	lsSum[ij_ldx]=lsSum[ij_ldx]+extracted[idx]._internal;

      }
    }}}
  }
  t1+=usecond();

  assert(mat.dimension(0) == Nmom);
  assert(mat.dimension(1) == Nt);
  t2-=usecond();
  // ld loop and local only??
  int pd = grid->_processors[orthogdim];
  int pc = grid->_processor_coor[orthogdim];
  parallel_for_nest2(int lt=0;lt<ld;lt++)
  {
    for(int pt=0;pt<pd;pt++){
      int t = lt + pt*ld;
      if (pt == pc){
	for(int i=0;i<Lblock;i++){
	  for(int j=0;j<Rblock;j++){
	    for(int m=0;m<Nmom;m++){
	      int ij_dx = m+Nmom*i + Nmom*Lblock * j + Nmom*Lblock * Rblock * lt;
	      mat(m,t,i,j) = lsSum[ij_dx];
	    }
	  }
	}
      } else { 
	const scalar_type zz(0.0);
	for(int i=0;i<Lblock;i++){
	  for(int j=0;j<Rblock;j++){
	    for(int m=0;m<Nmom;m++){
	      mat(m,t,i,j) =zz;
	    }
	  }
	}
      }
    }
  }
  t2+=usecond();

  t3-=usecond();
  grid->GlobalSumVector(&mat(0,0,0,0),Nmom*Nt*Lblock*Rblock);
  t3+=usecond();
}

///////////////////////////////////////////////////////////////////
//Meson 
// Interested in
//
//      sum_x,y Trace[ G S(x,tx,y,ty) G S(y,ty,x,tx) ]
//
// Conventional meson field:
//                 
//    = sum_x,y Trace[ sum_j G |v_j(y,ty)> <w_j(x,tx)|  G sum_i |v_i(x,tx) ><w_i(y,ty)| ]
//    = sum_ij sum_x,y < w_j(x,tx)| G |v_i(x,tx) > <w_i(y,ty) (x)|G| v_j(y,ty) >
//    = sum_ij PI_ji(tx) PI_ij(ty)
//
// G5-Hermiticity
//
//      sum_x,y Trace[ G S(x,tx,y,ty) G S(y,ty,x,tx) ]
//    = sum_x,y Trace[ G S(x,tx,y,ty) G g5 S^dag(x,tx,y,ty) g5 ]
//    = sum_x,y Trace[ g5 G sum_j |v_j(y,ty)> <w_j(x,tx)|  G g5 sum_i   (|v_j(y,ty)> <w_i(x,tx)|)^dag ]      --  (*)
//
// NB:  Dag applies to internal indices spin,colour,complex
//
//    = sum_ij sum_x,y Trace[ g5 G |v_j(y,ty)> <w_j(x,tx)|  G g5  |w_i(x,tx)> <v_i(y,ty)| ]
//    = sum_ij sum_x,y <v_i(y,ty)|g5 G |v_j(y,ty)> <w_j(x,tx)|  G g5 |w_i(x,tx)> 
//    = sum_ij  PionVV(ty) PionWW(tx)
//
// (*) is only correct estimator if w_i and w_j come from distinct noise sets to preserve the kronecker
//     expectation value. Otherwise biased.
////////////////////////////////////////////////////////////////////


template <typename FImpl>
void TA2APionField<FImpl>::PionFieldXX(Eigen::Tensor<ComplexD,3> &mat, 
				       const LatticeFermion *wi,
				       const LatticeFermion *vj,
				       int orthogdim,
				       int g5) 
{
  typedef typename FImpl::SiteSpinor vobj;

  typedef typename vobj::scalar_object sobj;
  typedef typename vobj::scalar_type scalar_type;
  typedef typename vobj::vector_type vector_type;

  typedef iSpinMatrix<vector_type> SpinMatrix_v;
  typedef iSpinMatrix<scalar_type> SpinMatrix_s;
  
  int Lblock = mat.dimension(1); 
  int Rblock = mat.dimension(2);

  GridBase *grid = wi[0]._grid;
  
  const int    nd = grid->_ndimension;
  const int Nsimd = grid->Nsimd();

  int Nt     = grid->GlobalDimensions()[orthogdim];

  int fd=grid->_fdimensions[orthogdim];
  int ld=grid->_ldimensions[orthogdim];
  int rd=grid->_rdimensions[orthogdim];

  // will locally sum vectors first
  // sum across these down to scalars
  // splitting the SIMD
  int MFrvol = rd*Lblock*Rblock;
  int MFlvol = ld*Lblock*Rblock;

  Vector<vector_type > lvSum(MFrvol);
  parallel_for (int r = 0; r < MFrvol; r++){
    lvSum[r] = zero;
  }

  Vector<scalar_type > lsSum(MFlvol);             
  parallel_for (int r = 0; r < MFlvol; r++){
    lsSum[r]=scalar_type(0.0);
  }

  int e1=    grid->_slice_nblock[orthogdim];
  int e2=    grid->_slice_block [orthogdim];
  int stride=grid->_slice_stride[orthogdim];

  parallel_for(int r=0;r<rd;r++){

    int so=r*grid->_ostride[orthogdim]; // base offset for start of plane 

    for(int n=0;n<e1;n++){
      for(int b=0;b<e2;b++){

	int ss= so+n*stride+b;

	for(int i=0;i<Lblock;i++){

	  auto w = conjugate(wi[i]._odata[ss]);

	  for(int j=0;j<Rblock;j++){

	    auto v = vj[j]._odata[ss];
	    auto vv = v()(0)(0);
	    if (g5) {
	         vv = w()(0)(0) * v()(0)(0)// Gamma5 Dirac basis explicitly written out
	      +       w()(0)(1) * v()(0)(1)
	      +       w()(0)(2) * v()(0)(2)
	      +       w()(1)(0) * v()(1)(0)
	      +       w()(1)(1) * v()(1)(1)
  	      +       w()(1)(2) * v()(1)(2)
	      -       w()(2)(0) * v()(2)(0)
	      -       w()(2)(1) * v()(2)(1)
	      -       w()(2)(2) * v()(2)(2)
	      -       w()(3)(0) * v()(3)(0)
	      -       w()(3)(1) * v()(3)(1)
	      -       w()(3)(2) * v()(3)(2);
	    } else {
	         vv = w()(0)(0) * v()(0)(0)// Gamma5 Dirac basis explicitly written out
	      +       w()(0)(1) * v()(0)(1)
	      +       w()(0)(2) * v()(0)(2)
	      +       w()(1)(0) * v()(1)(0)
	      +       w()(1)(1) * v()(1)(1)
	      +       w()(1)(2) * v()(1)(2)
	      +       w()(2)(0) * v()(2)(0)
	      +       w()(2)(1) * v()(2)(1)
	      +       w()(2)(2) * v()(2)(2)
	      +       w()(3)(0) * v()(3)(0)
	      +       w()(3)(1) * v()(3)(1)
	      +       w()(3)(2) * v()(3)(2);
	    }
	    
	    int idx = i+Lblock*j+Lblock*Rblock*r;
	    lvSum[idx] = lvSum[idx]+vv;
	  }
	}
      }
    }
  }

  // Sum across simd lanes in the plane, breaking out orthog dir.
  parallel_for(int rt=0;rt<rd;rt++){

    std::vector<int> icoor(nd);
    iScalar<vector_type> temp; 
    std::vector<iScalar<scalar_type> > extracted(Nsimd);               

    for(int i=0;i<Lblock;i++){
    for(int j=0;j<Rblock;j++){

      int ij_rdx = i+Lblock*j+Lblock*Rblock*rt;

      temp._internal =lvSum[ij_rdx];
      extract(temp,extracted);

      for(int idx=0;idx<Nsimd;idx++){

	grid->iCoorFromIindex(icoor,idx);

	int ldx    = rt+icoor[orthogdim]*rd;

	int ij_ldx =i+Lblock*j+Lblock*Rblock*ldx;

	lsSum[ij_ldx]=lsSum[ij_ldx]+extracted[idx]._internal;

      }
    }}
  }

  assert(mat.dimension(0) == Nt);
  // ld loop and local only??
  int pd = grid->_processors[orthogdim];
  int pc = grid->_processor_coor[orthogdim];
  parallel_for_nest2(int lt=0;lt<ld;lt++)
  {
    for(int pt=0;pt<pd;pt++){
      int t = lt + pt*ld;
      if (pt == pc){
	for(int i=0;i<Lblock;i++){
	  for(int j=0;j<Rblock;j++){
	    int ij_dx = i + Lblock * j + Lblock * Rblock * lt;
	    mat(t,i,j) = lsSum[ij_dx];
	  }
	}
      } else { 
	const scalar_type zz(0.0);
	for(int i=0;i<Lblock;i++){
	  for(int j=0;j<Rblock;j++){
	    mat(t,i,j) =zz;
	  }
	}
      }
    }
  }

  grid->GlobalSumVector(&mat(0,0,0),Nt*Lblock*Rblock);
}

template <typename FImpl>
void TA2APionField<FImpl>::PionFieldWV(Eigen::Tensor<ComplexD,3> &mat, 
				       const LatticeFermion *wi,
				       const LatticeFermion *vj,
				       int orthogdim) 
{
  const int g5=1;
  PionFieldXX(mat,wi,vj,orthogdim,g5);
}
template <typename FImpl>
void TA2APionField<FImpl>::PionFieldWW(Eigen::Tensor<ComplexD,3> &mat, 
				       const LatticeFermion *wi,
				       const LatticeFermion *wj,
				       int orthogdim) 
{
  const int nog5=0;
  PionFieldXX(mat,wi,wj,orthogdim,nog5);
}
template <typename FImpl>
void TA2APionField<FImpl>::PionFieldVV(Eigen::Tensor<ComplexD,3> &mat, 
				       const LatticeFermion *vi,
				       const LatticeFermion *vj,
				       int orthogdim) 
{
  const int nog5=0;
  PionFieldXX(mat,vi,vj,orthogdim,nog5);
}


/////////////////////////////////////////////////////////////////////////
// New dirac trace code needed for efficiency (I think).
// TODO: Ask Antonin to auto gen from Mathemetica
/////////////////////////////////////////////////////////////////////////
template<class vtype>
inline iScalar<vtype> traceGammaZ(const iMatrix<vtype, Ns> &rhs)
{
  iScalar<vtype> ret;
  ret() = timesI(rhs(2,0)) + timesMinusI(rhs(3,1)) + timesMinusI(rhs(0,2)) + timesI(rhs(1,3));
  return ret;
};
template<class vtype>
inline iScalar<vtype> traceGammaZGamma5(const iMatrix<vtype, Ns> &rhs)
{
  iScalar<vtype> ret;
  ret() = timesMinusI(rhs(2,0)) + timesI(rhs(3,1)) + timesMinusI(rhs(0,2)) + timesI(rhs(1,3));
  return ret;
};
template<class vtype>
inline iScalar<vtype> traceGammaY(const iMatrix<vtype, Ns> &rhs)
{
  iScalar<vtype> ret;
  ret() = -rhs(3,0) + rhs(2,1) + rhs(1,2) -rhs(0,3);
  return ret;
};
template<class vtype>
inline iScalar<vtype> traceGammaYGamma5(const iMatrix<vtype, Ns> &rhs)
{
  iScalar<vtype> ret;
  ret() = rhs(3,0) - rhs(2,1) + rhs(1,2) -rhs(0,3);
  return ret;
};
template<class vtype>
inline iScalar<vtype> traceGamma5(const iMatrix<vtype, Ns> &rhs)
{
  iScalar<vtype> ret;
  ret() = rhs(0, 0) + rhs(1, 1) - rhs(2, 2) - rhs(3,3);
  return ret;
};
template<class vtype>
inline iScalar<vtype> traceGammaT(const iMatrix<vtype, Ns> &rhs)
{
  iScalar<vtype> ret;
  ret() = rhs(2, 0) + rhs(3, 1) + rhs(0, 2) + rhs(1,3);
  return ret;
};
template<class vtype>
inline iScalar<vtype> traceGammaTGamma5(const iMatrix<vtype, Ns> &rhs)
{
  iScalar<vtype> ret;
  ret() = -rhs(2, 0) - rhs(3, 1) + rhs(0, 2) + rhs(1,3);
  return ret;
};
template<class vtype>
inline iScalar<vtype> traceGammaX(const iMatrix<vtype, Ns> &rhs)
{
  iScalar<vtype> ret;
  ret() = timesMinusI(rhs(0, 3)) +timesMinusI(rhs(1, 2)) 
        + timesI(rhs(2, 1)) +timesI(rhs(3,0));
  return ret;
};
template<class vtype>
inline iScalar<vtype> traceGammaXGamma5(const iMatrix<vtype, Ns> &rhs)
{
  iScalar<vtype> ret;
  ret() = timesMinusI(rhs(0, 3)) +timesMinusI(rhs(1, 2)) 
        + timesMinusI(rhs(2, 1)) +timesMinusI(rhs(3,0));
  return ret;
};


////////////////////////////////////////////////////////////////////////
// DeltaF=2 contraction ; use exernal WW field for Kaon, anti Kaon sink
////////////////////////////////////////////////////////////////////////
//
// WW -- i vectors have adjoint, and j vectors not. 
//    -- Think of "i" as the strange quark, forward prop from 0
//    -- Think of "j" as the anti-down quark.
//
// WW_sd are w^dag_s w_d
//
// Hence VV vectors correspondingly are  v^dag_d,  v_s    from t=0
//                                  and  v^dag_d,  v_s    from t=dT
//
// There is an implicit g5 associated with each v^dag_d from use of g5 Hermiticity.
// The other gamma_5 lies in the WW external meson operator.
//
// From UKhadron wallbag.cc:
//
//   LatticePropagator anti_d0 =  adj( Gamma(G5) * Qd0 * Gamma(G5));
//   LatticePropagator anti_d1 =  adj( Gamma(G5) * Qd1 * Gamma(G5));
//
//   PR1 = Qs0 * Gamma(G5) * anti_d0;
//   PR2 = Qs1 * Gamma(G5) * anti_d1;
//
//   TR1 = trace( PR1 * G1 );
//   TR2 = trace( PR2 * G2 );
//   Wick1 = TR1 * TR2;
//
//   Wick2 = trace( PR1* G1 * PR2 * G2 );
//   // was      Wick2 = trace( Qs0 * Gamma(G5) * anti_d0 * G1 * Qs1 * Gamma(G5) * anti_d1 * G2 );
//
// TR TR(tx) = Wick1 = sum_x WW[t0]_sd < v^_d |g5 G| v_s>   WW[t1]_s'd' < v^_d' |g5 G| v_s'> |_{x,tx)
//           = sum_x [ Trace(WW[t0] VgV(t,x) )  x Trace( WW_[t1] VgV(t,x) ) ]
//           
//
// Calc all Nt Trace(WW VV) products at once, take Nt^2 products of these.
//
// Fig8(tx)  = Wick2 = sum_x WW[t0]_sd WW[t1]_s'd'  < v^_d |g5 G| v_s'> < v^_d' |g5 G| v_s> |_{x,tx}
//
//                   = sum_x Trace( WW[t0] VV[t,x] WW[t1] VV[t,x] )
//
// Might as well form Ns x Nj x Ngamma matrix 
// 
///////////////////////////////////////////////////////////////////////////////

template <typename FImpl>
void TA2APionField<FImpl>::DeltaFeq2(int dt_min,int dt_max,
				     Eigen::Tensor<ComplexD,2> &dF2_fig8,
				     Eigen::Tensor<ComplexD,2> &dF2_trtr,
				     Eigen::Tensor<ComplexD,3> &WW_sd, 
				     const LatticeFermion *vs,
				     const LatticeFermion *vd,
				     int orthogdim)
{
  LOG(Message) << "Computing A2A DeltaF=2 graph" << std::endl;

  int dt = dt_min; // HACK ; should loop over dt

  typedef typename FImpl::SiteSpinor vobj;

  typedef typename vobj::scalar_object sobj;
  typedef typename vobj::scalar_type scalar_type;
  typedef typename vobj::vector_type vector_type;

  typedef iSpinMatrix<vector_type> SpinMatrix_v;
  typedef iSpinMatrix<scalar_type> SpinMatrix_s;
  typedef iSinglet<vector_type> Scalar_v;
  typedef iSinglet<scalar_type> Scalar_s;
  
  int N_s = WW_sd.dimension(1); 
  int N_d = WW_sd.dimension(2);

  GridBase *grid = vs[0]._grid;
  
  const int    nd = grid->_ndimension;
  const int Nsimd = grid->Nsimd();
  int Nt          = grid->GlobalDimensions()[orthogdim];
  
  dF2_trtr.resize(Nt,16);
  dF2_fig8.resize(Nt,16);
  for(int t=0;t<Nt;t++){
    for(int g=0;g<dF2_trtr.dimension(1);g++) dF2_trtr(t,g)= ComplexD(0.0);
    for(int g=0;g<dF2_fig8.dimension(1);g++) dF2_fig8(t,g)= ComplexD(0.0);
  }

  //  std::vector<LatticePropagator> WWVV (Nt,grid)

  int fd=grid->_fdimensions[orthogdim];
  int ld=grid->_ldimensions[orthogdim];
  int rd=grid->_rdimensions[orthogdim];

  //////////////////////////////////////
  // will locally sum vectors first
  // sum across these down to scalars
  // splitting the SIMD
  //////////////////////////////////////
  int MFrvol = rd*N_s*N_d;
  int MFlvol = ld*N_s*N_d;

  Vector<vector_type > lvSum(MFrvol);
  parallel_for (int r = 0; r < MFrvol; r++){
    lvSum[r] = zero;
  }

  Vector<scalar_type > lsSum(MFlvol);             
  parallel_for (int r = 0; r < MFlvol; r++){
    lsSum[r]=scalar_type(0.0);
  }

  int e1=    grid->_slice_nblock[orthogdim];
  int e2=    grid->_slice_block [orthogdim];
  int stride=grid->_slice_stride[orthogdim];

  Eigen::Tensor<Scalar_v,3>  VgV_sd(N_s,N_d,16);  // trace with dirac structure
  Eigen::Tensor<Scalar_s,4>  VgV_sd_l(N_s,N_d,16,Nsimd); 
  int Ng;

  LOG(Message) << "Computing A2A DeltaF=2 graph entering site loop" << std::endl;

  double t_tot   =0;
  double t_vv    =0;
  double t_extr  =0;
  double t_transp=0;
  double t_WW    =0;
  double t_trtr  =0;
  double t_fig8  =0;

  t_tot -=usecond();

  for(int r=0;r<rd;r++){

    LOG(Message) << "Computing A2A DeltaF=2 timeslice "<< r << "/"<<rd << std::endl;
    int so=r*grid->_ostride[orthogdim]; // base offset for start of plane 

    for(int n=0;n<e1;n++){
      for(int b=0;b<e2;b++){

	int ss= so+n*stride+b;

	///////////////////////////////////
	//           _      _
	// O_VV+AA = s V_mu d s V_mu d
	//           _     _
	//         + s A_mu d s A_mu d
	///////////////////////////////////

	// assemble the v_s v_d spin outer, colour inner product matrix
	// the #vecs will be large
	t_vv -=usecond();
	parallel_for(int d=0;d<N_d;d++){

	  auto _vd = conjugate(vd[d]._odata[ss]);

	  for(int s=0;s<N_s;s++){

	    SpinMatrix_v vv;
	    auto _vs = vs[s]._odata[ss];

	    for(int s1=0;s1<Ns;s1++){
	    for(int s2=0;s2<Ns;s2++){
	      vv()(s1,s2)() = _vd()(s2)(0) * _vs()(s1)(0)
		+             _vd()(s2)(1) * _vs()(s1)(1)
		+             _vd()(s2)(2) * _vs()(s1)(2);
	    }}

	    int g=0;
	    //	    VgV_sd(s,d,g++) = trace(vv);      // S
	    //	    VgV_sd(s,d,g++)() = traceGamma5(vv());// P
	    VgV_sd(s,d,g++)() = traceGammaX(vv());// Vmu
	    VgV_sd(s,d,g++)() = traceGammaY(vv());
	    VgV_sd(s,d,g++)() = traceGammaZ(vv());
	    VgV_sd(s,d,g++)() = traceGammaT(vv());

	    VgV_sd(s,d,g++)() = traceGammaXGamma5(vv());// Amu
	    VgV_sd(s,d,g++)() = traceGammaYGamma5(vv());
	    VgV_sd(s,d,g++)() = traceGammaZGamma5(vv());
	    VgV_sd(s,d,g++)() = traceGammaTGamma5(vv());

	    /*
	    VgV_sd(s,d,g++)() = traceSigmaXY(vv);// Sigma_munu unimplemented
	    VgV_sd(s,d,g++)() = traceSigmaXZ(vv);
	    VgV_sd(s,d,g++)() = traceSigmaXT(vv);
	    VgV_sd(s,d,g++)() = traceSigmaYZ(vv);
	    VgV_sd(s,d,g++)() = traceSigmaYT(vv);
	    VgV_sd(s,d,g++)() = traceSigmaZT(vv);
	    */
	    Ng = g;

	  }
	}
	t_vv +=usecond();

	/////////////////////////////////////////////////////////////////////
        // PLAN: Make use of the fact that the trace of the product of two
	// matrices is a_ij b_ji =   a_ij (b^T)_ij and write as a "ZDOT"
	// with one of the matrices transposed.
	//
	// Further since we want to take the trace-product with nt 
	// other matrices can do these in a single pass and gain another 2x.
	// Can do this data parallel on all SIMD lanes of VgV_sd, with a 
	// broadcast of WW[t0] and WW[t1].
	//
	// 
	// Strategies:
	//
	// Wick1:  TR TR(tx) = sum_x [ Trace(WW[t0] VgV(t,x) )  x Trace( WW_[t1] VgV(t,x) ) ]
	//            
	//         Take Tr (WW[t'] . VgV^T) for all "t'", store in array.
	// 
	//         Accumulate (WW[t0].VgV). (WW[t1].VgV) for dt and (t-t0+nt) % nt.
	//           
	// Wick2:  
	//
	// Fig8(tx)  = sum_x Trace( VV[t,x] WW[t0] VV[t,x] WW[t1] )
	//
	//       for(t0) 
	//         form (VV^T WW VV^T)^T = VV WW^T VV = M0
	//         for(t1) 
	//           Accumulate Tr(M0 . WW[t1] in dt and (t-t0+nt)%nt
	//
	/////////////////////////////////////////////////////////////////////
	// Loop over t0 accumulate the dT matrix elements.
	for(int g=0;g<Ng;g++){
	  
	  /////////////////////////////////////////////////
	  // break out into Nsimd sites worth of scalar code.
	  /////////////////////////////////////////////////
	  t_extr -=usecond();
	  parallel_for(int d=0;d<N_d;d++){
	    std::vector<Scalar_s> extracted(Nsimd);               
	    Scalar_v temp; 
	    for(int s=0;s<N_s;s++){
	      for(int l=0;l<Nsimd;l++){

		temp = VgV_sd(s,d,g);
		extract(temp,extracted);
		for(int l=0;l<Nsimd;l++){
		  VgV_sd_l(s,d,g,l) = extracted[l]()()();
		}
		    
	      }// lane
	    }
	  }// s,d
	  t_extr +=usecond();

	  ////////////////////////////////////////////
	  // Work on a series of scalar problems
	  ////////////////////////////////////////////
	  for(int l=0;l<Nsimd;l++){

	    std::vector<int> icoor(nd);
	    grid->iCoorFromIindex(icoor,l);
	    int ttt = r+icoor[orthogdim]*rd;

	    Eigen::MatrixXcd VgV(N_d,N_s);
	    Eigen::MatrixXcd VgV_T(N_s,N_d);
	    Eigen::MatrixXcd WW0(N_s,N_d);
	    Eigen::MatrixXcd WW1(N_s,N_d);

	    /////////////////////////////////////////
	    // Single site VgV , scalar order
	    /////////////////////////////////////////
	    t_transp -=usecond();
	    parallel_for(int d=0;d<N_d;d++){ for(int s=0;s<N_s;s++){
		// Note pre-transpose VgV in the copy out
		VgV(d,s)   = VgV_sd_l(s,d,g,l);
		VgV_T(s,d) = VgV(d,s);
	    }}
	    t_transp +=usecond();
	    
	    /////////////////////////////////////////
	    // loop over time planes of meson
	    /////////////////////////////////////////
	    for(int t0=0;t0<Nt;t0++){	    
	      int t1 = (t0+dt)%Nt;       // Future loop over dT 
	      int tt  = (ttt-t0+Nt)%Nt;  // Time of this site relative to t0
	      /////////////////////////////////////////////////////////////////
	      // Extract this pair of WW matrices infuture loop over dT
	      /////////////////////////////////////////////////////////////////
	      t_WW -=usecond();
	      parallel_for(int d=0;d<N_d;d++){ for(int s=0;s<N_s;s++){
		  WW0(s,d) = WW_sd(t0,s,d); 
		  WW1(s,d) = WW_sd(t1,s,d); 
	      }}
	      t_WW +=usecond();

	      /////////////////////////////////////////
	      //	    Wick1 -- transpose
	      /////////////////////////////////////////
	      //                     VgV_ds WW_sd
	      t_trtr -=usecond();
	      ComplexD trWW0VV = (VgV_T.array()*WW0.array()).sum();
	      ComplexD trWW1VV = (VgV_T.array()*WW1.array()).sum();
	      dF2_trtr(tt,g) +=  trWW0VV * trWW1VV;
	      t_trtr +=usecond();

	      /////////////////////////////////////////
	      //	    Wick2 -- transpose
	      /////////////////////////////////////////
	      //             VgV_ds WW_sd VgV_d's'
	      //
    	      // This is the time consuming loop.
	      t_fig8 -=usecond();
	      Eigen::MatrixXcd VVWW0VV = VgV * WW0 * VgV ;
	      Eigen::MatrixXcd VVWW0VV_T =VVWW0VV.transpose();
	      auto trVVWW0VVWW1 =(VVWW0VV_T.array()*WW1.array()).sum();
	      
	      dF2_fig8(tt,g) += trVVWW0VVWW1;
	      t_fig8 +=usecond();

	    }// t0 loop

	  }// l loop
	  
	}// gamma
	
      }// close loop over time plane and block stride loop within timeplane
    }
  }

  grid->GlobalSumVector(&dF2_fig8(0),Nt*Ng);
  grid->GlobalSumVector(&dF2_trtr(0),Nt*Ng);
  t_tot +=usecond();

  LOG(Message) << "Computing A2A DeltaF=2 graph t_tot    " << t_tot    << " us "<< std::endl;
  LOG(Message) << "Computing A2A DeltaF=2 graph t_vv     " << t_vv     << " us "<< std::endl;
  LOG(Message) << "Computing A2A DeltaF=2 graph t_extr   " << t_extr   << " us "<< std::endl;
  LOG(Message) << "Computing A2A DeltaF=2 graph t_transp " << t_transp << " us "<< std::endl;
  LOG(Message) << "Computing A2A DeltaF=2 graph t_WW     " << t_WW     << " us "<< std::endl;
  LOG(Message) << "Computing A2A DeltaF=2 graph t_trtr   " << t_trtr   << " us "<< std::endl;
  LOG(Message) << "Computing A2A DeltaF=2 graph t_fig8   " << t_fig8   << " us "<< std::endl;

}

// WW is w_s^dag (x) w_d       (G5 implicitly absorbed)
//
// WWVV will have spin-col (x) spin-col tensor.
//
// Want this to look like a strange fwd prop, anti-d prop.
// 
// Take WW_sd v^dag_d (x) v_s
// 
// 

template <typename FImpl>
void TA2APionField<FImpl>::DeltaFeq2_alt(int dt_min,int dt_max,
					 Eigen::Tensor<ComplexD,2> &dF2_fig8,
					 Eigen::Tensor<ComplexD,2> &dF2_trtr,
					 Eigen::Tensor<ComplexD,1> &den0,
					 Eigen::Tensor<ComplexD,1> &den1,
					 Eigen::Tensor<ComplexD,3> &WW_sd, 
					 const LatticeFermion *vs,
					 const LatticeFermion *vd,
					 int orthogdim)
{
  LOG(Message) << "Computing A2A DeltaF=2 graph" << std::endl;

  int dt = dt_min; // HACK ; should loop over dt

  auto G5 = Gamma(Gamma::Algebra::Gamma5);

  typedef typename FImpl::SiteSpinor vobj;

  typedef typename vobj::scalar_object sobj;
  typedef typename vobj::scalar_type scalar_type;
  typedef typename vobj::vector_type vector_type;

  typedef iSpinMatrix<vector_type> SpinMatrix_v;
  typedef iSpinMatrix<scalar_type> SpinMatrix_s;
  typedef iSinglet<vector_type> Scalar_v;
  typedef iSinglet<scalar_type> Scalar_s;
  
  int N_s = WW_sd.dimension(1); 
  int N_d = WW_sd.dimension(2);

  GridBase *grid = vs[0]._grid;
  
  const int    nd = grid->_ndimension;
  const int Nsimd = grid->Nsimd();
  int N_t          = grid->GlobalDimensions()[orthogdim];
  double vol         = 1.0;
  for(int dim=0;dim<nd;dim++){
    vol = vol * grid->GlobalDimensions()[dim];
  }
  double nodes = grid->NodeCount();
  dF2_trtr.resize(N_t,16);
  dF2_fig8.resize(N_t,16);

  den0.resize(N_t);
  den1.resize(N_t);
  for(int t=0;t<N_t;t++){
    for(int g=0;g<dF2_trtr.dimension(1);g++) dF2_trtr(t,g)= ComplexD(0.0);
    for(int g=0;g<dF2_fig8.dimension(1);g++) dF2_fig8(t,g)= ComplexD(0.0);
    den0(t) =ComplexD(0.0);
    den1(t) =ComplexD(0.0);
  }


  LatticeComplex D0(grid); // <P|A0> correlator from each wall
  LatticeComplex D1(grid);

  LatticeComplex O1_trtr(grid);
  LatticeComplex O2_trtr(grid);
  LatticeComplex O3_trtr(grid);
  LatticeComplex O4_trtr(grid);
  LatticeComplex O5_trtr(grid);

  LatticeComplex O1_fig8(grid);
  LatticeComplex O2_fig8(grid);
  LatticeComplex O3_fig8(grid);
  LatticeComplex O4_fig8(grid);
  LatticeComplex O5_fig8(grid);

  O1_trtr = zero;
  O2_trtr = zero;
  O3_trtr = zero;
  O4_trtr = zero;
  O5_trtr = zero;

  O1_fig8 = zero;
  O2_fig8 = zero;
  O3_fig8 = zero;
  O4_fig8 = zero;
  O5_fig8 = zero;

  D0 = zero;
  D1 = zero;

  double t_tot  = -usecond();
  std::vector<LatticePropagator> WWVV (N_t,grid);
  for(int t=0;t<N_t;t++){
    WWVV[t] = zero;
  }

  //////////////////////////////////////////////////////////
  // Easily cache blocked and unrolled.
  //////////////////////////////////////////////////////////
  //
  // Ways to speed up?? 20 min in 8^3 x 16 local volume with 400 modes.
  //
  // Flops are  (cmul=6 + 1 (add))* 12*12  * N_s * N_d
  // Bytes are read Nc x Ns x 2  + read/write Nc^2 x Ns^2 complex 
  //
  double t_outer= -usecond();
  double t_outer_0 =0.0;
  double t_outer_1 =0.0;


#undef  METHOD_1
#define METHOD_5

#ifdef METHOD_1
  LOG(Message) << "METHOD_1" << std::endl;
  // Method-1
  // i) calculate flop rate and data rate. 17 GF/s and 80 GB/s
  //    Dominated by accumulating the N_t summands
  parallel_for(int ss=0;ss<grid->oSites();ss++){
    for(int s=0;s<N_s;s++){
      for(int d=0;d<N_d;d++){

	double t0=-usecond();
	// RHS is conjugated
	auto tmp = outerProduct(vs[s]._odata[ss],vd[d]._odata[ss]); 
	t0+=usecond();

	double t1=-usecond();
	for(int t=0;t<N_t;t++) { 
	  WWVV[t]._odata[ss] += WW_sd(t,s,d) * tmp;
	}
	t1+=usecond();

#pragma omp critical
	{ t_outer_0+=t0; t_outer_1+=t1;}
	
      }
    }
  }
#endif

  double b_outer = (12+12+2*12*12 ) * N_s * N_d * vol * sizeof(Complex);
  double f_outer = 8 * 12*12 * N_s * N_d * vol;

#ifdef METHOD_5
  LOG(Message) << "METHOD_5" << std::endl;
  // Method-5
  parallel_for(int ss=0;ss<grid->oSites();ss++){
    for(int d_o=0;d_o<N_d;d_o+=d_unroll){
      for(int t=0;t<N_t;t++){
      for(int s=0;s<N_s;s++){
	auto tmp1 = vs[s]._odata[ss];
	vobj tmp2 = zero;
	// Surprisingly slow
	for(int d=d_o;d<MIN(d_o+d_unroll,N_d);d++){
	  Scalar_v coeff = WW_sd(t,s,d);
	  mac(&tmp2 ,& coeff, & vd[d]._odata[ss]);
	}

	// Outer product of tmp1 with a sum of terms suppressed by d_unroll
	tmp2 = conjugate(tmp2);
	for(int s1=0;s1<Ns;s1++){
	for(int s2=0;s2<Ns;s2++){
	  WWVV[t]._odata[ss]()(s1,s2)(0,0) += tmp1()(s1)(0)*tmp2()(s2)(0);
	  WWVV[t]._odata[ss]()(s1,s2)(0,1) += tmp1()(s1)(0)*tmp2()(s2)(1);
	  WWVV[t]._odata[ss]()(s1,s2)(0,2) += tmp1()(s1)(0)*tmp2()(s2)(2);
	  WWVV[t]._odata[ss]()(s1,s2)(1,0) += tmp1()(s1)(1)*tmp2()(s2)(0);
	  WWVV[t]._odata[ss]()(s1,s2)(1,1) += tmp1()(s1)(1)*tmp2()(s2)(1);
	  WWVV[t]._odata[ss]()(s1,s2)(1,2) += tmp1()(s1)(1)*tmp2()(s2)(2);
	  WWVV[t]._odata[ss]()(s1,s2)(2,0) += tmp1()(s1)(2)*tmp2()(s2)(0);
	  WWVV[t]._odata[ss]()(s1,s2)(2,1) += tmp1()(s1)(2)*tmp2()(s2)(1);
	  WWVV[t]._odata[ss]()(s1,s2)(2,2) += tmp1()(s1)(2)*tmp2()(s2)(2);
	}}

      }}
    }
  }
#endif

  t_outer+=usecond();

  //////////////////////////////
  // Implicit gamma-5 
  //////////////////////////////
  for(int t=0;t<N_t;t++){
    WWVV[t] = WWVV[t]* G5 ; 
  }


  //////////////////////////////////////////////////
  // Used to store appropriate correlation funcs
  //////////////////////////////////////////////////
  std::vector<TComplex>  C1;
  std::vector<TComplex>  C2;
  std::vector<TComplex>  C3;
  std::vector<TComplex>  C4;
  std::vector<TComplex>  C5;

  //////////////////////////////////////////////////////////
  // Could do AA, VV, SS, PP, TT and form linear combinations later.
  // Almost 2x. but for many modes, the above loop dominates.
  //////////////////////////////////////////////////////////
  double t_contr= -usecond();
  for(int t0=0;t0<N_t;t0++){
    // No loop over t1
    // Cost is trivial to add given cost of outer product
    {
      int t1 = (t0+dt)%N_t;

      parallel_for(int ss=0;ss<grid->oSites();ss++){

	auto VV0= WWVV[t0]._odata[ss];
	auto VV1= WWVV[t1]._odata[ss];

	// Tr Tr Wick contraction
	auto VX = Gamma(Gamma::Algebra::GammaX);
	auto VY = Gamma(Gamma::Algebra::GammaY);
	auto VZ = Gamma(Gamma::Algebra::GammaZ);
	auto VT = Gamma(Gamma::Algebra::GammaT);

	auto AX = Gamma(Gamma::Algebra::GammaXGamma5);
	auto AY = Gamma(Gamma::Algebra::GammaYGamma5);
	auto AZ = Gamma(Gamma::Algebra::GammaZGamma5);
	auto AT = Gamma(Gamma::Algebra::GammaTGamma5);

	auto S  = Gamma(Gamma::Algebra::Identity);
	auto P  = Gamma(Gamma::Algebra::Gamma5);

	auto T0 = Gamma(Gamma::Algebra::SigmaXY);
	auto T1 = Gamma(Gamma::Algebra::SigmaXZ);
	auto T2 = Gamma(Gamma::Algebra::SigmaXT);
	auto T3 = Gamma(Gamma::Algebra::SigmaYZ);
	auto T4 = Gamma(Gamma::Algebra::SigmaYT);
	auto T5 = Gamma(Gamma::Algebra::SigmaZT);

	O1_trtr._odata[ss] = trace(VX*VV0) * trace(VX*VV1)
	  +                  trace(VY*VV0) * trace(VY*VV1)
	  +                  trace(VZ*VV0) * trace(VZ*VV1)
	  +                  trace(VT*VV0) * trace(VT*VV1)
          +                  trace(AX*VV0) * trace(AX*VV1)
	  +                  trace(AY*VV0) * trace(AY*VV1)
	  +                  trace(AZ*VV0) * trace(AZ*VV1)
	  +                  trace(AT*VV0) * trace(AT*VV1);

	O2_trtr._odata[ss] = trace(VX*VV0) * trace(VX*VV1)
	  +                  trace(VY*VV0) * trace(VY*VV1)
	  +                  trace(VZ*VV0) * trace(VZ*VV1)
	  +                  trace(VT*VV0) * trace(VT*VV1)
          -                  trace(AX*VV0) * trace(AX*VV1)
	  -                  trace(AY*VV0) * trace(AY*VV1)
	  -                  trace(AZ*VV0) * trace(AZ*VV1)
	  -                  trace(AT*VV0) * trace(AT*VV1);

	O3_trtr._odata[ss] = trace(S*VV0) * trace(S*VV1)
	  +                  trace(P*VV0) * trace(P*VV1);

	O4_trtr._odata[ss] = trace(S*VV0) * trace(S*VV1)
	  -                  trace(P*VV0) * trace(P*VV1);

	O5_trtr._odata[ss] = trace(T0*VV0) * trace(T0*VV1)
	  +                  trace(T1*VV0) * trace(T1*VV1)
	  +                  trace(T2*VV0) * trace(T2*VV1)
	  +                  trace(T3*VV0) * trace(T3*VV1)
	  +                  trace(T4*VV0) * trace(T4*VV1)
	  +                  trace(T5*VV0) * trace(T5*VV1);


	////////////////////////////////////
	// Fig8 Wick contraction
	////////////////////////////////////
        // was (UKhadron) Wick2 = trace( Qs0 * Gamma(G5) * anti_d0 * G * Qs1 * Gamma(G5) * anti_d1 * G );
	//
	//                      = trace( [ Vs WW_sd'[t0] Vd'^dag ]  G [ Vs' WW_s'd[t1] Vd^dag ] G  )
	//                      = trace( WWVV * G * WWVV * G )
	// 
	// Can do VV and AA seperately and then add/sub later cheaper

	O1_fig8._odata[ss]  = trace (VV0 * VX * VV1 * VX)
	  +                   trace (VV0 * VY * VV1 * VY)
	  +                   trace (VV0 * VZ * VV1 * VZ)
	  +                   trace (VV0 * VT * VV1 * VT)
          +                   trace (VV0 * AX * VV1 * AX)
	  +                   trace (VV0 * AY * VV1 * AY)
	  +                   trace (VV0 * AZ * VV1 * AZ)
	  +                   trace (VV0 * AT * VV1 * AT);

	O2_fig8._odata[ss]  = trace (VV0 * VX * VV1 * VX)
	  +                   trace (VV0 * VY * VV1 * VY)
	  +                   trace (VV0 * VZ * VV1 * VZ)
	  +                   trace (VV0 * VT * VV1 * VT)
          -                   trace (VV0 * AX * VV1 * AX)
	  -                   trace (VV0 * AY * VV1 * AY)
	  -                   trace (VV0 * AZ * VV1 * AZ)
	  -                   trace (VV0 * AT * VV1 * AT);

	O3_fig8._odata[ss]  = trace (VV0 * S * VV1 * S)
	  +                   trace (VV0 * P * VV1 * P);

	O4_fig8._odata[ss]  = trace (VV0 * S * VV1 * S)
	  -                   trace (VV0 * P * VV1 * P);

	O5_fig8._odata[ss] =  trace (VV0 * T0 * VV1 * T0)
	  +                   trace (VV0 * T1 * VV1 * T1)
	  +                   trace (VV0 * T2 * VV1 * T2)
	  +                   trace (VV0 * T3 * VV1 * T3)
	  +                   trace (VV0 * T4 * VV1 * T4)
	  +                   trace (VV0 * T5 * VV1 * T5);


	// Hack force PP correlator
	//	D0._odata[ss] = trace(P*VV0); // These match signed off
	//	D1._odata[ss] = trace(P*VV1);
	D0._odata[ss] = trace(AT*VV0);
	D1._odata[ss] = trace(AT*VV1);

      }

      sliceSum(O1_trtr,C1, orthogdim);
      sliceSum(O2_trtr,C2, orthogdim);
      sliceSum(O3_trtr,C3, orthogdim);
      sliceSum(O4_trtr,C4, orthogdim);
      sliceSum(O5_trtr,C5, orthogdim);

      for(int t=0;t<N_t;t++){// 2x from Wick contraction reordering
	dF2_trtr(t,0)+= 2.0*C1[(t+t0)%N_t]()()()/vol;
	dF2_trtr(t,1)+= 2.0*C2[(t+t0)%N_t]()()()/vol;
	dF2_trtr(t,2)+= 2.0*C3[(t+t0)%N_t]()()()/vol;
	dF2_trtr(t,3)+= 2.0*C4[(t+t0)%N_t]()()()/vol;
	dF2_trtr(t,4)+= 2.0*C5[(t+t0)%N_t]()()()/vol;
      }

      sliceSum(O1_fig8,C1, orthogdim);
      sliceSum(O2_fig8,C2, orthogdim);
      sliceSum(O3_fig8,C3, orthogdim);
      sliceSum(O4_fig8,C4, orthogdim);
      sliceSum(O5_fig8,C5, orthogdim);

      for(int t=0;t<N_t;t++){
	dF2_fig8(t,0)= 2.0*C1[(t+t0)%N_t]()()()/vol;
	dF2_fig8(t,1)= 2.0*C2[(t+t0)%N_t]()()()/vol;
	dF2_fig8(t,2)= 2.0*C3[(t+t0)%N_t]()()()/vol;
	dF2_fig8(t,3)= 2.0*C4[(t+t0)%N_t]()()()/vol;
	dF2_fig8(t,4)= 2.0*C5[(t+t0)%N_t]()()()/vol;
      }

      sliceSum(D0,C1, orthogdim);
      sliceSum(D1,C2, orthogdim);
      for(int t=0;t<N_t;t++){
	den0(t)+=C1[(t+t0)%N_t]()()()/vol;
	den1(t)+=C2[(t+t0)%N_t]()()()/vol;
      }

    }
  }
  t_contr +=usecond();
  
  t_tot+=usecond();
  double million=1.0e6;
  LOG(Message) << "Computing A2A DeltaF=2 graph t_tot      " << t_tot      /million << " s "<< std::endl;
  LOG(Message) << "Computing A2A DeltaF=2 graph t_outer    " << t_outer    /million << " s "<< std::endl;
  LOG(Message) << "Computing A2A DeltaF=2 graph  t_outer_0 " << t_outer_0  /million << " s "<< std::endl;
  LOG(Message) << "Computing A2A DeltaF=2 graph  t_outer_1 " << t_outer_1  /million << " s "<< std::endl;
  LOG(Message) << "Computing A2A DeltaF=2 graph t_contr    " << t_contr    /million << " s "<< std::endl;

  LOG(Message) << "Computing A2A DeltaF=2 graph mflops/s outer  " << f_outer/t_outer << " MF/s "<< std::endl;
  LOG(Message) << "Computing A2A DeltaF=2 graph MB/s     outer  " << b_outer/t_outer << " MB/s "<< std::endl;
  LOG(Message) << "Computing A2A DeltaF=2 graph mflops/s outer/node  " << f_outer/t_outer/nodes << " MF/s "<< std::endl;
  LOG(Message) << "Computing A2A DeltaF=2 graph MB/s     outer/node  " << b_outer/t_outer/nodes << " MB/s "<< std::endl;

}



// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2APionField<FImpl>::execute(void)
{
    LOG(Message) << "Computing A2A Pion fields" << std::endl;

    auto &a2a_i = envGet(A2ABase, par().A2A_i + "_class");
    auto &a2a_j = envGet(A2ABase, par().A2A_j + "_class");

    ///////////////////////////////////////////////
    // Square assumption for now Nl = Nr = N
    ///////////////////////////////////////////////
    int nt = env().getDim(Tp);
    int nx = env().getDim(Xp);
    int ny = env().getDim(Yp);
    int nz = env().getDim(Zp);

    //    int N_i  = a2a_i.par().N;
    //    int N_j  = a2a_j.par().N;
    int N_i  = a2a_i.getN();
    int N_j  = a2a_j.getN();

    int nmom=par().Nmom;

    int schurBlock = par().schurBlock;
    int cacheBlock = par().cacheBlock;


    ///////////////////////////////////////////////
    // Momentum setup
    ///////////////////////////////////////////////
    GridBase *grid = env().getGrid(1);
    std::vector<LatticeComplex> phases(nmom,grid);
    for(int m=0;m<nmom;m++){
      phases[m] = Complex(1.0);    // All zero momentum for now
    }

    ///////////////////////////////////////////////////////////////////////
    // i and j represent different flavours, hits, with different ranks.
    // in general non-square case.
    ///////////////////////////////////////////////////////////////////////
    Eigen::Tensor<ComplexD,4> pionFieldWVmom_ij     (nmom,nt,N_i,N_j);    
    Eigen::Tensor<ComplexD,3> pionFieldWV_ij        (nt,N_i,N_j);    

    Eigen::Tensor<ComplexD,4> pionFieldWVmom_ji     (nmom,nt,N_j,N_i);    
    Eigen::Tensor<ComplexD,3> pionFieldWV_ji        (nt,N_j,N_i);    


    LOG(Message) << "Rank for A2A PionField is " << N_i << " x "<<N_j << std::endl;

    envGetTmp(std::vector<FermionField>, wi);
    envGetTmp(std::vector<FermionField>, vi);

    envGetTmp(std::vector<FermionField>, wj);
    envGetTmp(std::vector<FermionField>, vj);
    envGetTmp(FermionField, tmp_5d);

    LOG(Message) << "Finding v and w vectors " << std::endl;

    //////////////////////////////////////////////////////////////////////////
    // i,j   is first  loop over SchurBlock factors reusing 5D matrices
    // ii,jj is second loop over cacheBlock factors for high perf contractoin
    // iii,jjj are loops within cacheBlock
    // Total index is sum of these  i+ii+iii etc...
    //////////////////////////////////////////////////////////////////////////
    
    double flops = 0.0;
    double bytes = 0.0;
    double vol   = nx*ny*nz*nt;
    double vol3  = nx*ny*nz;
    double t_schur=0;
    double t_contr_vwm=0;
    double t_contr_vw=0;
    double t_contr_ww=0;
    double t_contr_vv=0;

    double tt0 = usecond();
    for(int i=0;i<N_i;i+=schurBlock){ //loop over SchurBlocking to suppress 5D matrix overhead
    for(int j=0;j<N_j;j+=schurBlock){

      ///////////////////////////////////////////////////////////////
      // Get the W and V vectors for this schurBlock^2 set of terms
      ///////////////////////////////////////////////////////////////
      int N_ii = MIN(N_i-i,schurBlock);
      int N_jj = MIN(N_j-j,schurBlock);

      t_schur-=usecond();
      for(int ii =0;ii < N_ii;ii++) a2a_i.return_w(i+ii, tmp_5d, wi[ii]);
      for(int jj =0;jj < N_jj;jj++) a2a_j.return_w(j+jj, tmp_5d, wj[jj]);

      for(int ii =0;ii < N_ii;ii++) a2a_i.return_v(i+ii, tmp_5d, vi[ii]);
      for(int jj =0;jj < N_jj;jj++) a2a_j.return_v(j+jj, tmp_5d, vj[jj]);
      t_schur+=usecond();

      LOG(Message) << "Found i w&v vectors " << i <<" .. " << i+N_ii-1 << std::endl;
      LOG(Message) << "Found j w&v vectors " << j <<" .. " << j+N_jj-1 << std::endl;

      ///////////////////////////////////////////////////////////////
      // Series of cache blocked chunks of the contractions within this SchurBlock
      /////////////////////////////////////////////////////////////// 
      for(int ii=0;ii<N_ii;ii+=cacheBlock){
      for(int jj=0;jj<N_jj;jj+=cacheBlock){

	int N_iii = MIN(N_ii-ii,cacheBlock);
	int N_jjj = MIN(N_jj-jj,cacheBlock);

	Eigen::Tensor<ComplexD,4> pionFieldWVmomB_ij(nmom,nt,N_iii,N_jjj);    
	Eigen::Tensor<ComplexD,4> pionFieldWVmomB_ji(nmom,nt,N_jjj,N_iii);    

	Eigen::Tensor<ComplexD,3> pionFieldWVB_ij(nt,N_iii,N_jjj);    
	Eigen::Tensor<ComplexD,3> pionFieldWVB_ji(nt,N_jjj,N_iii);    

	t_contr_vwm-=usecond();
	PionFieldWVmom(pionFieldWVmomB_ij, &wi[ii], &vj[jj], phases,Tp);
	PionFieldWVmom(pionFieldWVmomB_ji, &wj[jj], &vi[ii], phases,Tp);
	t_contr_vwm+=usecond();

	t_contr_vw-=usecond();
	PionFieldWV(pionFieldWVB_ij, &wi[ii], &vj[jj],Tp);
	PionFieldWV(pionFieldWVB_ji, &wj[jj], &vi[ii],Tp);
	t_contr_vw+=usecond();


	flops += vol * ( 2 * 8.0 + 6.0 + 8.0*nmom) * N_iii*N_jjj;

	bytes  += vol * (12.0 * sizeof(Complex) ) * N_iii*N_jjj
	       +  vol * ( 2.0 * sizeof(Complex) *nmom ) * N_iii*N_jjj;

	///////////////////////////////////////////////////////////////
	// Copy back to full meson field tensor
	/////////////////////////////////////////////////////////////// 
	parallel_for_nest2(int iii=0;iii< N_iii;iii++) {
        for(int jjj=0;jjj< N_jjj;jjj++) {

	  for(int m =0;m< nmom;m++) {
          for(int t =0;t< nt;t++) {
	    pionFieldWVmom_ij(m,t,i+ii+iii,j+jj+jjj) = pionFieldWVmomB_ij(m,t,iii,jjj);
	    pionFieldWVmom_ji(m,t,j+jj+jjj,i+ii+iii) = pionFieldWVmomB_ji(m,t,jjj,iii);
	  }}

          for(int t =0;t< nt;t++) {
	    pionFieldWV_ij(t,i+ii+iii,j+jj+jjj) = pionFieldWVB_ij(t,iii,jjj);
	    pionFieldWV_ji(t,j+jj+jjj,i+ii+iii) = pionFieldWVB_ji(t,jjj,iii);
	  }


	}}
      }}
    }}

    double nodes=grid->NodeCount();
    double tt1 = usecond();
    LOG(Message) << " Contraction of PionFields took "<<(tt1-tt0)/1.0e6<< " seconds "  << std::endl;
    LOG(Message) << " Schur "<<(t_schur)/1.0e6<< " seconds "  << std::endl;
    LOG(Message) << " Contr WVmom "<<(t_contr_vwm)/1.0e6<< " seconds "  << std::endl;
    LOG(Message) << " Contr WV    "<<(t_contr_vw)/1.0e6<< " seconds "  << std::endl;

    double t_kernel = t_contr_vwm;
    LOG(Message) << " Arith "<<flops/(t_kernel)/1.0e3/nodes<< " Gflop/s / node "  << std::endl;
    LOG(Message) << " Arith "<<bytes/(t_kernel)/1.0e3/nodes<< " GB/s /node "  << std::endl;

    /////////////////////////////////////////////////////////////////////////
    // Test: Build the pion correlator (two end)
    // < PI_ij(t0) PI_ji (t0+t) >
    /////////////////////////////////////////////////////////////////////////
    std::vector<ComplexD> corrMom(nt,ComplexD(0.0));

    for(int i=0;i<N_i;i++){
    for(int j=0;j<N_j;j++){
      int m=0; // first momentum
      for(int t0=0;t0<nt;t0++){
      for(int t=0;t<nt;t++){
	int tt = (t0+t)%nt;
	corrMom[t] += pionFieldWVmom_ij(m,t0,i,j)* pionFieldWVmom_ji(m,tt,j,i);
      }}
    }}    
    for(int t=0;t<nt;t++) corrMom[t] = corrMom[t]/ (double)nt;

    for(int t=0;t<nt;t++) LOG(Message) << " C_vwm " << t << " " << corrMom[t]<<std::endl;


    /////////////////////////////////////////////////////////////////////////
    // Test: Build the pion correlator (two end) from zero mom contraction
    // < PI_ij(t0) PI_ji (t0+t) >
    /////////////////////////////////////////////////////////////////////////
    std::vector<ComplexD> corr(nt,ComplexD(0.0));

    for(int i=0;i<N_i;i++){
    for(int j=0;j<N_j;j++){
      for(int t0=0;t0<nt;t0++){
      for(int t=0;t<nt;t++){
	int tt = (t0+t)%nt;
	corr[t] += pionFieldWV_ij(t0,i,j)* pionFieldWV_ji(tt,j,i);
      }}
    }}    
    for(int t=0;t<nt;t++) corr[t] = corr[t]/ (double)nt;

    for(int t=0;t<nt;t++) LOG(Message) << " C_vw " << t << " " << corr[t]<<std::endl;


    /////////////////////////////////////////////////////////////////////////
    // Test: Build the pion correlator from zero mom contraction with revers 
    // charge flow
    /////////////////////////////////////////////////////////////////////////
    std::vector<ComplexD> corr_wwvv(nt,ComplexD(0.0));

    wi.resize(N_i,grid);
    vi.resize(N_i,grid);
    wj.resize(N_j,grid);
    vj.resize(N_j,grid);

    for(int i =0;i < N_i;i++) a2a_i.return_v(i, tmp_5d, vi[i]);
    for(int i =0;i < N_i;i++) a2a_i.return_w(i, tmp_5d, wi[i]);
    for(int j =0;j < N_j;j++) a2a_j.return_v(j, tmp_5d, vj[j]);
    for(int j =0;j < N_j;j++) a2a_j.return_w(j, tmp_5d, wj[j]);

    Eigen::Tensor<ComplexD,3> pionFieldWW_ij        (nt,N_i,N_j);    
    Eigen::Tensor<ComplexD,3> pionFieldVV_ji        (nt,N_j,N_i);    
    Eigen::Tensor<ComplexD,3> pionFieldWW_ji        (nt,N_j,N_i);    
    Eigen::Tensor<ComplexD,3> pionFieldVV_ij        (nt,N_i,N_j);    

    PionFieldWW(pionFieldWW_ij, &wi[0], &wj[0],Tp);
    PionFieldVV(pionFieldVV_ji, &vj[0], &vi[0],Tp);
    PionFieldWW(pionFieldWW_ji, &wj[0], &wi[0],Tp);
    PionFieldVV(pionFieldVV_ij, &vi[0], &vj[0],Tp);


    for(int i=0;i<N_i;i++){
    for(int j=0;j<N_j;j++){
      for(int t0=0;t0<nt;t0++){
      for(int t=0;t<nt;t++){
	int tt = (t0+t)%nt;
	corr_wwvv[t] += pionFieldWW_ij(t0,i,j)* pionFieldVV_ji(tt,j,i);
	corr_wwvv[t] += pionFieldWW_ji(t0,j,i)* pionFieldVV_ij(tt,i,j);
      }}
    }}    
    for(int t=0;t<nt;t++) corr_wwvv[t] = corr_wwvv[t] / vol /2.0 ; // (ij+ji noise contribs if i!=j ).

    for(int t=0;t<nt;t++) LOG(Message) << " C_wwvv " << t << " " << corr_wwvv[t]<<std::endl;


    /////////////////////////////////////////////////////////////////////////
    // This is only correct if there are NO low modes
    // Use the "ii" case to construct possible Z wall source one end trick
    /////////////////////////////////////////////////////////////////////////
    std::vector<ComplexD> corr_z2(nt,ComplexD(0.0));
    Eigen::Tensor<ComplexD,3> pionFieldWW        (nt,N_i,N_i);    
    Eigen::Tensor<ComplexD,3> pionFieldVV        (nt,N_i,N_i);    


    PionFieldWW(pionFieldWW, &wi[0], &wi[0],Tp);
    PionFieldVV(pionFieldVV, &vi[0], &vi[0],Tp);
    for(int i=0;i<N_i;i++){
      for(int t0=0;t0<nt;t0++){
      for(int t=0;t<nt;t++){
	int tt = (t0+t)%nt;
	corr_z2[t] += pionFieldWW(t0,i,i) * pionFieldVV(tt,i,i) /vol ;
      }}
    }

    LOG(Message) << " C_z2 WARNING only correct if Nl == 0 "<<std::endl;
    for(int t=0;t<nt;t++) LOG(Message) << " C_z2 " << t << " " << corr_z2[t]<<std::endl;

    /////////////////////////////////////////////////////////////////////////
    // Test: Build a bag contraction
    /////////////////////////////////////////////////////////////////////////
    Eigen::Tensor<ComplexD,2> DeltaF2_fig8  (nt,16);
    Eigen::Tensor<ComplexD,2> DeltaF2_trtr  (nt,16);
    Eigen::Tensor<ComplexD,1> denom0 (nt);
    Eigen::Tensor<ComplexD,1> denom1 (nt);
    
    const int dT=16;

    DeltaFeq2_alt  (dT,dT,DeltaF2_fig8,DeltaF2_trtr,
		    denom0,denom1,
		    pionFieldWW_ij,&vi[0],&vj[0],Tp);
    for(int t=0;t<nt;t++) LOG(Message) << " denom0 [" << t << "]  " << denom0(t)<<std::endl;
    for(int t=0;t<nt;t++) LOG(Message) << " denom1 [" << t << "]  " << denom1(t)<<std::endl;
    for(int g=0;g<4;g++){
      for(int t=0;t<nt;t++) LOG(Message) << " DeltaF2_fig8 [" << t << ","<<g<<"]  " << DeltaF2_fig8(t,g)<<std::endl;
      for(int t=0;t<nt;t++) LOG(Message) << " DeltaF2_trtr [" << t << ","<<g<<"]  " << DeltaF2_trtr(t,g)<<std::endl;
    }
    for(int g=0;g<4;g++){
      for(int t=0;t<nt;t++)
	LOG(Message) << " Bag [" << t << ","<<g<<"]  " 
		     << (DeltaF2_fig8(t,g)+DeltaF2_trtr(t,g)) 
	  /             ( 8.0/3.0 * denom0[t]*denom1[t])
		     <<std::endl;
    }

    /////////////////////////////////////////////////////////////////////////
    // Test: Build a bag contraction the Z2 way
    // Build a wall bag comparison assuming no low modes
    /////////////////////////////////////////////////////////////////////////
    LOG(Message) << " Bag_z2 WARNING only correct if Nl == 0 "<<std::endl;

    int t0=0;
    int t1=dT;
    int Nl=0;
    LatticePropagator Qd0(grid);
    LatticePropagator Qd1(grid);
    LatticePropagator Qs0(grid);
    LatticePropagator Qs1(grid);
    for(int s=0;s<4;s++){
      for(int c=0;c<3;c++){
	int idx0 = Nl+t0*12+s*3+c;
	int idx1 = Nl+t1*12+s*3+c;
	FermToProp<FImpl>(Qd0, vi[idx0], s, c);
	FermToProp<FImpl>(Qd1, vi[idx1], s, c);
	FermToProp<FImpl>(Qs0, vj[idx0], s, c);
	FermToProp<FImpl>(Qs1, vj[idx1], s, c);
      }
    }

    std::vector<Gamma::Algebra> gammas ( {
	  Gamma::Algebra::GammaX,
	  Gamma::Algebra::GammaY,
	  Gamma::Algebra::GammaZ,
	  Gamma::Algebra::GammaT,
	  Gamma::Algebra::GammaXGamma5,
	  Gamma::Algebra::GammaYGamma5,
	  Gamma::Algebra::GammaZGamma5,
	  Gamma::Algebra::GammaTGamma5,
	  Gamma::Algebra::Identity,    
          Gamma::Algebra::Gamma5,
	  Gamma::Algebra::SigmaXY,
	  Gamma::Algebra::SigmaXZ,
	  Gamma::Algebra::SigmaXT,
	  Gamma::Algebra::SigmaYZ,
	  Gamma::Algebra::SigmaYT,
	  Gamma::Algebra::SigmaZT
    });

    auto G5 = Gamma::Algebra::Gamma5;
    LatticePropagator anti_d0 =  adj( Gamma(G5) * Qd0 * Gamma(G5));
    LatticePropagator anti_d1 =  adj( Gamma(G5) * Qd1 * Gamma(G5));
    LatticeComplex TR1(grid);
    LatticeComplex TR2(grid);
    LatticeComplex Wick1(grid);
    LatticeComplex Wick2(grid);

    LatticePropagator PR1(grid);
    LatticePropagator PR2(grid);
    PR1 = Qs0 * Gamma(G5) * anti_d0;
    PR2 = Qs1 * Gamma(G5) * anti_d1;

    for(int g=0;g<Nd*Nd;g++){
      auto g1 = gammas[g];
      Gamma G1 (g1);
      TR1 = trace( PR1 * G1 );
      TR2 = trace( PR2 * G1 );
      Wick1 = TR1*TR2;
      Wick2 = trace( PR1* G1 * PR2 * G1 );
      
      std::vector<TComplex>  C1;
      std::vector<TComplex>  C2;
      std::vector<TComplex>  C3;
      sliceSum(Wick1,C1, Tp);
      sliceSum(Wick2,C2, Tp);
      sliceSum(TR1  ,C3, Tp);
      
      if(g<5){
	for(int t=0;t<C1.size();t++){
	  LOG(Message) << " Wick1["<<g<<","<<t<< "] "<< C1[t]<<std::endl; 
	}
	for(int t=0;t<C2.size();t++){
	  LOG(Message) << " Wick2["<<g<<","<<t<< "] "<< C2[t]<<std::endl; 
	}
      }
      if( (g==9) || (g==7) ){
	for(int t=0;t<C3.size();t++){
	  LOG(Message) << " <G|P>["<<g<<","<<t<< "] "<< C3[t]<<std::endl; 
	}
      } 
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_A2APionField_hpp_
  
