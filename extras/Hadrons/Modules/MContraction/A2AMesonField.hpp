#ifndef Hadrons_MContraction_A2AMesonField_hpp_
#define Hadrons_MContraction_A2AMesonField_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>
#include <Grid/Hadrons/AllToAllVectors.hpp>

#include <unsupported/Eigen/CXX11/Tensor>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         A2AMesonField                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

typedef std::pair<Gamma::Algebra, Gamma::Algebra> GammaPair;


class A2AMesonFieldPar : Serializable
{
  public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2AMesonFieldPar,
				    int, cacheBlock,
				    int, schurBlock,
				    int, N,
				    int, Nl,
                                    std::string, A2A,
                                    std::string, output);
};

template <typename FImpl>
class TA2AMesonField : public Module<A2AMesonFieldPar>
{
  public:
    FERM_TYPE_ALIASES(FImpl, );
    SOLVER_TYPE_ALIASES(FImpl, );

    typedef A2AModesSchurDiagTwo<typename FImpl::FermionField, FMat, Solver> A2ABase;

  public:
    // constructor
    TA2AMesonField(const std::string name);
    // destructor
    virtual ~TA2AMesonField(void){};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);

    // Arithmetic help. Move to Grid??
    virtual void MesonField(Eigen::Tensor<ComplexD,5> &mat, 
			    const std::vector<LatticeFermion > &lhs,
			    const std::vector<LatticeFermion > &rhs,
			    std::vector<Gamma::Algebra> gammas,
			    const std::vector<LatticeComplex > &mom,
			    int orthogdim) ;

};

MODULE_REGISTER(A2AMesonField, ARG(TA2AMesonField<FIMPL>), MContraction);
MODULE_REGISTER(ZA2AMesonField, ARG(TA2AMesonField<ZFIMPL>), MContraction);

/******************************************************************************
*                  TA2AMesonField implementation                             *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TA2AMesonField<FImpl>::TA2AMesonField(const std::string name)
    : Module<A2AMesonFieldPar>(name)
{
}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TA2AMesonField<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().A2A + "_class"};
    in.push_back(par().A2A + "_w_high_4d");
    in.push_back(par().A2A + "_v_high_4d");

    return in;
}

template <typename FImpl>
std::vector<std::string> TA2AMesonField<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {};

    return out;
}


// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2AMesonField<FImpl>::setup(void)
{
    auto &a2a = envGet(A2ABase, par().A2A + "_class");
    int nt = env().getDim(Tp);
    int Nl = par().Nl;
    int N  = par().N;
    int Ls_ = env().getObjectLs(par().A2A + "_class");

    // Four D fields
    envTmp(std::vector<FermionField>, "w", 1, par().schurBlock, FermionField(env().getGrid(1)));
    envTmp(std::vector<FermionField>, "v", 1, par().schurBlock, FermionField(env().getGrid(1)));

    // 5D tmp
    envTmpLat(FermionField, "tmp_5d", Ls_);
}


//////////////////////////////////////////////////////////////////////////////////
// Cache blocked arithmetic routine
// Could move to Grid ???
//////////////////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2AMesonField<FImpl>::MesonField(Eigen::Tensor<ComplexD,5> &mat, 
				       const std::vector<LatticeFermion > &lhs,
				       const std::vector<LatticeFermion > &rhs,
				       std::vector<Gamma::Algebra> gammas,
				       const std::vector<LatticeComplex > &mom,
				       int orthogdim) 
{
  typedef typename FImpl::SiteSpinor vobj;

  typedef typename vobj::scalar_object sobj;
  typedef typename vobj::scalar_type scalar_type;
  typedef typename vobj::vector_type vector_type;

  typedef iSpinMatrix<vector_type> SpinMatrix_v;
  typedef iSpinMatrix<scalar_type> SpinMatrix_s;
  
  int Lblock = lhs.size();
  int Rblock = rhs.size();

  GridBase *grid = lhs[0]._grid;
  
  const int    Nd = grid->_ndimension;
  const int Nsimd = grid->Nsimd();

  int Nt     = grid->GlobalDimensions()[orthogdim];
  int Ngamma = gammas.size();
  int Nmom   = mom.size();

  int fd=grid->_fdimensions[orthogdim];
  int ld=grid->_ldimensions[orthogdim];
  int rd=grid->_rdimensions[orthogdim];

  // will locally sum vectors first
  // sum across these down to scalars
  // splitting the SIMD
  int MFrvol = rd*Lblock*Rblock*Nmom;
  int MFlvol = ld*Lblock*Rblock*Nmom;

  Vector<SpinMatrix_v > lvSum(MFrvol);
  parallel_for (int r = 0; r < MFrvol; r++){
    lvSum[r] = zero;
  }

  Vector<SpinMatrix_s > lsSum(MFlvol);             
  parallel_for (int r = 0; r < MFlvol; r++){
    lsSum[r]=scalar_type(0.0);
  }

  int e1=    grid->_slice_nblock[orthogdim];
  int e2=    grid->_slice_block [orthogdim];
  int stride=grid->_slice_stride[orthogdim];
  
  std::cout << GridLogMessage << " Entering first parallel loop "<<std::endl;

  // Parallelise over t-direction doesn't expose as much parallelism as needed for KNL
  parallel_for(int r=0;r<rd;r++){

    int so=r*grid->_ostride[orthogdim]; // base offset for start of plane 

    for(int n=0;n<e1;n++){
      for(int b=0;b<e2;b++){
	int ss= so+n*stride+b;

	Vector<iSinglet<vector_type> > phase(Nmom);

	for(int m=0;m<Nmom;m++) phase[m] = mom[m]._odata[ss];

	for(int i=0;i<Lblock;i++){

	  auto left = conjugate(lhs[i]._odata[ss]);
	  for(int j=0;j<Rblock;j++){

	    SpinMatrix_v vv;
	    auto right = rhs[j]._odata[ss];
	    for(int s1=0;s1<Ns;s1++){
	    for(int s2=0;s2<Ns;s2++){
	      vv()(s1,s2)() = left()(s1)(0) * right()(s2)(0)
		+             left()(s1)(1) * right()(s2)(1)
		+             left()(s1)(2) * right()(s2)(2);
	    }}
	    
	    // After getting the sitewise product do the mom phase loop
	    for ( int m=0;m<Nmom;m++){
	      int idx = m+Nmom*i+Nmom*Lblock*j+Nmom*Lblock*Rblock*r;
	      lvSum[idx]=lvSum[idx]+vv*phase[m];
	    }
	  
	  }
	}
      }
    }
  }

  // Sum across simd lanes in the plane, breaking out orthog dir.
  parallel_for(int rt=0;rt<rd;rt++){

    std::vector<int> icoor(Nd);
    std::vector<SpinMatrix_s> extracted(Nsimd);               


    for(int i=0;i<Lblock;i++){
    for(int j=0;j<Rblock;j++){
    for(int m=0;m<Nmom;m++){

      int ij_rdx = m+Nmom*i+Nmom*Lblock*j+Nmom*Lblock*Rblock*rt;

      extract(lvSum[ij_rdx],extracted);

      for(int idx=0;idx<Nsimd;idx++){

	grid->iCoorFromIindex(icoor,idx);

	int ldx    = rt+icoor[orthogdim]*rd;

	int ij_ldx = m+Nmom*i+Nmom*Lblock*j+Nmom*Lblock*Rblock*ldx;

	lsSum[ij_ldx]=lsSum[ij_ldx]+extracted[idx];

      }
    }}}
  }

  assert(mat.dimension(0) == Nt);
  assert(mat.dimension(1) == Nmom);
  assert(mat.dimension(2) == Ngamma);
  assert(mat.dimension(3) == Lblock);
  assert(mat.dimension(4) == Rblock);
  mat.setZero();
  parallel_for(int t=0;t<fd;t++)
  {
    int pt = t / ld; // processor plane
    int lt = t % ld;
    if (pt == grid->_processor_coor[orthogdim]){
      for(int i=0;i<Lblock;i++){
	for(int j=0;j<Rblock;j++){
	  for(int m=0;m<Nmom;m++){
	    int ij_dx = m+Nmom*i + Nmom*Lblock * j + Nmom*Lblock * Rblock * lt;
	    for(int mu=0;mu<Ngamma;mu++){
	      mat(t,m,mu,i,j) = trace(lsSum[ij_dx]*Gamma(gammas[mu]));
	    }
	  }
	}
      }
    }
  }
  grid->GlobalSumVector(&mat(0,0,0,0,0),Nmom*Rblock*Lblock*Nt*Ngamma);

  return;
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2AMesonField<FImpl>::execute(void)
{
    LOG(Message) << "Computing A2A meson field" << std::endl;

    auto &a2a = envGet(A2ABase, par().A2A + "_class");
    
    // 2+6+4+4 = 16 gammas
    // Ordering defined here
    std::vector<Gamma::Algebra> gammas ( {
	  Gamma::Algebra::Identity,    
          Gamma::Algebra::Gamma5,
	  Gamma::Algebra::GammaX,
	  Gamma::Algebra::GammaY,
	  Gamma::Algebra::GammaZ,
	  Gamma::Algebra::GammaT,
	  Gamma::Algebra::GammaXGamma5,
	  Gamma::Algebra::GammaYGamma5,
	  Gamma::Algebra::GammaZGamma5,
	  Gamma::Algebra::GammaTGamma5,
	  Gamma::Algebra::SigmaXY,
	  Gamma::Algebra::SigmaXZ,
	  Gamma::Algebra::SigmaXT,
	  Gamma::Algebra::SigmaYZ,
	  Gamma::Algebra::SigmaYT,
	  Gamma::Algebra::SigmaZT
    });

    ///////////////////////////////////////////////
    // Square assumption for now Nl = Nr = N
    ///////////////////////////////////////////////
    int nt = env().getDim(Tp);
    int N  = par().N;
    int Nl = par().Nl;
    int ngamma = gammas.size();

    ///////////////////////////////////////////////
    // Momentum setup
    ///////////////////////////////////////////////
    std::vector<LatticeComplex> phases(1,env().getGrid(1));
    int nmom = phases.size();
    phases[0] = Complex(1.0);

    Eigen::Tensor<ComplexD,5> mesonField       (nmom,ngamma,nt,N,N);    

    LOG(Message) << "N = Nh+Nl for A2A MesonField is " << N << std::endl;

    envGetTmp(std::vector<FermionField>, w);
    envGetTmp(std::vector<FermionField>, v);
    envGetTmp(FermionField, tmp_5d);

    LOG(Message) << "Finding v and w vectors for N =  " << N << std::endl;

    int schurBlock = par().schurBlock;
    int cacheBlock = par().cacheBlock;
    for(int i_base=0;i_base<N;i_base+=schurBlock){
    for(int j_base=0;j_base<N;j_base+=schurBlock){

      ///////////////////////////////////////////////////////////////
      // Get the W and V vectors for this schurBlock^2 set of terms
      ///////////////////////////////////////////////////////////////
      int i_max = MIN(N,i_base+schurBlock);
      int j_max = MIN(N,j_base+schurBlock);

      int N_i   = i_max-i_base;
      int N_j   = j_max-j_base;

      for(int ii =0;ii+i_base< i_max;ii++) a2a.return_v(i_base+ii, tmp_5d, v[ii]);
      for(int jj =0;jj+j_base< j_max;jj++) a2a.return_w(j_base+jj, tmp_5d, w[jj]);

      LOG(Message) << "Found v vectors " << i_base <<" .. " << i_max-1 << std::endl;
      LOG(Message) << "Found w vectors " << j_base <<" .. " << j_max-1 << std::endl;

      ///////////////////////////////////////////////////////////////
      // Do a cache blocked chunk of the contractions
      /////////////////////////////////////////////////////////////// 
      Eigen::Tensor<ComplexD,5> mesonFieldBlocked(nmom,ngamma,nt,N_i,N_j);    

      MesonField(mesonFieldBlocked, w, v, gammas, phases,Tp);

      ///////////////////////////////////////////////////////////////
      // Copy out to full meson field tensor
      /////////////////////////////////////////////////////////////// 
      for(int ii =0;ii< N_i;ii++) {
      for(int jj =0;jj< N_j;jj++) {
      for(int m =0;m< nmom;m++) {
      for(int g =0;g< ngamma;g++) {
      for(int t =0;t< nt;t++) {
	mesonField(m,g,t,i_base+ii,j_base+jj) = mesonFieldBlocked(m,g,t,ii,jj);
      }}}}}

      LOG(Message) << "Contracted MesonFields " <<std::endl;

    }}

    //    saveResult(par().output, "meson", result);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_A2AMesonField_hpp_
