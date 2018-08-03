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
				    int, Nmom,
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
			    const LatticeFermion *lhs,
			    const LatticeFermion *rhs,
			    std::vector<Gamma::Algebra> gammas,
			    const std::vector<LatticeComplex > &mom,
			    int orthogdim,
			    double &t0,
			    double &t1,
			    double &t2,
			    double &t3);      
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
					 const LatticeFermion *lhs_wi,
					 const LatticeFermion *rhs_vj,
					 std::vector<Gamma::Algebra> gammas,
					 const std::vector<LatticeComplex > &mom,
					 int orthogdim,
					 double &t0,
					 double &t1,
					 double &t2,
					 double &t3) 
{
  typedef typename FImpl::SiteSpinor vobj;

  typedef typename vobj::scalar_object sobj;
  typedef typename vobj::scalar_type scalar_type;
  typedef typename vobj::vector_type vector_type;

  typedef iSpinMatrix<vector_type> SpinMatrix_v;
  typedef iSpinMatrix<scalar_type> SpinMatrix_s;
  
  int Lblock = mat.dimension(3); 
  int Rblock = mat.dimension(4);

  GridBase *grid = lhs_wi[0]._grid;
  
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

  t0-=usecond();
  // Nested parallelism would be ok
  // Wasting cores here. Test case r
  parallel_for(int r=0;r<rd;r++){

    int so=r*grid->_ostride[orthogdim]; // base offset for start of plane 

    for(int n=0;n<e1;n++){
      for(int b=0;b<e2;b++){

	int ss= so+n*stride+b;

	for(int i=0;i<Lblock;i++){

	  auto left = conjugate(lhs_wi[i]._odata[ss]);

	  for(int j=0;j<Rblock;j++){

	    SpinMatrix_v vv;
	    auto right = rhs_vj[j]._odata[ss];
	    for(int s1=0;s1<Ns;s1++){
	    for(int s2=0;s2<Ns;s2++){
	      vv()(s1,s2)() = left()(s2)(0) * right()(s1)(0)
		+             left()(s2)(1) * right()(s1)(1)
		+             left()(s2)(2) * right()(s1)(2);
	    }}
	    
	    // After getting the sitewise product do the mom phase loop
	    int base = Nmom*i+Nmom*Lblock*j+Nmom*Lblock*Rblock*r;
	    for ( int m=0;m<Nmom;m++){
	      int idx = m+base;
	      auto phase = mom[m]._odata[ss];
	      mac(&lvSum[idx],&vv,&phase);
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
  t1+=usecond();

  assert(mat.dimension(0) == Nmom);
  assert(mat.dimension(1) == Ngamma);
  assert(mat.dimension(2) == Nt);
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
	      for(int mu=0;mu<Ngamma;mu++){
		// this is a bit slow
		mat(m,mu,t,i,j) = trace(lsSum[ij_dx]*Gamma(gammas[mu]));
	      }
	    }
	  }
	}
      } else { 
	const scalar_type zz(0.0);
	for(int i=0;i<Lblock;i++){
	  for(int j=0;j<Rblock;j++){
	    for(int mu=0;mu<Ngamma;mu++){
	      for(int m=0;m<Nmom;m++){
		mat(m,mu,t,i,j) =zz;
	      }
	    }
	  }
	}
      }
    }
  }
  t2+=usecond();
  ////////////////////////////////////////////////////////////////////
  // This global sum is taking as much as 50% of time on 16 nodes
  // Vector size is 7 x 16 x 32 x 16 x 16 x sizeof(complex) = 2MB - 60MB depending on volume
  // Healthy size that should suffice
  ////////////////////////////////////////////////////////////////////
  t3-=usecond();
  grid->GlobalSumVector(&mat(0,0,0,0,0),Nmom*Ngamma*Nt*Lblock*Rblock);
  t3+=usecond();
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
          Gamma::Algebra::Gamma5,
	  Gamma::Algebra::Identity,    
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
    int nx = env().getDim(Xp);
    int ny = env().getDim(Yp);
    int nz = env().getDim(Zp);
    int N  = par().N;
    int Nl = par().Nl;
    int ngamma = gammas.size();

    int schurBlock = par().schurBlock;
    int cacheBlock = par().cacheBlock;
    int nmom       = par().Nmom;

    ///////////////////////////////////////////////
    // Momentum setup
    ///////////////////////////////////////////////
    GridBase *grid = env().getGrid(1);
    std::vector<LatticeComplex> phases(nmom,grid);
    for(int m=0;m<nmom;m++){
      phases[m] = Complex(1.0);    // All zero momentum for now
    }

    Eigen::Tensor<ComplexD,5> mesonField       (nmom,ngamma,nt,N,N);    
    LOG(Message) << "N = Nh+Nl for A2A MesonField is " << N << std::endl;

    envGetTmp(std::vector<FermionField>, w);
    envGetTmp(std::vector<FermionField>, v);
    envGetTmp(FermionField, tmp_5d);

    LOG(Message) << "Finding v and w vectors for N =  " << N << std::endl;

    //////////////////////////////////////////////////////////////////////////
    // i,j   is first  loop over SchurBlock factors reusing 5D matrices
    // ii,jj is second loop over cacheBlock factors for high perf contractoin
    // iii,jjj are loops within cacheBlock
    // Total index is sum of these  i+ii+iii etc...
    //////////////////////////////////////////////////////////////////////////
    
    double flops = 0.0;
    double bytes = 0.0;
    double vol   = nx*ny*nz*nt;
    double t_schur=0;
    double t_contr=0;
    double t_int_0=0;
    double t_int_1=0;
    double t_int_2=0;
    double t_int_3=0;

    double t0 = usecond();
    int N_i = N;
    int N_j = N;
    for(int i=0;i<N_i;i+=schurBlock){ //loop over SchurBlocking to suppress 5D matrix overhead
    for(int j=0;j<N_j;j+=schurBlock){
      

      ///////////////////////////////////////////////////////////////
      // Get the W and V vectors for this schurBlock^2 set of terms
      ///////////////////////////////////////////////////////////////
      int N_ii = MIN(N_i-i,schurBlock);
      int N_jj = MIN(N_j-j,schurBlock);

      t_schur-=usecond();
      for(int ii =0;ii < N_ii;ii++) a2a.return_w(i+ii, tmp_5d, w[ii]);
      for(int jj =0;jj < N_jj;jj++) a2a.return_v(j+jj, tmp_5d, v[jj]);
      t_schur+=usecond();

      LOG(Message) << "Found w vectors " << i <<" .. " << i+N_ii-1 << std::endl;
      LOG(Message) << "Found v vectors " << j <<" .. " << j+N_jj-1 << std::endl;

      ///////////////////////////////////////////////////////////////
      // Series of cache blocked chunks of the contractions within this SchurBlock
      /////////////////////////////////////////////////////////////// 
      for(int ii=0;ii<N_ii;ii+=cacheBlock){
      for(int jj=0;jj<N_jj;jj+=cacheBlock){

	int N_iii = MIN(N_ii-ii,cacheBlock);
	int N_jjj = MIN(N_jj-jj,cacheBlock);

	Eigen::Tensor<ComplexD,5> mesonFieldBlocked(nmom,ngamma,nt,N_iii,N_jjj);    

	t_contr-=usecond();
	MesonField(mesonFieldBlocked, &w[ii], &v[jj], gammas, phases,Tp,
		   t_int_0,t_int_1,t_int_2,t_int_3);
	t_contr+=usecond();
	flops += vol * ( 2 * 8.0 + 6.0 + 8.0*nmom) * N_iii*N_jjj*ngamma;

	bytes  += vol * (12.0 * sizeof(Complex) ) * N_iii*N_jjj
               +  vol * ( 2.0 * sizeof(Complex) *nmom ) * N_iii*N_jjj* ngamma;

	///////////////////////////////////////////////////////////////
	// Copy back to full meson field tensor
	/////////////////////////////////////////////////////////////// 
	parallel_for_nest2(int iii=0;iii< N_iii;iii++) {
        for(int jjj=0;jjj< N_jjj;jjj++) {
	  for(int m =0;m< nmom;m++) {
	  for(int g =0;g< ngamma;g++) {
          for(int t =0;t< nt;t++) {
	    mesonField(m,g,t,i+ii+iii,j+jj+jjj) = mesonFieldBlocked(m,g,t,iii,jjj);
	  }}}

	}}
      }}
    }}


    double nodes=grid->NodeCount();
    double t1 = usecond();
    LOG(Message) << " Contraction of MesonFields took "<<(t1-t0)/1.0e6<< " seconds "  << std::endl;
    LOG(Message) << " Schur "<<(t_schur)/1.0e6<< " seconds "  << std::endl;
    LOG(Message) << " Contr "<<(t_contr)/1.0e6<< " seconds "  << std::endl;
    LOG(Message) << " Intern0 "<<(t_int_0)/1.0e6<< " seconds "  << std::endl;
    LOG(Message) << " Intern1 "<<(t_int_1)/1.0e6<< " seconds "  << std::endl;
    LOG(Message) << " Intern2 "<<(t_int_2)/1.0e6<< " seconds "  << std::endl;
    LOG(Message) << " Intern3 "<<(t_int_3)/1.0e6<< " seconds "  << std::endl;

    double t_kernel = t_int_0 + t_int_1;
    LOG(Message) << " Arith "<<flops/(t_kernel)/1.0e3/nodes<< " Gflop/s / node "  << std::endl;
    LOG(Message) << " Arith "<<bytes/(t_kernel)/1.0e3/nodes<< " GB/s /node "  << std::endl;

    /////////////////////////////////////////////////////////////////////////
    // Test: Build the pion correlator (two end)
    // < PI_ij(t0) PI_ji (t0+t) >
    /////////////////////////////////////////////////////////////////////////
    std::vector<ComplexD> corr(nt,ComplexD(0.0));

    for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
      int m=0; // first momentum
      int g=0; // first gamma in above ordering is gamma5 for pion
      for(int t0=0;t0<nt;t0++){
      for(int t=0;t<nt;t++){
	int tt = (t0+t)%nt;
	corr[t] += mesonField(m,g,t0,i,j)* mesonField(m,g,tt,j,i);
      }}
    }}    
    for(int t=0;t<nt;t++) corr[t] = corr[t]/ (double)nt;

    for(int t=0;t<nt;t++) LOG(Message) << " " << t << " " << corr[t]<<std::endl;

    //    saveResult(par().output, "meson", result);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_A2AMesonField_hpp_
