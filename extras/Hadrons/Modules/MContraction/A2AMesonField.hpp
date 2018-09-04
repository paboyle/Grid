#ifndef Hadrons_MContraction_A2AMesonField_hpp_
#define Hadrons_MContraction_A2AMesonField_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>
#include <Grid/Hadrons/AllToAllVectors.hpp>

#include <Grid/Hadrons/Modules/MContraction/A2Autils.hpp>

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
	A2Autils<FImpl>::MesonField(mesonFieldBlocked, 
				    &w[ii], 
				    &v[jj], gammas, phases,Tp);
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
