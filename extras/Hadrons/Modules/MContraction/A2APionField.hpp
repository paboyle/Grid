#ifndef Hadrons_MContraction_A2APionField_hpp_
#define Hadrons_MContraction_A2APionField_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>
#include <Grid/Hadrons/AllToAllVectors.hpp>
#include <Grid/Hadrons/Modules/MContraction/A2Autils.hpp>

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

  typedef typename FImpl::SiteSpinor vobj;

  typedef typename vobj::scalar_object sobj;
  typedef typename vobj::scalar_type scalar_type;
  typedef typename vobj::vector_type vector_type;
  
  typedef iSpinMatrix<vector_type> SpinMatrix_v;
  typedef iSpinMatrix<scalar_type> SpinMatrix_s;
  typedef iSinglet<vector_type> Scalar_v;
  typedef iSinglet<scalar_type> Scalar_s;

  typedef A2AModesSchurDiagTwo<typename FImpl::FermionField, FMat, Solver> A2ABase;

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
	A2Autils<FImpl>::PionFieldWVmom(pionFieldWVmomB_ij, &wi[ii], &vj[jj], phases,Tp);
	A2Autils<FImpl>::PionFieldWVmom(pionFieldWVmomB_ji, &wj[jj], &vi[ii], phases,Tp);
	t_contr_vwm+=usecond();

	t_contr_vw-=usecond();
	A2Autils<FImpl>::PionFieldWV(pionFieldWVB_ij, &wi[ii], &vj[jj],Tp);
	A2Autils<FImpl>::PionFieldWV(pionFieldWVB_ji, &wj[jj], &vi[ii],Tp);
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

    A2Autils<FImpl>::PionFieldWW(pionFieldWW_ij, &wi[0], &wj[0],Tp);
    A2Autils<FImpl>::PionFieldVV(pionFieldVV_ji, &vj[0], &vi[0],Tp);
    A2Autils<FImpl>::PionFieldWW(pionFieldWW_ji, &wj[0], &wi[0],Tp);
    A2Autils<FImpl>::PionFieldVV(pionFieldVV_ij, &vi[0], &vj[0],Tp);


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


    A2Autils<FImpl>::PionFieldWW(pionFieldWW, &wi[0], &wi[0],Tp);
    A2Autils<FImpl>::PionFieldVV(pionFieldVV, &vi[0], &vi[0],Tp);
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

    A2Autils<FImpl>::DeltaFeq2  (dT,dT,DeltaF2_fig8,DeltaF2_trtr,
				 denom0,denom1,
				 pionFieldWW_ij,&vi[0],&vj[0],Tp);
    
    { 
      int g=0; // O_{VV+AA}
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
      
      /*
      if(g<5){
	for(int t=0;t<C1.size();t++){
	  LOG(Message) << " Wick1["<<g<<","<<t<< "] "<< C1[t]<<std::endl; 
	}
	for(int t=0;t<C2.size();t++){
	  LOG(Message) << " Wick2["<<g<<","<<t<< "] "<< C2[t]<<std::endl; 
	}
      }
      if( (g==9) || (g==7) ){ // P and At in above ordering
	for(int t=0;t<C3.size();t++){
	  LOG(Message) << " <G|P>["<<g<<","<<t<< "] "<< C3[t]<<std::endl; 
	}
      } 
      */
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_A2APionField_hpp_
  
