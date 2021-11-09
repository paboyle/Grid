/*
 * Warning: This code illustrative only: not well tested, and not meant for production use
 * without regression / tests being applied
 */

#include <Grid/Grid.h>

using namespace std;
using namespace Grid;
typedef SpinColourMatrix Propagator;
typedef SpinColourVector Fermion;
typedef PeriodicGimplR   GimplR;

template<class Gimpl,class Field> class CovariantLaplacianCshift : public SparseMatrixBase<Field>
{
public:
  INHERIT_GIMPL_TYPES(Gimpl);

  GridBase *grid;
  GaugeField U;
  
  CovariantLaplacianCshift(GaugeField &_U)    :
    grid(_U.Grid()),
    U(_U) {  };

  virtual GridBase *Grid(void) { return grid; };

  virtual void  M    (const Field &in, Field &out)
  {
    out=Zero();
    for(int mu=0;mu<Nd-1;mu++) {
      GaugeLinkField Umu = PeekIndex<LorentzIndex>(U, mu); // NB: Inefficent
      out = out - Gimpl::CovShiftForward(Umu,mu,in);    
      out = out - Gimpl::CovShiftBackward(Umu,mu,in);    
      out = out + 2.0*in;
    }
  };
  virtual void  Mdag (const Field &in, Field &out) { M(in,out);}; // Laplacian is hermitian
  virtual  void Mdiag    (const Field &in, Field &out)                  {assert(0);}; // Unimplemented need only for multigrid
  virtual  void Mdir     (const Field &in, Field &out,int dir, int disp){assert(0);}; // Unimplemented need only for multigrid
  virtual  void MdirAll  (const Field &in, std::vector<Field> &out)     {assert(0);}; // Unimplemented need only for multigrid
};

void MakePhase(Coordinate mom,LatticeComplex &phase)
{
  GridBase *grid = phase.Grid();
  auto latt_size = grid->GlobalDimensions();
  ComplexD ci(0.0,1.0);
  phase=Zero();

  LatticeComplex coor(phase.Grid());
  for(int mu=0;mu<Nd;mu++){
    RealD TwoPiL =  M_PI * 2.0/ latt_size[mu];
    LatticeCoordinate(coor,mu);
    phase = phase + (TwoPiL * mom[mu]) * coor;
  }
  phase = exp(phase*ci);
}
void LinkSmear(int nstep, RealD rho,LatticeGaugeField &Uin,LatticeGaugeField &Usmr)
{
  Smear_Stout<GimplR> Stout(rho);
  LatticeGaugeField Utmp(Uin.Grid());
  Utmp = Uin;
  for(int i=0;i<nstep;i++){
    Stout.smear(Usmr,Utmp);
    Utmp = Usmr;
  }
}
void PointSource(Coordinate &coor,LatticePropagator &source)
{
  //  Coordinate coor({0,0,0,0});
  source=Zero();
  SpinColourMatrix kronecker; kronecker=1.0;
  pokeSite(kronecker,source,coor);
}
void GFWallSource(int tslice,LatticePropagator &source)
{
  GridBase *grid = source.Grid();
  LatticeComplex one(grid); one = ComplexD(1.0,0.0);
  LatticeComplex zz(grid); zz=Zero();
  LatticeInteger t(grid);
  LatticeCoordinate(t,Tdir);
  one = where(t==Integer(tslice), one, zz);
  source = 1.0;
  source = source * one;
}

void Z2WallSource(GridParallelRNG &RNG,int tslice,LatticePropagator &source)
{
  GridBase *grid = source.Grid();
  LatticeComplex noise(grid);
  LatticeComplex zz(grid); zz=Zero();
  LatticeInteger t(grid);

  RealD nrm=1.0/sqrt(2);
  bernoulli(RNG, noise); // 0,1 50:50

  noise = (2.*noise - Complex(1,1))*nrm;

  LatticeCoordinate(t,Tdir);
  noise = where(t==Integer(tslice), noise, zz);

  source = 1.0;
  source = source*noise;
  std::cout << " Z2 wall " << norm2(source) << std::endl;
}
void GaugeFix(LatticeGaugeField &U,LatticeGaugeField &Ufix)
{
  Real alpha=0.05;

  Real plaq=WilsonLoops<GimplR>::avgPlaquette(U);

  std::cout << " Initial plaquette "<<plaq << std::endl;

  LatticeColourMatrix   xform(U.Grid()); 
  Ufix = U;
  int orthog=Nd-1;
  FourierAcceleratedGaugeFixer<GimplR>::SteepestDescentGaugeFix(Ufix,xform,alpha,100000,1.0e-14, 1.0e-14,true,orthog);
  
  plaq=WilsonLoops<GimplR>::avgPlaquette(Ufix);

  std::cout << " Final plaquette "<<plaq << std::endl;
}
template<class Field>
void GaussianSmear(LatticeGaugeField &U,Field &unsmeared,Field &smeared)
{
  typedef CovariantLaplacianCshift <GimplR,Field> Laplacian_t;
  Laplacian_t Laplacian(U);

  Integer Iterations = 40;
  Real width = 2.0;
  Real coeff = (width*width) / Real(4*Iterations);

  Field tmp(U.Grid());
  smeared=unsmeared;
  //  chi = (1-p^2/2N)^N kronecker
  for(int n = 0; n < Iterations; ++n) {
    Laplacian.M(smeared,tmp);
    smeared = smeared - coeff*tmp;
    std::cout << " smear iter " << n<<" " <<norm2(smeared)<<std::endl;
  }
}
void GaussianSource(Coordinate &site,LatticeGaugeField &U,LatticePropagator &source)
{
  LatticePropagator tmp(source.Grid());
  PointSource(site,source);
  std::cout << " GaussianSource Kronecker "<< norm2(source)<<std::endl;
  tmp = source;
  GaussianSmear(U,tmp,source);
  std::cout << " GaussianSource Smeared "<< norm2(source)<<std::endl;
}
void GaussianWallSource(GridParallelRNG &RNG,int tslice,LatticeGaugeField &U,LatticePropagator &source)
{
  Z2WallSource(RNG,tslice,source);
  auto tmp = source;
  GaussianSmear(U,tmp,source);
}
void SequentialSource(int tslice,Coordinate &mom,LatticePropagator &spectator,LatticePropagator &source)
{
  assert(mom.size()==Nd);
  assert(mom[Tdir] == 0);

  GridBase * grid = spectator.Grid();

  LatticeInteger ts(grid);
  LatticeCoordinate(ts,Tdir);
  source = Zero();
  source = where(ts==Integer(tslice),spectator,source); // Stick in a slice of the spectator, zero everywhere else

  LatticeComplex phase(grid);
  MakePhase(mom,phase);

  source = source *phase;
}
template<class Action>
void Solve(Action &D,LatticePropagator &source,LatticePropagator &propagator)
{
  GridBase *UGrid = D.GaugeGrid();
  GridBase *FGrid = D.FermionGrid();

  LatticeFermion src4  (UGrid); 
  LatticeFermion src5  (FGrid); 
  LatticeFermion result5(FGrid);
  LatticeFermion result4(UGrid);
  
  ConjugateGradient<LatticeFermion> CG(1.0e-12,100000);
  SchurRedBlackDiagTwoSolve<LatticeFermion> schur(CG);
  ZeroGuesser<LatticeFermion> ZG; // Could be a DeflatedGuesser if have eigenvectors
  for(int s=0;s<Nd;s++){
    for(int c=0;c<Nc;c++){
      PropToFerm<Action>(src4,source,s,c);

      D.ImportPhysicalFermionSource(src4,src5);

      result5=Zero();
      schur(D,src5,result5,ZG);
      std::cout<<GridLogMessage
	       <<"spin "<<s<<" color "<<c
	       <<" norm2(src5d) "   <<norm2(src5)
               <<" norm2(result5d) "<<norm2(result5)<<std::endl;

      D.ExportPhysicalFermionSolution(result5,result4);

      FermToProp<Action>(propagator,result4,s,c);
    }
  }
}

class MesonFile: Serializable {
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(MesonFile, std::vector<std::vector<Complex> >, data);
};

void MesonTrace(std::string file,LatticePropagator &q1,LatticePropagator &q2,LatticeComplex &phase)
{
  const int nchannel=4;
  Gamma::Algebra Gammas[nchannel][2] = {
    {Gamma::Algebra::Gamma5      ,Gamma::Algebra::Gamma5},
    {Gamma::Algebra::GammaTGamma5,Gamma::Algebra::GammaTGamma5},
    {Gamma::Algebra::GammaTGamma5,Gamma::Algebra::Gamma5},
    {Gamma::Algebra::Gamma5      ,Gamma::Algebra::GammaTGamma5}
  };

  Gamma G5(Gamma::Algebra::Gamma5);

  LatticeComplex meson_CF(q1.Grid());
  MesonFile MF;

  for(int ch=0;ch<nchannel;ch++){

    Gamma Gsrc(Gammas[ch][0]);
    Gamma Gsnk(Gammas[ch][1]);

    meson_CF = trace(G5*adj(q1)*G5*Gsnk*q2*adj(Gsrc));

    std::vector<TComplex> meson_T;
    sliceSum(meson_CF,meson_T, Tdir);

    int nt=meson_T.size();

    std::vector<Complex> corr(nt);
    for(int t=0;t<nt;t++){
      corr[t] = TensorRemove(meson_T[t]); // Yes this is ugly, not figured a work around
      std::cout << " channel "<<ch<<" t "<<t<<" " <<corr[t]<<std::endl;
    }
    MF.data.push_back(corr);
  }

  {
    XmlWriter WR(file);
    write(WR,"MesonFile",MF);
  }
}

void Meson3pt(std::string file,LatticePropagator &q1,LatticePropagator &q2,LatticeComplex &phase)
{
  const int nchannel=4;
  Gamma::Algebra Gammas[nchannel][2] = {
    {Gamma::Algebra::Gamma5      ,Gamma::Algebra::GammaX},
    {Gamma::Algebra::Gamma5      ,Gamma::Algebra::GammaY},
    {Gamma::Algebra::Gamma5      ,Gamma::Algebra::GammaZ},
    {Gamma::Algebra::Gamma5      ,Gamma::Algebra::GammaT}
  };

  Gamma G5(Gamma::Algebra::Gamma5);

  LatticeComplex meson_CF(q1.Grid());
  MesonFile MF;

  for(int ch=0;ch<nchannel;ch++){

    Gamma Gsrc(Gammas[ch][0]);
    Gamma Gsnk(Gammas[ch][1]);

    meson_CF = trace(G5*adj(q1)*G5*Gsnk*q2*adj(Gsrc));

    std::vector<TComplex> meson_T;
    sliceSum(meson_CF,meson_T, Tdir);

    int nt=meson_T.size();

    std::vector<Complex> corr(nt);
    for(int t=0;t<nt;t++){
      corr[t] = TensorRemove(meson_T[t]); // Yes this is ugly, not figured a work around
      std::cout << " channel "<<ch<<" t "<<t<<" " <<corr[t]<<std::endl;
    }
    MF.data.push_back(corr);
  }

  {
    XmlWriter WR(file);
    write(WR,"MesonFile",MF);
  }
}


void WallSinkMesonTrace(std::string file,std::vector<Propagator> &q1,std::vector<Propagator> &q2)
{
  const int nchannel=4;
  Gamma::Algebra Gammas[nchannel][2] = {
    {Gamma::Algebra::Gamma5      ,Gamma::Algebra::Gamma5},
    {Gamma::Algebra::GammaTGamma5,Gamma::Algebra::GammaTGamma5},
    {Gamma::Algebra::GammaTGamma5,Gamma::Algebra::Gamma5},
    {Gamma::Algebra::Gamma5      ,Gamma::Algebra::GammaTGamma5}
  };

  Gamma G5(Gamma::Algebra::Gamma5);
  int nt=q1.size();
  std::vector<Complex> meson_CF(nt);
  MesonFile MF;

  for(int ch=0;ch<nchannel;ch++){

    Gamma Gsrc(Gammas[ch][0]);
    Gamma Gsnk(Gammas[ch][1]);

    std::vector<Complex> corr(nt);
    for(int t=0;t<nt;t++){
      meson_CF[t] = trace(G5*adj(q1[t])*G5*Gsnk*q2[t]*adj(Gsrc));
      corr[t] = TensorRemove(meson_CF[t]); // Yes this is ugly, not figured a work around
      std::cout << " channel "<<ch<<" t "<<t<<" " <<corr[t]<<std::endl;
    }
    MF.data.push_back(corr);
  }

  {
    XmlWriter WR(file);
    write(WR,"MesonFile",MF);
  }
}
int make_idx(int p, int m,int nmom)
{
  if (m==0) return p;
  assert(p==0);
  return nmom + m - 1;
}

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  // Double precision grids
  auto latt = GridDefaultLatt();
  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), 
								   GridDefaultSimd(Nd,vComplex::Nsimd()),
								   GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);


  LatticeGaugeField Umu(UGrid);
  LatticeGaugeField Utmp(UGrid);
  LatticeGaugeField Usmr(UGrid);
  std::string config;
  if( argc > 1 && argv[1][0] != '-' )
  {
    std::cout<<GridLogMessage <<"Loading configuration from "<<argv[1]<<std::endl;
    FieldMetaData header;
    NerscIO::readConfiguration(Umu, header, argv[1]);
    config=argv[1];
  }
  else
  {
    std::cout<<GridLogMessage <<"Using hot configuration"<<std::endl;
    SU<Nc>::ColdConfiguration(Umu);
    config="ColdConfig";
  }
  //  GaugeFix(Umu,Utmp);
  //  Umu=Utmp;

  int nsmr=3;
  RealD rho=0.1;
  LinkSmear(nsmr,rho,Umu,Usmr);


  std::vector<int>   smeared_link({ 0,0,1} ); 
  std::vector<RealD> masses({ 0.004,0.02477,0.447} ); // u/d, s, c ??
  std::vector<RealD> M5s   ({ 1.8,1.8,1.0} ); 
  std::vector<RealD> bs   ({ 1.0,1.0,1.5} );  // DDM
  std::vector<RealD> cs   ({ 0.0,0.0,0.5} );  // DDM
  std::vector<int>   Ls_s ({ 16,16,12} );
  std::vector<GridCartesian *> FGrids;
  std::vector<GridRedBlackCartesian *> FrbGrids;

  std::vector<Coordinate> momenta;
  momenta.push_back(Coordinate({0,0,0,0}));
  momenta.push_back(Coordinate({1,0,0,0}));
  momenta.push_back(Coordinate({2,0,0,0}));

  int nmass = masses.size();
  int nmom  = momenta.size();

  std::vector<MobiusFermionR *> FermActs;
  
  std::cout<<GridLogMessage <<"======================"<<std::endl;
  std::cout<<GridLogMessage <<"MobiusFermion action as Scaled Shamir kernel"<<std::endl;
  std::cout<<GridLogMessage <<"======================"<<std::endl;

  std::vector<Complex> boundary = {1,1,1,-1};
  typedef MobiusFermionR FermionAction;
  FermionAction::ImplParams Params(boundary);

  for(int m=0;m<masses.size();m++) {

    RealD mass = masses[m];
    RealD M5   = M5s[m];
    RealD b    = bs[m];
    RealD c    = cs[m];
    int   Ls   = Ls_s[m];

    if ( smeared_link[m] ) Utmp = Usmr;
    else                   Utmp = Umu;
    
    FGrids.push_back(SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid));
    FrbGrids.push_back(SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid));

    FermActs.push_back(new MobiusFermionR(Utmp,*FGrids[m],*FrbGrids[m],*UGrid,*UrbGrid,mass,M5,b,c,Params));
  }

  LatticePropagator z2wall_source(UGrid);
  LatticePropagator gfwall_source(UGrid);
  LatticePropagator phased_prop(UGrid);

  int tslice = 0;
  int tseq=(tslice+16)%latt[Nd-1];
  //////////////////////////////////////////////////////////////////////
  // RNG seeded for Z2 wall
  //////////////////////////////////////////////////////////////////////
  // You can manage seeds however you like.
  // Recommend SeedUniqueString.
  //////////////////////////////////////////////////////////////////////
  GridParallelRNG          RNG4(UGrid);  RNG4.SeedUniqueString("Study2-Source_Z2_p_0_0_0_t_0-880");
  Z2WallSource  (RNG4,tslice,z2wall_source);
  GFWallSource  (tslice,gfwall_source);

  std::vector<LatticeComplex> phase(nmom,UGrid);
  for(int m=0;m<nmom;m++){
    MakePhase(momenta[m],phase[m]);
  }

  std::vector<LatticePropagator> Z2Props   (nmom+nmass-1,UGrid);
  std::vector<LatticePropagator> GFProps   (nmom+nmass-1,UGrid);
  for(int p=0;p<nmom;p++) {
    int m=0;
    int idx = make_idx(p,m,nmom);
    phased_prop = z2wall_source * phase[p];
    Solve(*FermActs[m],phased_prop  ,Z2Props[idx]);

    phased_prop = gfwall_source * phase[p];
    Solve(*FermActs[m],phased_prop  ,GFProps[idx]);
  }
  for(int m=1;m<nmass;m++) {
    int p=0;
    int idx = make_idx(p,m,nmom);
    phased_prop = z2wall_source;
    Solve(*FermActs[m],phased_prop  ,Z2Props[idx]);

    phased_prop = gfwall_source;
    Solve(*FermActs[m],phased_prop  ,GFProps[idx]);
  }

  std::vector<std::vector<Propagator> > wsnk_z2Props(nmom+nmass-1);
  std::vector<std::vector<Propagator> > wsnk_gfProps(nmom+nmass-1);

  // Non-zero kaon and point and D two point
  // WW stick momentum on m1 (lighter)
  //     zero momentum on m2
  for(int m1=0;m1<nmass;m1++) {
  for(int m2=m1;m2<nmass;m2++) {
    int pmax = (m1==0)? nmom:1;
    for(int p=0;p<pmax;p++){

      std::stringstream ssg,ssz;
      std::stringstream wssg,wssz;

      int idx1 = make_idx(p,m1,nmom);
      int idx2 = make_idx(0,m2,nmom);

      /// Point sinks
      ssg<<config<<"_p"<<p<< "_m" << m1 << "_m"<< m2 << "_p_gf_meson.xml";
      ssz<<config<<"_p"<<p<< "_m" << m1 << "_m"<< m2 << "_p_z2_meson.xml";
      MesonTrace(ssz.str(),Z2Props[idx1],Z2Props[idx2],phase[p]); // Q1 is conjugated
      MesonTrace(ssg.str(),GFProps[idx1],GFProps[idx2],phase[p]); 
      
      /// Wall sinks
      wssg<<config<<"_p"<<p<< "_m" << m1 << "_m"<< m2 << "_w_gf_meson.xml";
      wssz<<config<<"_p"<<p<< "_m" << m1 << "_m"<< m2 << "_w_z2_meson.xml";
      
      phased_prop = GFProps[m2] * phase[p];
      sliceSum(phased_prop,wsnk_gfProps[m1],Tdir);
      sliceSum(GFProps[m1],wsnk_gfProps[m2],Tdir);
      WallSinkMesonTrace(wssg.str(),wsnk_gfProps[m1],wsnk_gfProps[m2]);

      phased_prop = Z2Props[m2] * phase[p];
      sliceSum(phased_prop,wsnk_gfProps[m1],Tdir);
      sliceSum(Z2Props[m1],wsnk_gfProps[m2],Tdir);
      WallSinkMesonTrace(wssz.str(),wsnk_z2Props[m1],wsnk_z2Props[m2]);
    }
  }}


  /////////////////////////////////////
  // Sequential solves
  /////////////////////////////////////
  LatticePropagator  seq_wsnk_z2src(UGrid);
  LatticePropagator  seq_wsnk_gfsrc(UGrid);
  LatticePropagator  seq_psnk_z2src(UGrid);
  LatticePropagator  seq_psnk_gfsrc(UGrid);
  LatticePropagator source(UGrid);
  for(int m=0;m<nmass-1;m++){
    int spect_idx = make_idx(0,m,nmom);
    int charm=nmass-1;

    SequentialSource(tseq,momenta[0],GFProps[spect_idx],source);
    Solve(*FermActs[charm],source,seq_psnk_gfsrc);
    
    SequentialSource(tseq,momenta[0],Z2Props[spect_idx],source);
    Solve(*FermActs[charm],source,seq_psnk_z2src);

    // Todo need wall sequential solve
    for(int p=0;p<nmom;p++){
      int active_idx = make_idx(p,0,nmom);
      std::stringstream seq_3pt_p_z2;
      std::stringstream seq_3pt_p_gf;
      std::stringstream seq_3pt_w_z2;
      std::stringstream seq_3pt_w_gf;
      seq_3pt_p_z2  <<config<<"_3pt_p"<<p<< "_m" << m << "_p_z2_meson.xml";
      seq_3pt_p_gf  <<config<<"_3pt_p"<<p<< "_m" << m << "_p_gf_meson.xml";
      seq_3pt_w_z2  <<config<<"_3pt_p"<<p<< "_m" << m << "_w_z2_meson.xml";
      seq_3pt_w_gf  <<config<<"_3pt_p"<<p<< "_m" << m << "_w_gf_meson.xml";
      Meson3pt(seq_3pt_p_gf.str(),GFProps[active_idx],seq_psnk_gfsrc,phase[p]);
      Meson3pt(seq_3pt_p_z2.str(),Z2Props[active_idx],seq_psnk_z2src,phase[p]);
    }    
  }
  
  Grid_finalize();
}



