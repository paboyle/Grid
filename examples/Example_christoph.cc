/*
 * Warning: This code illustrative only: not well tested, and not meant for production use
 * without regression / tests being applied
 */

#include <Grid/Grid.h>

using namespace std;
using namespace Grid;

RealD LLscale =1.0;
RealD LCscale =1.0;

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

void PointSource(Coordinate &coor,LatticePropagator &source)
{
  //  Coordinate coor({0,0,0,0});
  source=Zero();
  SpinColourMatrix kronecker; kronecker=1.0;
  pokeSite(kronecker,source,coor);
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
template<class Field>
void GaussianSmear(LatticeGaugeField &U,Field &unsmeared,Field &smeared)
{
  typedef CovariantLaplacianCshift <PeriodicGimplR,Field> Laplacian_t;
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
void MasslessFreePropagator(Action &D,LatticePropagator &source,LatticePropagator &propagator)
{			   
 GridBase *UGrid = source.Grid();
  GridBase *FGrid = D.FermionGrid();
  bool fiveD = true; //calculate 5d free propagator
  RealD mass = D.Mass();
  LatticeFermion src4  (UGrid);
  LatticeFermion result4  (UGrid);
  LatticeFermion result5(FGrid);
  LatticeFermion src5(FGrid);
  LatticePropagator prop5(FGrid);
  for(int s=0;s<Nd;s++){
    for(int c=0;c<Nc;c++){
 
      PropToFerm<Action>(src4,source,s,c);

      D.ImportPhysicalFermionSource(src4,src5);
      D.FreePropagator(src5,result5,mass,true);
      std::cout<<GridLogMessage
               <<"Free 5D prop spin "<<s<<" color "<<c
               <<" norm2(src5d) "   <<norm2(src5)
               <<" norm2(result5d) "<<norm2(result5)<<std::endl;

      D.ExportPhysicalFermionSolution(result5,result4);

      FermToProp<Action>(prop5,result5,s,c);
      FermToProp<Action>(propagator,result4,s,c);
    }
  }

  LatticePropagator Vector_mu(UGrid);
  LatticeComplex    VV (UGrid);
  std::vector<TComplex> sumVV;
  Gamma::Algebra GammaV[3] = {
    Gamma::Algebra::GammaX,
    Gamma::Algebra::GammaY,
    Gamma::Algebra::GammaZ
  };
  for( int mu=0;mu<3;mu++ ) {
    Gamma gV(GammaV[mu]);
    D.ContractConservedCurrent(prop5,prop5,Vector_mu,source,Current::Vector,mu);
    VV       = trace(gV*Vector_mu);     // (local) Vector-Vector conserved current
    sliceSum(VV,sumVV,Tdir);
    int Nt = sumVV.size();
    for(int t=0;t<Nt;t++){
      RealD Ct = real(TensorRemove(sumVV[t]))*LCscale;
      RealD Cont=0;
      if(t) Cont=1.0/(2 * M_PI *M_PI * t*t*t);
      std::cout<<GridLogMessage <<"VVc["<<mu<<"]["<<t<<"] "<< Ct
               << " 2 pi^2 t^3 C(t) "<< Ct/Cont << " delta Ct "<< Ct-Cont <<std::endl;
    }
  }
}
template<class Action>
void MasslessFreePropagator1(Action &D,LatticePropagator &source,LatticePropagator &propagator)
{			   
  bool fiveD = false; //calculate 4d free propagator
  RealD mass = D.Mass();
  GridBase *UGrid = source.Grid();
  LatticeFermion src4  (UGrid); 
  LatticeFermion result4  (UGrid); 
  for(int s=0;s<Nd;s++){
    for(int c=0;c<Nc;c++){
      PropToFerm<Action>(src4,source,s,c);
      D.FreePropagator(src4,result4,mass,false);
      FermToProp<Action>(propagator,result4,s,c);
    }
  }
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
  LatticePropagator prop5(FGrid);
  
  ConjugateGradient<LatticeFermion> CG(1.0e-7,100000);
  SchurRedBlackDiagMooeeSolve<LatticeFermion> schur(CG);
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

      FermToProp<Action>(prop5,result5,s,c);
      FermToProp<Action>(propagator,result4,s,c);
    }
  }
  LatticePropagator Axial_mu(UGrid); 
  LatticePropagator Vector_mu(UGrid); 

  LatticeComplex    PA (UGrid); 
  LatticeComplex    VV (UGrid); 
  LatticeComplex    PJ5q(UGrid);
  LatticeComplex    PP (UGrid);

  std::vector<TComplex> sumPA;
  std::vector<TComplex> sumVV;
  std::vector<TComplex> sumPP;
  std::vector<TComplex> sumPJ5q;

  Gamma g5(Gamma::Algebra::Gamma5);
  D.ContractConservedCurrent(prop5,prop5,Axial_mu,source,Current::Axial,Tdir);
  PA       = trace(g5*Axial_mu);      // Pseudoscalar-Axial conserved current
  sliceSum(PA,sumPA,Tdir);

  int Nt{static_cast<int>(sumPA.size())};

  for(int t=0;t<Nt;t++) std::cout<<GridLogMessage <<"PAc["<<t<<"] "<<real(TensorRemove(sumPA[t]))*LCscale<<std::endl;

  PP       = trace(adj(propagator)*propagator); // Pseudoscalar density
  sliceSum(PP,sumPP,Tdir);
  for(int t=0;t<Nt;t++) std::cout<<GridLogMessage <<"PP["<<t<<"] "<<real(TensorRemove(sumPP[t]))*LCscale<<std::endl;
  
  D.ContractJ5q(prop5,PJ5q);
  sliceSum(PJ5q,sumPJ5q,Tdir);
  for(int t=0;t<Nt;t++) std::cout<<GridLogMessage <<"PJ5q["<<t<<"] "<<real(TensorRemove(sumPJ5q[t]))<<std::endl;

  Gamma::Algebra GammaV[3] = {
    Gamma::Algebra::GammaX,
    Gamma::Algebra::GammaY,
    Gamma::Algebra::GammaZ
  };
  for( int mu=0;mu<3;mu++ ) {
    Gamma gV(GammaV[mu]);
    D.ContractConservedCurrent(prop5,prop5,Vector_mu,source,Current::Vector,mu);
    //    auto ss=sliceSum(Vector_mu,Tdir);
    //    for(int t=0;t<Nt;t++) std::cout<<GridLogMessage <<"ss["<<mu<<"]["<<t<<"] "<<ss[t]<<std::endl;
    VV       = trace(gV*Vector_mu);     // (local) Vector-Vector conserved current
    sliceSum(VV,sumVV,Tdir);
    for(int t=0;t<Nt;t++){
      RealD Ct = real(TensorRemove(sumVV[t]))*LCscale;
      RealD Cont=0;
      if(t) Cont=1.0/(2 * M_PI *M_PI * t*t*t);
      std::cout<<GridLogMessage <<"VVc["<<mu<<"]["<<t<<"] "<< Ct
               << " 2 pi^2 t^3 C(t) "<< Ct/Cont << " delta Ct "<< Ct-Cont <<std::endl;
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
    {Gamma::Algebra::GammaXGamma5,Gamma::Algebra::GammaXGamma5},
    {Gamma::Algebra::GammaYGamma5,Gamma::Algebra::GammaYGamma5},
    {Gamma::Algebra::GammaZGamma5,Gamma::Algebra::GammaZGamma5},
    {Gamma::Algebra::Identity,Gamma::Algebra::Identity}
  };

  LatticeComplex meson_CF(q1.Grid());
  MesonFile MF;

  for(int ch=0;ch<nchannel;ch++){

    Gamma Gsrc(Gammas[ch][0]);
    Gamma Gsnk(Gammas[ch][1]);

    meson_CF = trace(adj(q1)*Gsnk*q2*adj(Gsrc));

    std::vector<TComplex> meson_T;
    sliceSum(meson_CF,meson_T, Tdir);

    int nt=meson_T.size();

    std::vector<Complex> corr(nt);
    for(int t=0;t<nt;t++){
      corr[t] = TensorRemove(meson_T[t])*LLscale; // Yes this is ugly, not figured a work around
      RealD Ct = real(corr[t]);
      RealD Cont=0;
      if(t) Cont=1.0/(2 * M_PI *M_PI * t*t*t);
      std::cout << " channel "<<ch<<" t "<<t<<" " <<real(corr[t])<< " 2 pi^2 t^3 C(t) "<< 2 * M_PI *M_PI * t*t*t * Ct
		<< " deltaC " <<Ct-Cont<<std::endl;
    }
    MF.data.push_back(corr);
  }

  {
    XmlWriter WR(file);
    write(WR,"MesonFile",MF);
  }
}

int main (int argc, char ** argv)
{

  Grid_init(&argc,&argv);
  int Ls= atoi(getenv("Ls"));

  // Double precision grids
  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), 
								   GridDefaultSimd(Nd,vComplex::Nsimd()),
								   GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  //////////////////////////////////////////////////////////////////////
  // You can manage seeds however you like.
  // Recommend SeedUniqueString.
  //////////////////////////////////////////////////////////////////////
  //  std::vector<int> seeds4({1,2,3,4}); 
  //  GridParallelRNG          RNG4(UGrid);  RNG4.SeedFixedIntegers(seeds4);

  LatticeGaugeField Umu(UGrid);
  std::string config;
  RealD M5=atof(getenv("M5"));
  RealD mq = atof(getenv("mass"));
  int   point_x = atoi(getenv("point_x"));
  int   point_y = atoi(getenv("point_y"));
  int   point_z = atoi(getenv("point_z"));
  int   point_t = atoi(getenv("point_t"));
  std::vector<RealD> masses({ mq} ); // u/d, s, c ??
  if( argc > 1 && argv[1][0] != '-' )
  {
    std::cout<<GridLogMessage <<"Loading configuration from "<<argv[1]<<std::endl;
    FieldMetaData header;
    NerscIO::readConfiguration(Umu, header, argv[1]);
    config=argv[1];
    LLscale = 1.0;
    LCscale = 1.0;
  } else {
    printf("Expected a configuration");
    exit(0);
  }

  int nmass = masses.size();

  typedef MobiusFermionD FermionActionD;
  std::vector<FermionActionD *> FermActs;
  
  std::cout<<GridLogMessage <<"======================"<<std::endl;
  std::cout<<GridLogMessage <<"DomainWallFermion action"<<std::endl;
  std::cout<<GridLogMessage <<"======================"<<std::endl;

  for(auto mass: masses) {
    std::vector<Complex> boundary = {1,1,1,-1};
    FermionActionD::ImplParams Params(boundary);
    RealD b=1.5;
    RealD c=0.5;
    FermActs.push_back(new FermionActionD(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,b,c));
  }

  LatticePropagator point_source(UGrid);

  Coordinate Origin({point_x,point_y,point_z,point_t});
  PointSource   (Origin,point_source);
  
  std::vector<LatticePropagator> PointProps(nmass,UGrid);

  for(int m=0;m<nmass;m++) {
    Solve(*FermActs[m],point_source   ,PointProps[m]);
  }

  LatticeComplex phase(UGrid);
  Coordinate mom({0,0,0,0});
  MakePhase(mom,phase);
  
  for(int m1=0 ;m1<nmass;m1++) {
  for(int m2=m1;m2<nmass;m2++) {
    std::stringstream ssp,ssg,ssz;

    ssp<<config<< "_m" << m1 << "_m"<< m2 << "_point_meson.xml";
    ssz<<config<< "_m" << m1 << "_m"<< m2 << "_free_meson.xml";

    std::cout << "CG determined VV correlation function"<<std::endl;
    MesonTrace(ssp.str(),PointProps[m1],PointProps[m2],phase);
    
  }}

  Grid_finalize();
}



