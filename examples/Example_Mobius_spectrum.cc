/*
 * Warning: This code illustrative only: not well tested, and not meant for production use
 * without regression / tests being applied
 */

#include <Grid/Grid.h>

using namespace std;
using namespace Grid;

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
void Solve(Action &D,LatticePropagator &source,LatticePropagator &propagator)
{
  GridBase *UGrid = D.GaugeGrid();
  GridBase *FGrid = D.FermionGrid();

  LatticeFermion src4  (UGrid); 
  LatticeFermion src5  (FGrid); 
  LatticeFermion result5(FGrid);
  LatticeFermion result4(UGrid);
  
  ConjugateGradient<LatticeFermion> CG(1.0e-8,100000);
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

int main (int argc, char ** argv)
{
  const int Ls=8;

  Grid_init(&argc,&argv);

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
  std::vector<int> seeds4({1,2,3,4}); 
  GridParallelRNG          RNG4(UGrid);  RNG4.SeedFixedIntegers(seeds4);

  LatticeGaugeField Umu(UGrid);
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
    //    SU<Nc>::HotConfiguration(RNG4,Umu);
    config="HotConfig";
  }

  std::vector<RealD> masses({ 0.03,0.04,0.45} ); // u/d, s, c ??

  int nmass = masses.size();

  std::vector<MobiusFermionR *> FermActs;
  
  std::cout<<GridLogMessage <<"======================"<<std::endl;
  std::cout<<GridLogMessage <<"MobiusFermion action as Scaled Shamir kernel"<<std::endl;
  std::cout<<GridLogMessage <<"======================"<<std::endl;

  for(auto mass: masses) {

    RealD M5=1.0;
    RealD b=1.5;// Scale factor b+c=2, b-c=1
    RealD c=0.5;
    
    FermActs.push_back(new MobiusFermionR(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,b,c));
   
  }

  LatticePropagator point_source(UGrid);
  LatticePropagator wall_source(UGrid);
  LatticePropagator gaussian_source(UGrid);

  Coordinate Origin({0,0,0,0});
  PointSource   (Origin,point_source);
  Z2WallSource  (RNG4,0,wall_source);
  GaussianSource(Origin,Umu,gaussian_source);
  
  std::vector<LatticePropagator> PointProps(nmass,UGrid);
  std::vector<LatticePropagator> GaussProps(nmass,UGrid);
  std::vector<LatticePropagator> Z2Props   (nmass,UGrid);

  for(int m=0;m<nmass;m++) {
    
    Solve(*FermActs[m],point_source   ,PointProps[m]);
    Solve(*FermActs[m],gaussian_source,GaussProps[m]);
    Solve(*FermActs[m],wall_source    ,Z2Props[m]);
  
  }

  LatticeComplex phase(UGrid);
  Coordinate mom({0,0,0,0});
  MakePhase(mom,phase);
  
  for(int m1=0 ;m1<nmass;m1++) {
  for(int m2=m1;m2<nmass;m2++) {
    std::stringstream ssp,ssg,ssz;

    ssp<<config<< "_m" << m1 << "_m"<< m2 << "_point_meson.xml";
    ssg<<config<< "_m" << m1 << "_m"<< m2 << "_smeared_meson.xml";
    ssz<<config<< "_m" << m1 << "_m"<< m2 << "_wall_meson.xml";

    MesonTrace(ssp.str(),PointProps[m1],PointProps[m2],phase);
    MesonTrace(ssg.str(),GaussProps[m1],GaussProps[m2],phase);
    MesonTrace(ssz.str(),Z2Props[m1],Z2Props[m2],phase);
  }}

  Grid_finalize();
}



