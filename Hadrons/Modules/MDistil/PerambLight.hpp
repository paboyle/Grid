#ifndef Hadrons_MDistil_perambulator_l_hpp_
#define Hadrons_MDistil_perambulator_l_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         perambulator_l                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MDistil)

class perambulator_lPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(perambulator_lPar,
                                    unsigned int, i);
};

template <typename FImpl>
class Tperambulator_l: public Module<perambulator_lPar>
{
public:
    // constructor
    Tperambulator_l(const std::string name);
    // destructor
    virtual ~Tperambulator_l(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(perambulator_l, Tperambulator_l<FIMPL>, MDistil);

/******************************************************************************
 *                 Tperambulator_l implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
Tperambulator_l<FImpl>::Tperambulator_l(const std::string name)
: Module<perambulator_lPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> Tperambulator_l<FImpl>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename FImpl>
std::vector<std::string> Tperambulator_l<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void Tperambulator_l<FImpl>::setup(void)
{
/*
    std::cout << "Compute perambulator from timeslice " << tsrc << std::endl;

    LatticeSpinColourVector dist_source(grid4d);
    LatticeSpinColourVector tmp2(grid4d);
    LatticeSpinColourVector tmp3d(grid3d);
    LatticeColourVector tmp3d_nospin(grid3d);
    LatticeColourVector evec3d(grid3d);
    LatticeColourVector tmp_nospin(grid4d);

    LatticeColourVector result_tmp(grid3d);

    LatticeSpinVector peramb_tmp(grid4d);
    LatticeFermion result(grid4d); result=zero; //Fermion = SpinColourVector!!!
    LatticeFermion result_single_component(grid4d); result_single_component=zero; //Fermion = SpinColourVector!!!
    LatticeColourVector result_nospin(grid4d); result_nospin=zero; //Fermion = SpinColourVector!!!
    LatticeColourVector result_3d(grid3d); result_3d=zero; //Fermion = SpinColourVector!!!
    LatticeFermion result_test(grid3d); result_test=zero; //Fermion = SpinColourVector!!!


    Real mass=SPar.mass;    // TODO Infile
    Real M5  =SPar.M5;     // TODO Infile
    std::cout << "init RBG "  << std::endl;
    GridRedBlackCartesian RBGrid(grid4d);
    std::cout << "init RBG done"  << std::endl;

    GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(DPar.Ls,grid4d);
    GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(DPar.Ls,grid4d);

    typedef DomainWallFermionR FermionAction;

    FermionAction Dop(Umu,*FGrid,*FrbGrid,*grid4d,RBGrid,mass,M5);

    MdagMLinearOperator<FermionAction,LatticeFermion> HermOp(Dop);
    ConjugateGradient<LatticeFermion> CG(SPar.CGPrecision,SPar.MaxIterations);
    SchurRedBlackDiagMooeeSolve<LatticeFermion> SchurSolver(CG);
*/
    
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void Tperambulator_l<FImpl>::execute(void)
{
/*
    for (int inoise = 0; inoise < nnoise; inoise++) {
      for (int dk = 0; dk < LI; dk++) {
        for (int dt = 0; dt < Nt_inv; dt++) {
          if(full_tdil) dt=tsrc; //this works for now, as longs as tsrc=0, but will crash otherwise!
          for (int ds = 0; ds < Ns; ds++) {
            std::cout <<  "LapH source vector from noise " << inoise << " and dilution component (d_k,d_t,d_alpha) : (" << dk << ","<< dt << "," << ds << ")" << std::endl;
            dist_source = zero;
            tmp3d_nospin = zero;
            evec3d = zero;
            for (int it = dt; it < Nt; it += TI){
              if( it >= Ntfirst && it < Ntfirst + Ntlocal ) {
                for (int ik = dk; ik < nvec; ik += LI){
                  for (int is = ds; is < Ns; is += Ns){ //at the moment, full spin dilution is enforced
                    std::cout <<  "LapH source vector from noise " << it << " and dilution component (d_k,d_t,d_alpha) : (" << ik << ","<< is << ")" << std::endl;
                    ExtractSliceLocal(evec3d,eig4d.evec[ik],0,it,3);
                    tmp3d_nospin = evec3d * noises[inoise][it][ik]()(is)(); //noises do not have to be a spin vector
                    tmp3d=zero;
                    pokeSpin(tmp3d,tmp3d_nospin,is);
                    tmp2=zero;
#ifdef USE_LOCAL_SLICES
                    InsertSliceLocal(tmp3d,tmp2,0,it-Ntfirst,Grid::QCD::Tdir);
#else
                    InsertSlice(tmp3d,tmp2,it,Grid::QCD::Tdir);
#endif
                    dist_source += tmp2;
                  }
                }
              }
            }
            std::cout <<  "Inversion for noise " << inoise << " and dilution component (d_k,d_t,d_alpha) : (" << dk << ","<< dt << "," << ds << ")" << std::endl;
            result=zero;
            LatticeFermion src5(FGrid);
            LatticeFermion sol5(FGrid);
            Dop.ImportPhysicalFermionSource(dist_source,src5);
            SchurSolver(Dop,src5,sol5);
            Dop.ExportPhysicalFermionSolution(sol5,result); //These are the meson sinks
            if (compute_current_sink)
              current_sink[inoise+nnoise*(dk+LI*(dt+Nt_inv*ds))] = result;
            std::cout <<  "Contraction of perambulator from noise " << inoise << " and dilution component (d_k,d_t,d_alpha) : (" << dk << ","<< dt << "," << ds << ")" << std::endl;
            for (int is = 0; is < Ns; is++) {
              result_nospin = peekSpin(result,is);
              for (int t = Ntfirst; t < Ntfirst + Ntlocal; t++) {
#ifdef USE_LOCAL_SLICES
                ExtractSliceLocal(result_3d,result_nospin,0,t-Ntfirst,Grid::QCD::Tdir);
#else
                ExtractSlice(result_3d,result_nospin,t,3);
#endif
                for (int ivec = 0; ivec < nvec; ivec++) {
                  ExtractSliceLocal(evec3d,eig4d.evec[ivec],0,t,3);
                  pokeSpin(perambulator(t, ivec, dk, inoise,dt,ds),innerProduct(evec3d, result_3d),is);
                }
          }
        }
      }
    }
    std::cout <<  "perambulator done" << std::endl;
    perambulator.SliceShare( grid3d, grid4d );

    // THIS IS WHERE WE WANT TO SAVE THE PERAMBULATORS TO DISK
    perambulator.WriteTemporary(std::string(pszPerambPack));
*/
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MDistil_perambulator_l_hpp_
