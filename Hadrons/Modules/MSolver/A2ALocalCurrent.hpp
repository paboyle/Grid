#ifndef Hadrons_MSolver_A2ALocalCurrent_hpp_
#define Hadrons_MSolver_A2ALocalCurrent_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Solver.hpp>
#include <Hadrons/A2AVectors.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         A2ALocalCurrent                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSolver)

class A2ALocalCurrentPar : Serializable
{
  public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2ALocalCurrentPar,
                                    std::string, vector,
                                    unsigned int, tA,
                                    unsigned int, tB,
                                    Gamma::Algebra, gamma,
                                    std::string, solver,
                                    std::string, output,
                                    bool, multiFile);
};

template <typename FImpl>
class TA2ALocalCurrent : public Module<A2ALocalCurrentPar>
{
  public:
    FERM_TYPE_ALIASES(FImpl, );
    SOLVER_TYPE_ALIASES(FImpl, );

  public:
    // constructor
    TA2ALocalCurrent(const std::string name);
    // destructor
    virtual ~TA2ALocalCurrent(void){};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);

  private:
    unsigned int Ls_;
    std::string tName_;
};

MODULE_REGISTER_TMP(A2ALocalCurrent, TA2ALocalCurrent<FIMPL>, MSolver);

/******************************************************************************
 *                 TA2ALocalCurrent implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TA2ALocalCurrent<FImpl>::TA2ALocalCurrent(const std::string name)
    : Module<A2ALocalCurrentPar>(name)
{
}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TA2ALocalCurrent<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().vector, par().solver};

    return in;
}

template <typename FImpl>
std::vector<std::string> TA2ALocalCurrent<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2ALocalCurrent<FImpl>::setup(void)
{
    Ls_ = env().getObjectLs(par().solver);
    auto &vvector = envGet(std::vector<FermionField>, par().vector);
    unsigned int Nmodes = vvector.size();
    envCreate(std::vector<FermionField>, getName(), 1,
              Nmodes, envGetGrid(FermionField));

    envCache(Lattice<iScalar<vInteger>>, tName_, 1, envGetGrid(LatticeComplex));
    envTmpLat(FermionField, "v4dtmp");
    envTmpLat(FermionField, "v5dtmp", Ls_);
    envTmpLat(FermionField, "v5dtmp_sol", Ls_);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2ALocalCurrent<FImpl>::execute(void)
{

    if (par().tA == par().tB)
    {
        LOG(Message) << "Generating Sequential propagator: Gmu * v, with Gmu = "
                     << par().gamma
                     << " and the A2A vector "
                     << par().vector << " at t= " << par().tA << std::endl;
    }
    else
    {
        LOG(Message) << "Generating Sequential propagator: Gmu * v, with Gmu = "
                     << par().gamma
                     << " and the A2A vector "
                     << par().vector << " for "
                     << par().tA << " <= t <= " << par().tB << std::endl;
    }

    auto &solver = envGet(Solver, par().solver);
    auto &vvector = envGet(std::vector<FermionField>, par().vector);
    auto &gmuV = envGet(std::vector<FermionField>, getName());
    unsigned int Nmodes = vvector.size();
    auto &mat = solver.getFMat();
    auto &t = envGet(Lattice<iScalar<vInteger>>, tName_);
    envGetTmp(FermionField, v4dtmp);
    envGetTmp(FermionField, v5dtmp);
    envGetTmp(FermionField, v5dtmp_sol);
    Gamma gmu(par().gamma);

    startTimer("Seq Local Current");
    LOG(Message) << "Calculate Sequential propagator on Gmu * v with Gmu "
                 << par().gamma
                 << " and the A2A vector "
                 << par().vector << "." << std::endl;

    for (int i = 0; i < Nmodes; i++)
    {
        v4dtmp = zero;
        startTimer("Multiply Local Current");
        v4dtmp = where((t >= par().tA) and (t <= par().tB), gmu * vvector[i], 0. * vvector[i]);
        stopTimer("Multiply Local Current");

        LOG(Message) << "Gmu * v vector i =  " << i << std::endl;
        startTimer("Inversion");
        if (Ls_ == 1)
        {
            solver(gmuV[i], v4dtmp);
        }
        else
        {
            mat.ImportPhysicalFermionSource(v4dtmp, v5dtmp);
            solver(v5dtmp_sol, v5dtmp);
            mat.ExportPhysicalFermionSolution(v5dtmp_sol, v4dtmp);
            gmuV[i] = v4dtmp;
        }
        stopTimer("Inversion");
    }
    stopTimer("Seq Local Current");
    if (!par().output.empty())
    {
        startTimer("I/O");
        A2AVectorsIo::write(par().output, gmuV, par().multiFile, vm().getTrajectory());
        stopTimer("I/O");
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSolver_A2ALocalCurrent_hpp_
