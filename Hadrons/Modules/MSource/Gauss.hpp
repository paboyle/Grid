#ifndef Hadrons_MSource_Gauss_hpp_
#define Hadrons_MSource_Gauss_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         Gauss                                              *
 * result[n] = 1/(sqrt(2*pi)*width)^dim                                       *
 *            * exp(-|n-position|^2/(2*width^2))                              *
 *            * exp(i*2*pi/L*mom*n)                                           *
 * where:                                                                     *
 *   n=(n[0],n[1],...,n[dim-1])  (lattice coordinate)                         *
 *   dim=Nd-1                                                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSource)

class GaussPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(GaussPar,
                                    std::string, position,
                                    std::string, mom,
                                    Integer,     tA,
                                    Integer,     tB,
                                    double,      width);
};

template <typename FImpl>
class TGauss: public Module<GaussPar>
{
    BASIC_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TGauss(const std::string name);
    // destructor
    virtual ~TGauss(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    std::vector<int> position_;
    std::vector<int> mom_;
};

MODULE_REGISTER_TMP(Gauss, TGauss<FIMPL>, MSource);
MODULE_REGISTER_TMP(ScalarGauss, TGauss<ScalarImplCR>, MSource);

/******************************************************************************
 *                 TGauss implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TGauss<FImpl>::TGauss(const std::string name)
: Module<GaussPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TGauss<FImpl>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TGauss<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TGauss<FImpl>::setup(void)
{
    auto parse_vector = [](const std::string &vec, int dim,
            const std::string &desc)
    {
        std::vector<int> res = strToVec<int>(vec);
        if(res.size() != dim) {
            HADRONS_ERROR(Size, desc + " has "
                    + std::to_string(res.size()) + " instead of "
                    + std::to_string(dim) + " components");
        }
        return res;
    };
    position_ = parse_vector(par().position, env().getNd()-1, "position");
    mom_      = parse_vector(par().mom,      env().getNd(),   "momentum");

    envCreateLat(PropagatorField, getName());
    envTmpLat(ComplexField, "component");
    envTmpLat(ComplexField, "ScalarRho");
    envTmp(LatticeInteger, "compHelper", 1, envGetGrid(ComplexField));
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TGauss<FImpl>::execute(void)
{
    auto &rho = envGet(PropagatorField, getName());
    envGetTmp(ComplexField, component);
    envGetTmp(ComplexField, ScalarRho);
    envGetTmp(LatticeInteger, compHelper);
    const int dim=env().getNd()-1;
    const Real fact=-0.5/std::pow(par().width,2);
    const Complex i(0.0, 1.0);
    const Real Pi(M_PI);
    const SitePropagator idMat=[](){ SitePropagator s; s=1.; return s; }();

    ScalarRho=Zero();
    for(int mu=0; mu<dim; mu++) {
        assert(env().getDim(mu)%2==0);
        assert(position_[mu]>=0 && position_[mu]<env().getDim(mu));

        const int Lmu=env().getDim(mu);
        const int LmuHalf=Lmu/2;
        const int posMu=position_[mu];

        LatticeCoordinate(component, mu);
        LatticeCoordinate(compHelper, mu);

        //spatial dimensions of momentum phase
        ScalarRho+=(i*(mom_[mu]*2*Pi/Lmu))*component;

        //Gauss distribution
        component-=Complex(posMu);
        if(posMu<LmuHalf)
        {
            component=where((compHelper>Integer(posMu+LmuHalf)),
                    component-Complex(Lmu),
                    component);
        }
        else
        {
            component=where((compHelper<=Integer(posMu-LmuHalf)),
                    component+Complex(Lmu),
                    component);
        }
        ScalarRho+=component*component*fact;
    }

    //time component of momentum phase
    LatticeCoordinate(component, dim);
    ScalarRho+=(i*(mom_.at(dim)*2*Pi/env().getDim(dim)))*component;

    //compute scalar result
    ScalarRho=exp(ScalarRho)*Complex(std::pow(sqrt(2*Pi)*par().width,-dim));

    //select time slices
    LatticeCoordinate(compHelper, dim);
    ScalarRho=where((compHelper>=par().tA && compHelper<=par().tB),
          ScalarRho,
          0.*ScalarRho);

    //compute output field rho
    rho=ScalarRho*idMat;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSource_Gauss_hpp_
