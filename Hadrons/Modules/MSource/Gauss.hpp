#ifndef Hadrons_MSource_Gauss_hpp_
#define Hadrons_MSource_Gauss_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         Gauss                                              *
 * result[n] = 1/(sqrt(2*pi)*width)^dim                                       *
 *            * exp(-|n-position|^2/(2*width^2) + i p.n)                      *
 * where:                                                                     *
 *   n=(n[0],n[1],...,n[dim-1])  (lattice coordinate)                         *
 *   p[i]=2*pi/L[i]*mom[i]                                                    *
 *   dim=Nd-1                                                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSource)

class GaussPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(GaussPar,
                                    std::string, position,
                                    std::string, mom,
                                    double,      width);
};

template <typename FImpl>
class TGauss: public Module<GaussPar>
{
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
    auto parse_vector = [](const std::string &momentum, int dim,
            const std::string &desc)
    {
        std::vector<int> res = strToVec<int>(momentum);
        if(res.size() != dim) {
            HADRONS_ERROR(Size, desc + " has " + std::to_string(res.size())
                    + " instead of " + std::to_string(dim) + " components");
        }
        return res;
    };
    position_ = parse_vector(par().position, env().getNd()-1, "position");
    mom_      = parse_vector(par().mom,      env().getNd()-1, "momentum");

    envCreateLat(LatticeComplex, getName());
    envTmpLat(LatticeComplex, "component");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TGauss<FImpl>::execute(void)
{
    auto &rho = envGet(LatticeComplex, getName());
    envGetTmp(LatticeComplex, component);
    const int dim=env().getNd()-1;
    const double fact=-0.5/std::pow(par().width,2);
    const Complex i(0.0, 1.0);

    rho=zero;
    for(int mu=0; mu<dim; mu++) {
        LatticeCoordinate(component, mu);
        rho+=(i*(mom_[mu]*2*M_PI/env().getDim(mu)))*component;
        //FIXME: the next three lines are very inefficient...
        //       should not need any communication (Cshift) here
        assert(env().getDim(mu)%2==0);
        component-=Complex(env().getDim(mu)/2-1);
        component=Cshift(component, mu, env().getDim(mu)/2-1 -position_[mu]);
        rho+=component*component*fact;
    }
    rho=exp(rho);

    rho*=static_cast<Complex>(std::pow(sqrt(2*M_PI)*par().width,dim));
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSource_Gauss_hpp_
