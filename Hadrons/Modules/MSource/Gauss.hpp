#ifndef Hadrons_MSource_Gauss_hpp_
#define Hadrons_MSource_Gauss_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         Gauss                                              *
 * result[n] = 1/(sqrt(2*pi)*width)^dim*exp(-|n|^2/(2*width^2))               *
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
    
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TGauss<FImpl>::execute(void)
{
    auto &rho = envGet(LatticeComplex, getName());
    envGetTmp(LatticeComplex, component);
    const int dim=env().getNd()-1;
    const double fact=-0.5/std::pow(par().width,2);
    const std::vector<int> latt_size { env().getGrid()->FullDimensions() };

    //exp(fact*|n|^2)
    rho=zero;
    for(int mu=0; mu<dim; mu++) {
        LatticeCoordinate(component, mu);
        //FIXME: the next three lines are very inefficient...
        //       should not need any communication (Cshift) here
        assert(latt_size[mu]%2==0);
        component-=Complex(latt_size[mu]/2-1);
        component=Cshift(component, mu, latt_size[mu]/2-1 -position_[mu]);
        rho+=component*component*fact;
    }
    rho=exp(rho);

    rho*=static_cast<Complex>(std::pow(sqrt(2*M_PI)*par().width,dim));
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSource_Gauss_hpp_
