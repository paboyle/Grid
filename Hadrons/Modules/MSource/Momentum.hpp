#ifndef Hadrons_MomSource_hpp_
#define Hadrons_MomSource_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/* 
Plane Wave source
-----------------
src_x = e^i2pi/L * p *position
*/

/******************************************************************************
 *                                  TPlane                                     *
  ******************************************************************************/
  BEGIN_MODULE_NAMESPACE(MSource)

  class MomSourcePar: Serializable
  {
  public:
//What is meant by serializable in this context
GRID_SERIALIZABLE_CLASS_MEMBERS(MomSourcePar,
std::string, mom);
};


template <typename FImpl>
class TMomSource: public Module<MomSourcePar>
{
public:
FERM_TYPE_ALIASES(FImpl,);
public:
// constructor
TMomSource(const std::string name);
// destructor
virtual ~TMomSource(void) = default;
// dependency relation
virtual std::vector<std::string> getInput(void);
virtual std::vector<std::string> getOutput(void);
// setup
virtual void setup(void);
// execution
virtual void execute(void);
};

MODULE_REGISTER_TMP(MomSource, TMomSource<FIMPL>, MSource);
//MODULE_REGISTER_NS(MomSource, TMomSource, MSource);

/******************************************************************************
*                       TMomSource template implementation                       *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TMomSource<FImpl>::TMomSource(const std::string name)
: Module<MomSourcePar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TMomSource<FImpl>::getInput(void)
{
    std::vector<std::string> in;
    return in;
}

template <typename FImpl>
std::vector<std::string> TMomSource<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    return out;
}


// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TMomSource<FImpl>::setup(void)
{
    envCreateLat(PropagatorField, getName());
}


//execution//////////////////////////////////////////////////////////////////
template <typename FImpl>
void TMomSource<FImpl>::execute(void)
{
    LOG(Message) << "Generating planewave momentum source with momentum " << par().mom << std::endl;
    //what does this env do?
    PropagatorField &src = envGet(PropagatorField, getName());
    Lattice<iScalar<vInteger>> t(env().getGrid());
    LatticeComplex             C(env().getGrid()), coor(env().getGrid());
    std::vector<Real>          p;
    std::vector<Real> latt_size(GridDefaultLatt().begin(), GridDefaultLatt().end()); 
    Complex                    i(0.0,1.0);

    LOG(Message) << " " << std::endl;
    //get the momentum from parameters
    p  = strToVec<Real>(par().mom);
    C = zero;
    LOG(Message) << "momentum converted from string - " << std::to_string(p[0]) <<std::to_string(p[1]) <<std::to_string(p[2]) <<   std::to_string(p[3]) << std::endl;
    for(int mu=0;mu<4;mu++){
    Real TwoPiL =  M_PI * 2.0/ latt_size[mu];
    LatticeCoordinate(coor,mu);
    C = C +(TwoPiL * p[mu]) * coor;
    }
    C = exp(C*i);
    LOG(Message) << "exponential of pdotx taken " << std::endl;
    src = src + C;
    LOG(Message) << "source created" << std::endl;

}



END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MomSource_hpp_
