#ifndef Hadrons_MDistil_BContraction_hpp_
#define Hadrons_MDistil_BContraction_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Solver.hpp>
#include <Hadrons/EigenPack.hpp>
#include <Hadrons/A2AVectors.hpp>
#include <Hadrons/DilutedNoise.hpp>

// These are members of Distillation
#include <Hadrons/Modules/MDistil/Distil.hpp>
BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         BContraction                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MDistil)

class BContractionPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(BContractionPar,
		                    std::string, one,
		                    std::string, two,
		                    std::string, three,
		                    std::string, output,
                                    std::vector<std::string>, mom);
};

template <typename FImpl>
class TBContraction: public Module<BContractionPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TBContraction(const std::string name);
    // destructor
    virtual ~TBContraction(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    bool                               hasPhase_{false};
    std::string                        momphName_;
    std::vector<Gamma::Algebra>        gamma12_;
    std::vector<Gamma::Algebra>        gamma23_;
    std::vector<std::vector<Real>>     mom_;
protected:
    GridCartesian * grid4d;
    GridCartesian * grid3d;
};

MODULE_REGISTER_TMP(BContraction, TBContraction<FIMPL>, MDistil);

/******************************************************************************
 *                 TBContraction implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TBContraction<FImpl>::TBContraction(const std::string name)
: Module<BContractionPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TBContraction<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().one, par().two, par().three};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TBContraction<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TBContraction<FImpl>::setup(void)
{
    
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TBContraction<FImpl>::execute(void)
{
    auto &one   = envGet(std::vector<FermionField>, par().one);
    auto &two   = envGet(std::vector<FermionField>, par().two);
    auto &three = envGet(std::vector<FermionField>, par().three);

    int N_1     = one.size();
    int N_2     = two.size();
    int N_3     = three.size();

    LOG(Message) << "Computing distillation baryon fields" << std::endl;
    LOG(Message) << "One: '" << par().one << "' Two: '" << par().two  << "' Three: '" << par().three << "'" << std::endl;
    LOG(Message) << "Momenta:" << std::endl;
    for (auto &p: mom_)
    {
        LOG(Message) << "  " << p << std::endl;
    }

    grid4d = env().getGrid();
    grid3d = MakeLowerDimGrid(grid4d);    
    int Nmom=1;
    int Nt=64;
    std::vector<Complex> BField(Nmom*Nt*N_1*N_2*N_3);    
    int Bindex;
    int Nc=3; //Num colours

    FermionField tmp1(grid3d);
    FermionField tmp2(grid3d);
    FermionField tmp3(grid3d);
    //std::complex<double> * tmp33 = reinterpret_cast<std::complex<double> *>(&(tmp3[0]()(0)(0)));

    SpinColourVector * tmp11 = reinterpret_cast<SpinColourVector *>(&(tmp1[0]()(0)(0)));
    SpinColourVector * tmp22 = reinterpret_cast<SpinColourVector *>(&(tmp2[0]()(0)(0)));
    SpinColourVector * tmp33 = reinterpret_cast<SpinColourVector *>(&(tmp3[0]()(0)(0)));
   
    SpinVector tmp33s; 
    SpinColourVector * tmp333;
    SpinColourVector * diquark;
    SpinColourVector * tmp222;

    std::vector<std::vector<int>> epsilon = {{0,1,2},{1,2,0},{2,0,1},{0,2,1},{2,1,0},{1,0,2}};
    std::vector<int> epsilon_sgn = {1,1,1,-1,-1,-1};

    gamma12_ = {
       Gamma::Algebra::Identity, // I
       Gamma::Algebra::Gamma5, // gamma_5
       Gamma::Algebra::Identity, // I
    };
    gamma23_ = { // C = i gamma_2 gamma_4
       Gamma::Algebra::SigmaXZ, // C gamma_5 = -i gamma_1 gamma_3
       Gamma::Algebra::SigmaYT, // C = i gamma_2 gamma_4
       Gamma::Algebra::GammaYGamma5, // i gamma_4 C gamma_5 = i gamma_2 gamma_5
    };
    //std::vector<Complex> factor23 = {(0.,-1.),(0.,1.),(0.,1.)};

  //  SpinColourVector a;
  //  SpinVector b;
  //  b = peekColour(a,0);
    //b= a()(0)();
    //tmp33s = peekColour(a,0);

    for (int i1=0 ; i1 < N_1 ; i1++){
      for (int i2=0 ; i2 < N_2 ; i2++){
        for (int i3=0 ; i3 < N_3 ; i3++){
          for (int imom=0 ; imom < Nmom ; imom++){
            Bindex = i1 + N_1*(i2 + N_2*(i3 + N_3*imom));
            for (int t=0 ; t < Nt ; t++){
	      ExtractSliceLocal(tmp1,one[i1],0,t,3);
              parallel_for (unsigned int sU = 0; sU < grid3d->oSites(); ++sU)
              {
                for (int ie=0 ; ie < 6 ; ie++){
	          //tmp33s = peekColour(tmp33[sU],epsilon[ie][2]);
		  //tmp333 = Gamma(gamma23_[0])*tmp33s;
		  //tmp333 = Gamma(gamma23_[0])*tmp33[sU]()()(epsilon[ie][2]);
		  /*diquark = tmp22[sU]()()(epsilon[ie][1])*factor23[0]*tmp333;
		  tmp222 = Gamma(gamma12_[0])*diquark;
                  BField[Bindex]+=(double)epsilon_sgn[ie]*(tmp11[sU]()()(epsilon[ie][0])*tmp222);
                  */
  		  }
  	      }
            }
	  }
	}
      }
    }

}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MDistil_BContraction_hpp_
