#ifndef Hadrons_MDistil_Baryon2pt_hpp_
#define Hadrons_MDistil_Baryon2pt_hpp_

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
 *                         Baryon2pt                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MDistil)

class Baryon2ptPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(Baryon2ptPar,
		                    std::string, inputL,
		                    std::string, inputR,
		                    std::string, quarksL,
		                    std::string, quarksR,
		                    std::string, output
                                    );
};

template <typename FImpl>
class TBaryon2pt: public Module<Baryon2ptPar>
{
public:
    // constructor
    TBaryon2pt(const std::string name);
    // destructor
    virtual ~TBaryon2pt(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

class C2IO: Serializable{
public:
using C2Set = Eigen::Tensor<Complex, 2>;
  GRID_SERIALIZABLE_CLASS_MEMBERS(C2IO,
                                  C2Set, C2
		                  );
};

MODULE_REGISTER_TMP(Baryon2pt, TBaryon2pt<FIMPL>, MDistil);

/******************************************************************************
 *                 TBaryon2pt implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TBaryon2pt<FImpl>::TBaryon2pt(const std::string name)
: Module<Baryon2ptPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TBaryon2pt<FImpl>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TBaryon2pt<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TBaryon2pt<FImpl>::setup(void)
{
    
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TBaryon2pt<FImpl>::execute(void)
{

    const std::string &inputL{par().inputL};
    const std::string &inputR{par().inputR};
    const std::string &quarksL{par().quarksL};
    const std::string &quarksR{par().quarksR};
    const std::string &output{par().output};

    int Nmom=1;
    int Nt=64;
    int Nc=3; //Num colours
    std::vector<std::vector<int>> epsilon = {{0,1,2},{1,2,0},{2,0,1},{0,2,1},{2,1,0},{1,0,2}};
    std::vector<int> epsilon_sgn = {1,1,1,-1,-1,-1};
    int Ngamma=3;

    int N_1=20;
    int N_2=20;
    int N_3=20;

  //  using BaryonTensorSet = Eigen::Tensor<Complex, 7>;

    BFieldIO BFieldL;
    BFieldL.BField.resize(Nmom,Ngamma,Nt,4,N_1,N_2,N_3);

    std::string filenameL ="./" + inputL + ".h5"; 
    std::cout << "Reading from file " << filenameL << std::endl;
    Hdf5Reader readerL(filenameL);
    read(readerL,"BaryonField",BFieldL.BField);

    BFieldIO BFieldR;
    BFieldR.BField.resize(Nmom,Ngamma,Nt,4,N_1,N_2,N_3);

    std::string filenameR ="./" + inputR + ".h5"; 
    std::cout << "Reading from file " << filenameR << std::endl;
    Hdf5Reader readerR(filenameR);
    read(readerR,"BaryonField",BFieldR.BField);

    Eigen::Tensor<Complex, 2> corr(Nmom,Nt); 

    int Npairs = 0;
    char left[] = "uud";
    char right[] = "uud";
    std::vector<int> pairs(6);
    for (int ie=0, i=0 ; ie < 6 ; ie++){
      if (left[0] == right[epsilon[ie][0]] && left[1] == right[epsilon[ie][1]] && left[2] == right[epsilon[ie][2]]){
        pairs[i] = ie;
        i++;
        Npairs++;
      }
    }
    pairs.resize(Npairs);
    std::cout << Npairs << " pairs: " << pairs << std::endl;
    for (int imom=0 ; imom < Nmom ; imom++){
      for (int t=0 ; t < Nt ; t++){
        corr(imom,t) = 0.;
      }
    }

    int tsrc=0;

    for (int ipair=0 ; ipair < Npairs ; ipair++){
      Eigen::array<Eigen::IndexPair<int>, 3> product_dims = { Eigen::IndexPair<int>(0,epsilon[pairs[ipair]][0]),Eigen::IndexPair<int>(1,epsilon[pairs[ipair]][1]) ,Eigen::IndexPair<int>(2,epsilon[pairs[ipair]][2])  };
      for (int imom=0 ; imom < Nmom ; imom++){
        std::cout << imom << std::endl;
        Eigen::Tensor<Complex,6> B6L = BFieldL.BField.chip(imom,0);
        Eigen::Tensor<Complex,6> B6R = BFieldR.BField.chip(imom,0);
        for (int ig=0 ; ig < Ngamma ; ig++){
          Eigen::Tensor<Complex,5> B5L = B6L.chip(ig,0);
          Eigen::Tensor<Complex,5> B5R = B6R.chip(ig,0);
          for (int t=0 ; t < Nt ; t++){
            Eigen::Tensor<Complex,4> B4L = B5L.chip(t,0);
            Eigen::Tensor<Complex,4> B4R = B5R.chip(tsrc,0);
            for (int is=0 ; is < 4 ; is++){
              Eigen::Tensor<Complex,3> B3L = B4L.chip(is,0);
              Eigen::Tensor<Complex,3> B3R = B4R.chip(is,0);
              Eigen::Tensor<Complex,0> C2 = B3L.contract(B3R,product_dims);
              corr(imom,t) += static_cast<Real>(epsilon_sgn[pairs[ipair]])*C2(0);
            }
          }
        }
      }
    }
    for (int t=0 ; t < Nt ; t++){
      std::cout << "C2(t=" << t << ") = " << corr(0,t) << std::endl;
    }    

    C2IO C2_save;
    C2_save.C2 = corr;

    std::string filename ="./" + output + ".h5"; 
    std::cout << "Writing to file " << filename << std::endl;
    Hdf5Writer writer(filename);
    write(writer,"C2",C2_save.C2);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MDistil_Baryon2pt_hpp_
