/*************************************************************************************
 
 Grid physics library, www.github.com/paboyle/Grid
 
 Source file: Hadrons/Modules/MDistil/BContraction.hpp
 
 Copyright (C) 2019
 
 Author: Felix Erben <ferben@ed.ac.uk>
 Author: Michael Marshall <Michael.Marshall@ed.ac.uk>
 
 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License along
 with this program; if not, write to the Free Software Foundation, Inc.,
 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 
 See the full license in the file "LICENSE" in the top level distribution directory
 *************************************************************************************/
/*  END LEGAL */

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
#include <Hadrons/Distil.hpp>
BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         BContraction                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MDistil)

    // general baryon tensor set based on Eigen tensors and Grid-allocated memory
    // Dimensions:
    //   0 - ext - external field (momentum, EM field, ...)
    //   1 - str - dirac structure
    //   2 - t   - timeslice
    //   3 - s   - free spin index
    //   4 - i   - left  distillation mode index
    //   5 - j   - middle  distillation mode index
    //   6 - k   - left  distillation mode index
    // template <typename T>
    // using BaryonTensorSet = Eigen::TensorMap<Eigen::Tensor<T, 7, Eigen::RowMajor>>;


class BContractionPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(BContractionPar,
		                    std::string, one,
		                    std::string, two,
		                    std::string, three,
		                    std::string, output,
				    int, parity,
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

/*class BFieldIO: Serializable{
public:
  using BaryonTensorSet = Eigen::Tensor<Complex, 7>;
  GRID_SERIALIZABLE_CLASS_MEMBERS(BFieldIO,
                                  BaryonTensorSet, BField
		                  );
};*/

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

    int parity = par().parity;
    const std::string &output{par().output};

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
   // std::vector<Complex> BField(Nmom*Nt*N_1*N_2*N_3);    
    int Bindex;
    int Nc=3; //Num colours

    FermionField tmp1(grid3d);
    FermionField tmp2(grid3d);
    FermionField tmp3(grid3d);
    //std::complex<double> * tmp33 = reinterpret_cast<std::complex<double> *>(&(tmp3[0]()(0)(0)));

    SpinColourVector * tmp11 = reinterpret_cast<SpinColourVector *>(&(tmp1[0]()(0)(0)));
    SpinColourVector * tmp22 = reinterpret_cast<SpinColourVector *>(&(tmp2[0]()(0)(0)));
    SpinColourVector * tmp33 = reinterpret_cast<SpinColourVector *>(&(tmp3[0]()(0)(0)));
   
    SpinVector tmp11s; 
    SpinVector tmp22s; 
    SpinVector tmp33s; 
    SpinVector tmp333;
    SpinMatrix diquark;
    SpinMatrix g_diquark;
    SpinVector tmp222;
    SpinVector tmp111;


    assert(parity == 1 || parity == -1);

    std::vector<std::vector<int>> epsilon = {{0,1,2},{1,2,0},{2,0,1},{0,2,1},{2,1,0},{1,0,2}};
    std::vector<int> epsilon_sgn = {1,1,1,-1,-1,-1};

    Gamma g4(Gamma::Algebra::GammaT);

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
  std::vector<Complex> factor23{{0.,-1.},{0.,1.},{0.,1.}}; 
  using BaryonTensorSet = Eigen::Tensor<Complex, 7>;
  int Ngamma=3;
  BaryonTensorSet BField3(Nmom,Ngamma,Nt,4,N_1,N_2,N_3); 

  Eigen::Tensor<Complex, 3> corr(Nmom,4,Nt); 

   
    Complex diquark2;
    for (int i1=0 ; i1 < N_1 ; i1++){
      for (int i2=0 ; i2 < N_2 ; i2++){
        for (int i3=0 ; i3 < N_3 ; i3++){
          for (int imom=0 ; imom < Nmom ; imom++){
            for (int t=0 ; t < Nt ; t++){
              Bindex = i1 + N_1*(i2 + N_2*(i3 + N_3*(imom+Nmom*t)));
	      ExtractSliceLocal(tmp1,one[i1],0,t,3);
	      ExtractSliceLocal(tmp2,two[i2],0,t,3);
	      ExtractSliceLocal(tmp3,three[i3],0,t,3);
              parallel_for (unsigned int sU = 0; sU < grid3d->oSites(); ++sU)
              {
                for (int ie=0 ; ie < 6 ; ie++){
		  // Why does peekColour not work????
                  for (int is=0 ; is < 4 ; is++){
	            tmp11s()(is)() = tmp11[sU]()(is)(epsilon[ie][0]);
	            tmp22s()(is)() = tmp22[sU]()(is)(epsilon[ie][1]);
	            tmp33s()(is)() = tmp33[sU]()(is)(epsilon[ie][2]);
		  }
                  for (int ig=0 ; ig < Ngamma ; ig++){
		    tmp333 = Gamma(gamma23_[ig])*tmp33s;
		    tmp111 = Gamma(gamma12_[ig])*tmp11s;
		    tmp222 = g4*tmp111;
		    tmp111 = 0.5*(double)parity*(tmp111 + tmp222); // P_\pm * ...
                    diquark2 = factor23[0]*innerProduct(tmp22s,tmp333);
                    for (int is=0 ; is < 4 ; is++){
                      BField3(imom,ig,t,is,i1,i2,i3)+=static_cast<Real>(epsilon_sgn[ie])*tmp111()(is)()*diquark2;
                    }
  		  }
		}
  	      }
            }
	  }
	}
      }
    }
    for (int is=0 ; is < 4 ; is++){
      for (int t=0 ; t < Nt ; t++){
        std::cout << "BaryonField(is=" << is << ",t=" << t << ") = " << BField3(0,0,t,is,0,0,0) << std::endl;
      }
    }

    BFieldIO BField_save;
#ifdef FELIX_ISSUE
    BField_save.BField = BField3;
#endif
  std::string filename ="./" + output + ".h5"; 
  std::cout << "Writing to file " << filename << std::endl;
  Hdf5Writer writer(filename);
  write(writer,"BaryonField",BField_save.BField);

}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MDistil_BContraction_hpp_
