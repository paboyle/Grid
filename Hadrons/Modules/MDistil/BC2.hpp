/*************************************************************************************
 
 Grid physics library, www.github.com/paboyle/Grid
 
 Source file: Hadrons/Modules/MDistil/BC2.hpp
 
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

#ifndef Hadrons_MDistil_BC2_hpp_
#define Hadrons_MDistil_BC2_hpp_

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
 *                         BC2                                 *
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
    //   6 - k   - right  distillation mode index
    // template <typename T>
    // using BaryonTensorSet = Eigen::TensorMap<Eigen::Tensor<T, 7, Eigen::RowMajor>>;

class BC2Par: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(BC2Par,
		                    std::string, one,
		                    std::string, two,
		                    std::string, three,
		                    std::string, output,
				    int, parity,
                                    std::vector<std::string>, mom);
};

template <typename FImpl>
class TBC2: public Module<BC2Par>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TBC2(const std::string name);
    // destructor
    virtual ~TBC2(void) {};
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

MODULE_REGISTER_TMP(BC2, TBC2<FIMPL>, MDistil);

/******************************************************************************
 *                 TBC2 implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TBC2<FImpl>::TBC2(const std::string name)
: Module<BC2Par>(name)
, momphName_(name + "_momph")
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TBC2<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().one, par().two, par().three};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TBC2<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TBC2<FImpl>::setup(void)
{
  if(!mom_.size()) {
       for (auto &pstr: par().mom)
       {
         auto p = strToVec<Real>(pstr);

         if (p.size() != env().getNd() - 1)
         {
           HADRONS_ERROR(Size, "Momentum has " + std::to_string(p.size()) + " components instead of " + std::to_string(env().getNd() - 1));
	 }
	 mom_.push_back(p);
    }
  }
    //envCache(std::vector<ComplexField>, momphName_, 1, par().mom.size(), envGetGrid(ComplexField));
  static GridCartesian * MyGrid{env().getGrid()};
  if( MyGrid == envGetGrid(ComplexField) )
    LOG(Message) << "envGetGrid(ComplexField) == env().getGrid()" << std::endl;
  else
    LOG(Message) << "envGetGrid(ComplexField) != env().getGrid()" << std::endl;
  envTmp(std::vector<ComplexField>, "ph", 1, std::vector<ComplexField>());
  envGetTmp(std::vector<ComplexField>, ph);
  if(!ph.size()) {
    for (unsigned int j = 0; j < par().mom.size(); ++j)
      ph.push_back(ComplexField(MyGrid));
  }

    envTmpLat(ComplexField, "coor");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TBC2<FImpl>::execute(void)
{
    auto &one   = envGet(std::vector<FermionField>, par().one);
    auto &two   = envGet(std::vector<FermionField>, par().two);
    auto &three = envGet(std::vector<FermionField>, par().three);
    const std::string &output{par().output};

    int N_1 = static_cast<int>(one.size());
    int N_2 = static_cast<int>(two.size());
    int N_3 = static_cast<int>(three.size());

    LOG(Message) << "Computing distillation baryon fields" << std::endl;
    LOG(Message) << "One: '" << par().one << "' Two: '" << par().two  << "' Three: '" << par().three << "'" << std::endl;
    LOG(Message) << "Momenta:" << std::endl;
    for (auto &p: mom_)
    {
        LOG(Message) << "  " << p << std::endl;
    }


    int Nmom = static_cast<int>(mom_.size());
    const int Nt{env().getDim(Tdir)};

    int parity = 1;
    int orthogDim=3;

    //auto &ph = envGet(std::vector<ComplexField>, momphName_);
  envGetTmp(std::vector<ComplexField>, ph);

    if (!hasPhase_)
    {
        startTimer("Momentum phases");
        for (unsigned int j = 0; j < Nmom; ++j)
        {
            Complex           i(0.0,1.0);
            std::vector<Real> p;

            envGetTmp(ComplexField, coor);
            ph[j] = zero;
            for(unsigned int mu = 0; mu < mom_[j].size(); mu++)
            {
                LatticeCoordinate(coor, mu);
                ph[j] = ph[j] + (mom_[j][mu]/env().getDim(mu))*coor;
            }
            ph[j] = exp((Real)(2*M_PI)*i*ph[j]);
        }
        hasPhase_ = true;
        stopTimer("Momentum phases");
    }
    //envCache(std::vector<ComplexField>, momphName_, 1, mom_.size(), envGetGrid(ComplexField));

    Eigen::Tensor<ComplexD, 6> m(Nmom,Nt,N_1,N_2,N_3,4);
    //A2Autils<FImpl>::NucleonFieldMom(m, &one[0], &two[0], &three[0], ph, parity, orthogDim);
    A2Autils<FImpl>::NucleonFieldMom(m, one, two, three, ph, parity, orthogDim);
    for (int is=0 ; is < 4 ; is++){
      for (int t=0 ; t < Nt ; t++){
        std::cout << "BaryonField(is=" << is << ",t=" << t << ") = " << m(0,t,0,0,0,is) << std::endl;
      }
    }
    
    BFieldIO BField_save;
    BField_save.BField = m;


    GridCartesian * grid = env().getGrid();
    if(grid->IsBoss()) {
      std::string filename ="./" + output + ".h5"; 
      std::cout << "Writing to file " << filename << std::endl;
      Grid::Hdf5Writer writer(filename);
      write(writer,"BaryonField",BField_save.BField);
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MDistil_BC2_hpp_
