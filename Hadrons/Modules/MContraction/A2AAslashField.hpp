/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MContraction/A2AAslashField.hpp

Copyright (C) 2015-2018

Author: Antonin Portelli <antonin.portelli@me.com>

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
#ifndef Hadrons_MContraction_A2AAslashField_hpp_
#define Hadrons_MContraction_A2AAslashField_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/A2AMatrix.hpp>

#ifndef ASF_IO_TYPE
#define ASF_IO_TYPE ComplexF
#endif

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         A2AAslashField                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class A2AAslashFieldPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2AAslashFieldPar,
                                    int, cacheBlock,
                                    int, block,
                                    std::string, left,
                                    std::string, right,
                                    std::string, output,
                                    std::vector<std::string>, emField);
};

class A2AAslashFieldMetadata: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2AAslashFieldMetadata,
                                    std::string, emFieldName);
};

template <typename T, typename FImpl>
class AslashFieldKernel: public A2AKernel<T, typename FImpl::FermionField>
{
public:
    typedef typename FImpl::FermionField FermionField;
public:
    AslashFieldKernel(const std::vector<LatticeComplex> &emB0,
                      const std::vector<LatticeComplex> &emB1,
                      GridBase *grid)
    : emB0_(emB0), emB1_(emB1), grid_(grid)
    {
        vol_ = 1.;
        for (auto &d: grid_->GlobalDimensions())
        {
            vol_ *= d;
        }
    }

    virtual ~AslashFieldKernel(void) = default;
    virtual void operator()(A2AMatrixSet<T> &m, const FermionField *left, 
                            const FermionField *right,
                            const unsigned int orthogDim, double &t)
    {
        A2Autils<FImpl>::AslashField(m, left, right, emB0_, emB1_, orthogDim, &t);
    }

    virtual double flops(const unsigned int blockSizei, const unsigned int blockSizej)
    {
        return 0.;
    }

    virtual double bytes(const unsigned int blockSizei, const unsigned int blockSizej)
    {
        return 0.;
    }
private:
    const std::vector<LatticeComplex> &emB0_, &emB1_;
    GridBase                          *grid_;
    double                            vol_;
};

template <typename FImpl, typename PhotonImpl>
class TA2AAslashField: public Module<A2AAslashFieldPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    typedef typename PhotonImpl::GaugeField EmField;
    typedef A2AMatrixBlockComputation<Complex, 
                                      FermionField, 
                                      A2AAslashFieldMetadata, 
                                      ASF_IO_TYPE> Computation;
    typedef AslashFieldKernel<Complex, FImpl> Kernel;
public:
    // constructor
    TA2AAslashField(const std::string name);
    // destructor
    virtual ~TA2AAslashField(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(A2AAslashField, ARG(TA2AAslashField<FIMPL, PhotonR>), MContraction);

/******************************************************************************
 *                 TA2AAslashField implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename PhotonImpl>
TA2AAslashField<FImpl, PhotonImpl>::TA2AAslashField(const std::string name)
: Module<A2AAslashFieldPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename PhotonImpl>
std::vector<std::string> TA2AAslashField<FImpl, PhotonImpl>::getInput(void)
{
    std::vector<std::string> in = par().emField;
    
    in.push_back(par().left);
    in.push_back(par().right);

    return in;
}

template <typename FImpl, typename PhotonImpl>
std::vector<std::string> TA2AAslashField<FImpl, PhotonImpl>::getOutput(void)
{
    std::vector<std::string> out = {};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, typename PhotonImpl>
void TA2AAslashField<FImpl, PhotonImpl>::setup(void)
{
    envTmp(Computation, "computation", 1, envGetGrid(FermionField), 
           env().getNd() - 1, par().emField.size(), 1, par().block, 
           par().cacheBlock, this);
    envTmp(std::vector<ComplexField>, "B0", 1, 
           par().emField.size(), envGetGrid(ComplexField));
    envTmp(std::vector<ComplexField>, "B1", 1, 
           par().emField.size(), envGetGrid(ComplexField));
    envTmpLat(ComplexField, "Amu");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, typename PhotonImpl>
void TA2AAslashField<FImpl, PhotonImpl>::execute(void)
{
    auto &left  = envGet(std::vector<FermionField>, par().left);
    auto &right = envGet(std::vector<FermionField>, par().right);

    int nt         = env().getDim().back();
    int N_i        = left.size();
    int N_j        = right.size();
    int nem        = par().emField.size();
    int block      = par().block;
    int cacheBlock = par().cacheBlock;

    LOG(Message) << "Computing all-to-all A-slash fields" << std::endl;
    LOG(Message) << "Left: '" << par().left << "' Right: '" << par().right << "'" << std::endl;
    LOG(Message) << "EM fields:" << std::endl;
    for (auto &name: par().emField)
    {
        LOG(Message) << "  " << name << std::endl;
    }
    LOG(Message) << "A-slash field size: " << nt << "*" << N_i << "*" << N_j 
                 << " (filesize " << sizeString(nt*N_i*N_j*sizeof(ASF_IO_TYPE)) 
                 << "/EM field)" << std::endl;
    
    // preparing "B" complexified fields
    startTimer("Complexify EM fields");
    envGetTmp(std::vector<ComplexField>, B0);
    envGetTmp(std::vector<ComplexField>, B1);
    for (unsigned int i = 0; i < par().emField.size(); ++i)
    {
        auto &A = envGet(EmField, par().emField[i]);
        envGetTmp(ComplexField, Amu);

        B0[i]  = peekLorentz(A, 0);
        B0[i] += timesI(peekLorentz(A, 1));
        B1[i]  = peekLorentz(A, 2);
        B1[i] += timesI(peekLorentz(A, 3));
    }
    stopTimer("Complexify EM fields");

    // I/O name & metadata lambdas
    auto ionameFn = [this](const unsigned int em, const unsigned int dummy)
    {
        return par().emField[em];
    };

    auto filenameFn = [this, &ionameFn](const unsigned int em, const unsigned int dummy)
    {
        return par().output + "." + std::to_string(vm().getTrajectory()) 
               + "/" + ionameFn(em, dummy) + ".h5";
    };

    auto metadataFn = [this](const unsigned int em, const unsigned int dummy)
    {
        A2AAslashFieldMetadata md;

        md.emFieldName = par().emField[em];
        
        return md;
    };

    // executing computation
    Kernel kernel(B0, B1, envGetGrid(FermionField));

    envGetTmp(Computation, computation);
    computation.execute(left, right, kernel, ionameFn, filenameFn, metadataFn);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_A2AAslashField_hpp_
