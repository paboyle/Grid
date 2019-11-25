/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MFermion/GaugeSpectProp.hpp

Copyright (C) 2015-2019

Author: T. Blum

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

#ifndef Hadrons_MFermion_GaugeSpectProp_hpp_
#define Hadrons_MFermion_GaugeSpectProp_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Solver.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                                GaugeSpectProp                                   *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MFermion)

class GaugeSpectPropPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(GaugeSpectPropPar,
                                    std::string, source,
                                    std::string, solver,
                                    std::string, left,
                                    std::string, right);
};

template <typename FImpl>
class TStagGaugeSpectProp: public Module<GaugeSpectPropPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    SOLVER_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TStagGaugeSpectProp(const std::string name);
    // destructor
    virtual ~TStagGaugeSpectProp(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    unsigned int Ls_;
    Solver       *solver_{nullptr};
};

MODULE_REGISTER_TMP(StagGaugeSpectProp, TStagGaugeSpectProp<STAGIMPL>, MFermion);

/******************************************************************************
 *                      TStagGaugeSpectProp implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TStagGaugeSpectProp<FImpl>::TStagGaugeSpectProp(const std::string name)
: Module<GaugeSpectPropPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TStagGaugeSpectProp<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().source,
        par().solver,
        par().left,
        par().right};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TStagGaugeSpectProp<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName(), getName() + "_5d"};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TStagGaugeSpectProp<FImpl>::setup(void)
{
    Ls_ = env().getObjectLs(par().solver);
    envCreateLat(PropagatorField, getName());
    envTmpLat(FermionField, "tmp");
    if (Ls_ > 1)
    {
        envTmpLat(FermionField, "source", Ls_);
        envTmpLat(FermionField, "sol", Ls_);
        envCreateLat(PropagatorField, getName() + "_5d", Ls_);
    }
    else
    {
        envTmpLat(FermionField, "source");
        envTmpLat(FermionField, "sol");
    }
}

//#include <Grid/tensors/Tensor_inner.h>
// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TStagGaugeSpectProp<FImpl>::execute(void)
{
    LOG(Message) << "Computing quark propagator from Spectral Decompostion'" << getName() << "'"
    << std::endl;
    
    std::string propName = (Ls_ == 1) ? getName() : (getName() + "_5d");

    auto &left  = envGet(std::vector<FermionField>, par().left);
    auto &right = envGet(std::vector<FermionField>, par().right);
    //int nt         = env().getDim().back();
    int N_i        = left.size();
    int N_j        = right.size();
    
    auto        &prop    = envGet(PropagatorField, propName);
    auto        &fullSrc = envGet(PropagatorField, par().source);
    auto        &solver  = envGet(Solver, par().solver);
    auto        &mat     = solver.getFMat();
    envGetTmp(FermionField, source);
    envGetTmp(FermionField, sol);
    envGetTmp(FermionField, tmp);
    
    LOG(Message) << "Using a2a vecs for prop '" << "' on source '"
    << par().source << "'" << std::endl;
    
    // loop over source colors
    for (unsigned int c = 0; c < FImpl::Dimension; ++c)
    {
        LOG(Message) << "Inversion for color= " << c << std::endl;
        // source conversion for 4D sources
        LOG(Message) << "Import source" << std::endl;
        
        if (!env().isObject5d(par().source))
        {
            if (Ls_ == 1)
            {
                PropToFerm<FImpl>(source, fullSrc, c);
            }
            else
            {
                PropToFerm<FImpl>(tmp, fullSrc, c);
                mat.ImportPhysicalFermionSource(tmp, source);
            }
        }
        // source conversion for 5D sources
        else
        {
            if (Ls_ != env().getObjectLs(par().source))
            {
                HADRONS_ERROR(Size, "Ls mismatch between quark action and source");
            }
            else
            {
                PropToFerm<FImpl>(source, fullSrc, c);
            }
        }
        
        sol = Zero();
        auto N = left.size();
        LOG(Message) << "Spectral Decomposition using " << N << " evecs" << std::endl;
        for (int i=0;i<N;i++) {
            const FermionField& tmp = left[i];
            const FermionField& tmp2 = right[i];
            axpy(sol,TensorRemove(innerProduct(tmp2,source)),tmp,sol);
        }
        
        LOG(Message) << "Export solution" << std::endl;
        FermToProp<FImpl>(prop, sol, c);
        std::cout<< "color " << c << " sol= " << sol << std::endl;
        // create 4D propagators from 5D one if necessary
        if (Ls_ > 1)
        {
            PropagatorField &p4d = envGet(PropagatorField, getName());
            mat.ExportPhysicalFermionSolution(sol, tmp);
            FermToProp<FImpl>(p4d, tmp, c);
        }
    }
}


END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MFermion_GaugeSpectProp_hpp_
