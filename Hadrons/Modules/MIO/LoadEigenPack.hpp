/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MIO/LoadEigenPack.hpp

Copyright (C) 2015-2019

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
#ifndef Hadrons_MIO_LoadEigenPack_hpp_
#define Hadrons_MIO_LoadEigenPack_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/EigenPack.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                   Load eigen vectors/values package                        *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MIO)

class LoadEigenPackPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(LoadEigenPackPar,
                                    std::string, filestem,
                                    bool, multiFile,
                                    unsigned int, size,
                                    unsigned int, Ls,
                                    std::string, gaugeXform);
};

template <typename Pack, typename GImpl>
class TLoadEigenPack: public Module<LoadEigenPackPar>
{
public:
    typedef typename Pack::Field   Field;
    typedef typename Pack::FieldIo FieldIo;
    typedef BaseEigenPack<Field>   BasePack;

public:
    GAUGE_TYPE_ALIASES(GImpl, );
    typedef typename GImpl::GaugeLinkField GaugeMat;
public:
    // constructor
    TLoadEigenPack(const std::string name);
    // destructor
    virtual ~TLoadEigenPack(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(LoadFermionEigenPack, ARG(TLoadEigenPack<FermionEigenPack<FIMPL>, GIMPL>), MIO);
#ifdef GRID_DEFAULT_PRECISION_DOUBLE
MODULE_REGISTER_TMP(LoadFermionEigenPackIo32, ARG(TLoadEigenPack<FermionEigenPack<FIMPL, FIMPLF>, GIMPL>), MIO);
#endif

/******************************************************************************
 *                    TLoadEigenPack implementation                           *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename Pack, typename GImpl>
TLoadEigenPack<Pack, GImpl>::TLoadEigenPack(const std::string name)
: Module<LoadEigenPackPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename Pack, typename GImpl>
std::vector<std::string> TLoadEigenPack<Pack, GImpl>::getInput(void)
{
    std::vector<std::string> in;

    if (!par().gaugeXform.empty())
    {
        in = {par().gaugeXform};
    }
    
    return in;
}

template <typename Pack, typename GImpl>
std::vector<std::string> TLoadEigenPack<Pack, GImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename Pack, typename GImpl>
void TLoadEigenPack<Pack, GImpl>::setup(void)
{
    GridBase *gridIo = nullptr;

    if (typeHash<Field>() != typeHash<FieldIo>())
    {
        gridIo = envGetRbGrid(FieldIo, par().Ls);
    }
    envCreateDerived(BasePack, Pack, getName(), par().Ls, par().size, 
                     envGetRbGrid(Field, par().Ls), gridIo);

    if (!par().gaugeXform.empty())
    {
        if (par().Ls > 1)
        {
            LOG(Message) << "Setup 5d GaugeMat for Ls = " << par().Ls << std::endl;
            envTmp(GaugeMat,    "tmpXform", par().Ls, envGetGrid5(Field, par().Ls));
            envTmp(GaugeMat, "tmpXformOdd", par().Ls, envGetRbGrid5(Field, par().Ls));
        }
        else
        {
            LOG(Message) << "Setup 4d GaugeMat for Ls = " << par().Ls << std::endl;
            envTmp(GaugeMat,    "tmpXform", par().Ls, envGetGrid(Field));
            envTmp(GaugeMat, "tmpXformOdd", par().Ls, envGetRbGrid(Field));
        }
        
    }
}

// execution ///////////////////////////////////////////////////////////////////
template <typename Pack, typename GImpl>
void TLoadEigenPack<Pack, GImpl>::execute(void)
{
    auto &epack = envGetDerived(BasePack, Pack, getName());

    epack.read(par().filestem, par().multiFile, vm().getTrajectory());
    epack.eval.resize(par().size);

    if (!par().gaugeXform.empty())
    {

        LOG(Message) << "Applying gauge transformation to eigenvectors " << getName()
                     << " using " << par().gaugeXform << std::endl;
        auto &xform = envGet(GaugeMat, par().gaugeXform);
        envGetTmp(GaugeMat,    tmpXform);
        envGetTmp(GaugeMat, tmpXformOdd);

        if (par().Ls > 1) 
        {
            LOG(Message) << "Creating 5d GaugeMat from " << par().gaugeXform << std::endl;
            startTimer("5-d gauge transform creation");
            for (unsigned int j = 0; j < par().Ls; j++)
            {
                InsertSlice(xform, tmpXform, j, 0);
            }
            stopTimer("5-d gauge transform creation");
        }

        pickCheckerboard(Odd, tmpXformOdd, tmpXform);
        startTimer("Transform application");
        for (unsigned int i = 0; i < par().size; i++)
        {
            LOG(Message) << "Applying gauge transformation to eigenvector i = " << i << "/" << par().size << std::endl;
            epack.evec[i].checkerboard = Odd;
            epack.evec[i] = tmpXformOdd * epack.evec[i];
        }
        stopTimer("Transform application");
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MIO_LoadEigenPack_hpp_
