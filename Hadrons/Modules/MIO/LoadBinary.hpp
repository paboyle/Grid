/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MIO/LoadBinary.hpp

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
#ifndef Hadrons_MIO_LoadBinary_hpp_
#define Hadrons_MIO_LoadBinary_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                       Load a binary configurations                         *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MIO)

class LoadBinaryPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(LoadBinaryPar,
                                    std::string, file,
                                    std::string, format);
};

template <typename Impl>
class TLoadBinary: public Module<LoadBinaryPar>
{
public:
    typedef typename Impl::Field                  Field;
    typedef typename Impl::Simd                   Simd;
    typedef typename Field::vector_object         vobj;
    typedef typename vobj::scalar_object          sobj;
    typedef typename sobj::DoublePrecision        sobj_double;
    typedef BinarySimpleMunger<sobj_double, sobj> Munger;
public:
    // constructor
    TLoadBinary(const std::string name);
    // destructor
    virtual ~TLoadBinary(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(LoadBinary, TLoadBinary<GIMPL>, MIO);
MODULE_REGISTER_TMP(LoadBinaryScalarSU2, TLoadBinary<ScalarNxNAdjImplR<2>>, MIO);
MODULE_REGISTER_TMP(LoadBinaryScalarSU3, TLoadBinary<ScalarNxNAdjImplR<3>>, MIO);
MODULE_REGISTER_TMP(LoadBinaryScalarSU4, TLoadBinary<ScalarNxNAdjImplR<4>>, MIO);
MODULE_REGISTER_TMP(LoadBinaryScalarSU5, TLoadBinary<ScalarNxNAdjImplR<5>>, MIO);
MODULE_REGISTER_TMP(LoadBinaryScalarSU6, TLoadBinary<ScalarNxNAdjImplR<6>>, MIO);

/******************************************************************************
 *                         TLoadBinary implementation                         *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename Impl>
TLoadBinary<Impl>::TLoadBinary(const std::string name)
: Module<LoadBinaryPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename Impl>
std::vector<std::string> TLoadBinary<Impl>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename Impl>
std::vector<std::string> TLoadBinary<Impl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename Impl>
void TLoadBinary<Impl>::setup(void)
{
    envCreateLat(Field, getName());
}

// execution ///////////////////////////////////////////////////////////////////
template <typename Impl>
void TLoadBinary<Impl>::execute(void)
{
    Munger      munge;
    uint32_t    nersc_csum, scidac_csuma, scidac_csumb;
    auto        &U = envGet(Field, getName());
    std::string filename = par().file + "."
                           + std::to_string(vm().getTrajectory());

    LOG(Message) << "Loading " << par().format 
                 << " binary configuration from file '" << filename
                 << "'" << std::endl;
    BinaryIO::readLatticeObject<vobj, sobj_double>(U, filename, munge, 0, 
                                                   par().format, nersc_csum,
                                                   scidac_csuma, scidac_csumb);
    LOG(Message) << "Checksums:" << std::endl;
    LOG(Message) << "  NERSC    " << nersc_csum << std::endl;
    LOG(Message) << "  SciDAC A " << scidac_csuma << std::endl;
    LOG(Message) << "  SciDAC B " << scidac_csumb << std::endl;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MIO_LoadBinary_hpp_
