/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Archive/Modules/WeakHamiltonian.hpp

Copyright (C) 2015-2019

Author: Antonin Portelli <antonin.portelli@me.com>
Author: Lanny91 <andrew.lawson@gmail.com>

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

#ifndef Hadrons_MContraction_WeakHamiltonian_hpp_
#define Hadrons_MContraction_WeakHamiltonian_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         WeakHamiltonian                                    *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

/*******************************************************************************
 * Utilities for contractions involving the Weak Hamiltonian.
 ******************************************************************************/
//// Sum and store correlator.
#define MAKE_DIAG(exp, buf, res, n)\
sliceSum(exp, buf, Tp);\
res.name = (n);\
res.corr.resize(buf.size());\
for (unsigned int t = 0; t < buf.size(); ++t)\
{\
    res.corr[t] = TensorRemove(buf[t]);\
}

//// Contraction of mu index: use 'mu' variable in exp.
#define SUM_MU(buf,exp)\
buf = zero;\
for (unsigned int mu = 0; mu < ndim; ++mu)\
{\
    buf += exp;\
}

enum 
{
  i_V = 0,
  i_A = 1,
  n_i = 2
};

class WeakHamiltonianPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(WeakHamiltonianPar,
                                    std::string, q1,
                                    std::string, q2,
                                    std::string, q3,
                                    std::string, q4,
                                    unsigned int, tSnk,
                                    std::string, output);
};

#define MAKE_WEAK_MODULE(modname)\
class T##modname: public Module<WeakHamiltonianPar>\
{\
public:\
    FERM_TYPE_ALIASES(FIMPL,)\
    class Result: Serializable\
    {\
    public:\
        GRID_SERIALIZABLE_CLASS_MEMBERS(Result,\
                                        std::string, name,\
                                        std::vector<Complex>, corr);\
    };\
public:\
    /* constructor */ \
    T##modname(const std::string name);\
    /* destructor */ \
    virtual ~T##modname(void) {};\
    /* dependency relation */ \
    virtual std::vector<std::string> getInput(void);\
    virtual std::vector<std::string> getOutput(void);\
public:\
    std::vector<std::string> VA_label = {"V", "A"};\
protected:\
    /* setup */ \
    virtual void setup(void);\
    /* execution */ \
    virtual void execute(void);\
};\
MODULE_REGISTER(modname, T##modname, MContraction);

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_WeakHamiltonian_hpp_
