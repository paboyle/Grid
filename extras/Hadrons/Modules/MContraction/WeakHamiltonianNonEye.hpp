/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: extras/Hadrons/Modules/MContraction/WeakHamiltonianNonEye.hpp

Copyright (C) 2017

Author: Andrew Lawson    <andrew.lawson1991@gmail.com>

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

#ifndef Hadrons_WeakHamiltonianNonEye_hpp_
#define Hadrons_WeakHamiltonianNonEye_hpp_

#include <Grid/Hadrons/Modules/MContraction/WeakHamiltonian.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         WeakHamiltonianNonEye                              *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

enum
{
    W_diag = 0,
    C_diag = 1,
    n_noneye_diag = 2
};

// Wing and Connected subdiagram contractions
#define MAKE_CW_SUBDIAG(Q_1, Q_2, gamma) (Q_1*adj(Q_2)*(g5*gamma))

class TWeakHamiltonianNonEye: public Module<WeakHamiltonianPar>
{
public:
    TYPE_ALIASES(FIMPL,)
    class Result: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                                        std::string, name,
                                        std::vector<Complex>, corr);
    };
public:
    // constructor
    TWeakHamiltonianNonEye(const std::string name);
    // destructor
    virtual ~TWeakHamiltonianNonEye(void) = default;
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_NS(WeakHamiltonianNonEye, TWeakHamiltonianNonEye, MContraction);

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_WeakHamiltonianNonEye_hpp_
