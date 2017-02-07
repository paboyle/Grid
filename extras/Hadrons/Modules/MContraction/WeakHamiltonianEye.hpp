/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: extras/Hadrons/Modules/MContraction/WeakHamiltonianEye.hpp

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

#ifndef Hadrons_WeakHamiltonianEye_hpp_
#define Hadrons_WeakHamiltonianEye_hpp_

#include <Grid/Hadrons/Modules/MContraction/WeakHamiltonian.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         WeakHamiltonianEye                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

enum
{
    S_diag = 0,
    E_diag = 1,
    n_eye_diag = 2
};

// Saucer and Eye subdiagram contractions.
#define MAKE_SE_BODY(Q_1, Q_2, Q_3, gamma) (Q_3*g5*Q_1*adj(Q_2)*(g5*gamma))
#define MAKE_SE_LOOP(Q_loop, gamma) (Q_loop*gamma)

class TWeakHamiltonianEye: public Module<WeakHamiltonianPar>
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
    TWeakHamiltonianEye(const std::string name);
    // destructor
    virtual ~TWeakHamiltonianEye(void) = default;
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_NS(WeakHamiltonianEye, TWeakHamiltonianEye, MContraction);

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_WeakHamiltonianEye_hpp_
