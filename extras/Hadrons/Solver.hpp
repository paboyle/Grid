/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: extras/Hadrons/EigenPack.hpp

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
#ifndef Hadrons_Solver_hpp_
#define Hadrons_Solver_hpp_

#include <Grid/Hadrons/Global.hpp>

BEGIN_HADRONS_NAMESPACE

template <typename FImpl>
class Solver
{
public:
    typedef typename FImpl::FermionField                      FermionField;
    typedef FermionOperator<FImpl>                            FMat; 
    typedef std::function<void(FermionField &, 
                               const FermionField &)>         SolverFn;
public:
    Solver(SolverFn fn, FMat &mat): mat_(mat), fn_(fn) {}

    void operator()(FermionField &sol, const FermionField &src)
    {
        fn_(sol, src);
    }

    FMat & getFMat(void)
    {
        return mat_;
    }
private:
    FMat     &mat_;
    SolverFn fn_;
};

END_HADRONS_NAMESPACE

#endif // Hadrons_Solver_hpp_
