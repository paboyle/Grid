/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: extras/Hadrons/A2AVectors.hpp

Copyright (C) 2015-2018

Author: Antonin Portelli <antonin.portelli@me.com>
Author: fionnoh <fionnoh@gmail.com>

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
#ifndef A2A_Vectors_hpp_
#define A2A_Vectors_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Environment.hpp>
#include <Grid/Hadrons/Solver.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *               Classes to generate V & W all-to-all vectors                 *
 ******************************************************************************/
template <typename FImpl>
class A2AVectorsSchurDiagTwo
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    SOLVER_TYPE_ALIASES(FImpl,);
public:
    A2AVectorsSchurDiagTwo(FMat &action, Solver &solver);
    virtual ~A2AVectorsSchurDiagTwo(void) = default;
    void makeLowModeV(FermionField &vout, 
                      const FermionField &evec, const Real &eval);
    void makeLowModeV5D(FermionField &vout_4d, FermionField &vout_5d, 
                        const FermionField &evec, const Real &eval);
    void makeLowModeW(FermionField &wout, 
                      const FermionField &evec, const Real &eval);
    void makeLowModeW5D(FermionField &wout_4d, FermionField &wout_5d, 
                        const FermionField &evec, const Real &eval);
    void makeHighModeV(FermionField &vout, const FermionField &noise);
    void makeHighModeV5D(FermionField &vout_4d, FermionField &vout_5d, 
                         const FermionField &noise_5d);
    void makeHighModeW(FermionField &wout, const FermionField &noise);
    void makeHighModeW5D(FermionField &vout_5d, FermionField &wout_5d, 
                         const FermionField &noise_5d);
private:
    FMat                                     &action_;
    Solver                                   &solver_;
    GridBase                                 *fGrid_, *frbGrid_, *gGrid_;
    bool                                     is5d_;
    FermionField                             src_o_, sol_e_, sol_o_, tmp_, tmp5_;
    SchurDiagTwoOperator<FMat, FermionField> op_;
};

/******************************************************************************
 *               A2AVectorsSchurDiagTwo template implementation               *
 ******************************************************************************/
template <typename FImpl>
A2AVectorsSchurDiagTwo<FImpl>::A2AVectorsSchurDiagTwo(FMat &action, Solver &solver)
: action_(action)
, solver_(solver)
, fGrid_(action_.FermionGrid())
, frbGrid_(action_.FermionRedBlackGrid())
, gGrid_(action_.GaugeGrid())
, src_o_(frbGrid_)
, sol_e_(frbGrid_)
, sol_o_(frbGrid_)
, tmp_(frbGrid_)
, tmp5_(fGrid_)
, op_(action_)
{}

template <typename FImpl>
void A2AVectorsSchurDiagTwo<FImpl>::makeLowModeV(FermionField &vout, const FermionField &evec, const Real &eval)
{
    src_o_ = evec;
    src_o_.checkerboard = Odd;
    pickCheckerboard(Even, sol_e_, vout);
    pickCheckerboard(Odd, sol_o_, vout);

    /////////////////////////////////////////////////////
    // v_ie = -(1/eval_i) * MeeInv Meo MooInv evec_i
    /////////////////////////////////////////////////////
    action_.MooeeInv(src_o_, tmp_);
    assert(tmp_.checkerboard == Odd);
    action_.Meooe(tmp_, sol_e_);
    assert(sol_e_.checkerboard == Even);
    action_.MooeeInv(sol_e_, tmp_);
    assert(tmp_.checkerboard == Even);
    sol_e_ = (-1.0 / eval) * tmp_;
    assert(sol_e_.checkerboard == Even);

    /////////////////////////////////////////////////////
    // v_io = (1/eval_i) * MooInv evec_i
    /////////////////////////////////////////////////////
    action_.MooeeInv(src_o_, tmp_);
    assert(tmp_.checkerboard == Odd);
    sol_o_ = (1.0 / eval) * tmp_;
    assert(sol_o_.checkerboard == Odd);
    setCheckerboard(vout, sol_e_);
    assert(sol_e_.checkerboard == Even);
    setCheckerboard(vout, sol_o_);
    assert(sol_o_.checkerboard == Odd);
}

template <typename FImpl>
void A2AVectorsSchurDiagTwo<FImpl>::makeLowModeV5D(FermionField &vout_4d, FermionField &vout_5d, const FermionField &evec, const Real &eval)
{
    makeLowModeV(vout_5d, evec, eval);
    action_.ExportPhysicalFermionSolution(vout_5d, vout_4d);
}

template <typename FImpl>
void A2AVectorsSchurDiagTwo<FImpl>::makeLowModeW(FermionField &wout, const FermionField &evec, const Real &eval)
{
    src_o_ = evec;
    src_o_.checkerboard = Odd;
    pickCheckerboard(Even, sol_e_, wout);
    pickCheckerboard(Odd, sol_o_, wout);

    /////////////////////////////////////////////////////
    // w_ie = - MeeInvDag MoeDag Doo evec_i
    /////////////////////////////////////////////////////
    op_.Mpc(src_o_, tmp_);
    assert(tmp_.checkerboard == Odd);
    action_.MeooeDag(tmp_, sol_e_);
    assert(sol_e_.checkerboard == Even);
    action_.MooeeInvDag(sol_e_, tmp_);
    assert(tmp_.checkerboard == Even);
    sol_e_ = (-1.0) * tmp_;

    /////////////////////////////////////////////////////
    // w_io = Doo evec_i
    /////////////////////////////////////////////////////
    op_.Mpc(src_o_, sol_o_);
    assert(sol_o_.checkerboard == Odd);
    setCheckerboard(wout, sol_e_);
    assert(sol_e_.checkerboard == Even);
    setCheckerboard(wout, sol_o_);
    assert(sol_o_.checkerboard == Odd);
}

template <typename FImpl>
void A2AVectorsSchurDiagTwo<FImpl>::makeLowModeW5D(FermionField &wout_4d, 
                                                   FermionField &wout_5d, 
                                                   const FermionField &evec, 
                                                   const Real &eval)
{
    makeLowModeW(tmp5_, evec, eval);
    action_.DminusDag(tmp5_, wout_5d);
    action_.ExportPhysicalFermionSource(wout_5d, wout_4d);
}

template <typename FImpl>
void A2AVectorsSchurDiagTwo<FImpl>::makeHighModeV(FermionField &vout, 
                                                  const FermionField &noise)
{
    solver_(vout, noise);
}

template <typename FImpl>
void A2AVectorsSchurDiagTwo<FImpl>::makeHighModeV5D(FermionField &vout_4d, 
                                                    FermionField &vout_5d, 
                                                    const FermionField &noise)
{
    if (noise._grid->Dimensions() == fGrid_->Dimensions() - 1)
    {
        action_.ImportPhysicalFermionSource(noise, tmp5_);
    }
    else
    {
        tmp5_ = noise;
    }
    makeHighModeV(vout_5d, tmp5_);
    action_.ExportPhysicalFermionSolution(vout_5d, vout_4d);
}

template <typename FImpl>
void A2AVectorsSchurDiagTwo<FImpl>::makeHighModeW(FermionField &wout, 
                                                  const FermionField &noise)
{
    wout = noise;
}

template <typename FImpl>
void A2AVectorsSchurDiagTwo<FImpl>::makeHighModeW5D(FermionField &wout_4d, 
                                                    FermionField &wout_5d, 
                                                    const FermionField &noise)
{
    if (noise._grid->Dimensions() == fGrid_->Dimensions() - 1)
    {
        action_.ImportUnphysicalFermion(noise, wout_5d);
        wout_4d = noise;
    }
    else
    {
        wout_5d = noise;
        action_.ExportPhysicalFermionSource(wout_5d, wout_4d);
    }
}

END_HADRONS_NAMESPACE

#endif // A2A_Vectors_hpp_