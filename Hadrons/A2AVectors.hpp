/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/A2AVectors.hpp

Copyright (C) 2015-2019

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

#include <Hadrons/Global.hpp>
#include <Hadrons/Environment.hpp>
#include <Hadrons/Solver.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                 Class to generate V & W all-to-all vectors                 *
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

template <typename FImpl>
class A2AVectorsSchurStaggered
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    SOLVER_TYPE_ALIASES(FImpl,);
public:
    A2AVectorsSchurStaggered(FMat &action, Solver &solver);
    virtual ~A2AVectorsSchurStaggered(void) = default;
    void makeLowModeV(FermionField &vout,
                      const FermionField &evec, const ComplexD eval, const int sign=0);
    void makeLowModeV5D(FermionField &vout_4d, FermionField &vout_5d,
                        const FermionField &evec, const ComplexD eval, const int sign=0);
    void makeLowModeW(FermionField &wout,
                      const FermionField &evec, const ComplexD eval, const int sign=0);
    void makeLowModeW5D(FermionField &wout_4d, FermionField &wout_5d,
                        const FermionField &evec, const ComplexD eval, const int sign=0);
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
    SchurStaggeredOperator<FMat, FermionField> op_;
};

/******************************************************************************
 *                  Methods for V & W all-to-all vectors I/O                  *
 ******************************************************************************/
class A2AVectorsIo
{
public:
    struct Record: Serializable
    {
        GRID_SERIALIZABLE_CLASS_MEMBERS(Record,
                                        unsigned int, index);
        Record(void): index(0) {}
    };
public:
    template <typename Field>
    static void write(const std::string fileStem, std::vector<Field> &vec, 
                      const bool multiFile, const int trajectory = -1);
    template <typename Field>
    static void read(std::vector<Field> &vec, const std::string fileStem,
                     const bool multiFile, const int trajectory = -1);
private:
    static inline std::string vecFilename(const std::string stem, const int traj, 
                                          const bool multiFile)
    {
        std::string t = (traj < 0) ? "" : ("." + std::to_string(traj));

        if (multiFile)
        {
            return stem + t;
        }
        else
        {
            return stem + t + ".bin";
        }
    }
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
    src_o_.Checkerboard() = Odd;
    pickCheckerboard(Even, sol_e_, vout);
    pickCheckerboard(Odd, sol_o_, vout);

    /////////////////////////////////////////////////////
    // v_ie = -(1/eval_i) * MeeInv Meo MooInv evec_i
    /////////////////////////////////////////////////////
    action_.MooeeInv(src_o_, tmp_);
    assert(tmp_.Checkerboard() == Odd);
    action_.Meooe(tmp_, sol_e_);
    assert(sol_e_.Checkerboard() == Even);
    action_.MooeeInv(sol_e_, tmp_);
    assert(tmp_.Checkerboard() == Even);
    sol_e_ = (-1.0 / eval) * tmp_;
    assert(sol_e_.Checkerboard() == Even);

    /////////////////////////////////////////////////////
    // v_io = (1/eval_i) * MooInv evec_i
    /////////////////////////////////////////////////////
    action_.MooeeInv(src_o_, tmp_);
    assert(tmp_.Checkerboard() == Odd);
    sol_o_ = (1.0 / eval) * tmp_;
    assert(sol_o_.Checkerboard() == Odd);
    setCheckerboard(vout, sol_e_);
    assert(sol_e_.Checkerboard() == Even);
    setCheckerboard(vout, sol_o_);
    assert(sol_o_.Checkerboard() == Odd);
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
    src_o_.Checkerboard() = Odd;
    pickCheckerboard(Even, sol_e_, wout);
    pickCheckerboard(Odd, sol_o_, wout);

    /////////////////////////////////////////////////////
    // w_ie = - MeeInvDag MoeDag Doo evec_i
    /////////////////////////////////////////////////////
    op_.Mpc(src_o_, tmp_);
    assert(tmp_.Checkerboard() == Odd);
    action_.MeooeDag(tmp_, sol_e_);
    assert(sol_e_.Checkerboard() == Even);
    action_.MooeeInvDag(sol_e_, tmp_);
    assert(tmp_.Checkerboard() == Even);
    sol_e_ = (-1.0) * tmp_;

    /////////////////////////////////////////////////////
    // w_io = Doo evec_i
    /////////////////////////////////////////////////////
    op_.Mpc(src_o_, sol_o_);
    assert(sol_o_.Checkerboard() == Odd);
    setCheckerboard(wout, sol_e_);
    assert(sol_e_.Checkerboard() == Even);
    setCheckerboard(wout, sol_o_);
    assert(sol_o_.Checkerboard() == Odd);
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
    if (noise.Grid()->Dimensions() == fGrid_->Dimensions() - 1)
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
    if (noise.Grid()->Dimensions() == fGrid_->Dimensions() - 1)
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

/******************************************************************************
 *               A2AVectorsSchurStaggered template implementation               *
 ******************************************************************************/
template <typename FImpl>
A2AVectorsSchurStaggered<FImpl>::A2AVectorsSchurStaggered(FMat &action, Solver &solver)
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
void A2AVectorsSchurStaggered<FImpl>::makeLowModeV(FermionField &vout,
                                                   const FermionField &evec,
                                                   const ComplexD eval,
                                                   const int sign)
{
    ComplexD eval_ = eval;
    // evec_o = -evec_o ?
    if(sign){eval_=conjugate(eval);}
    src_o_ = evec;
    src_o_.Checkerboard() = Odd;
    pickCheckerboard(Even, sol_e_, vout);
    pickCheckerboard(Odd, sol_o_, vout);
    
    /////////////////////////////////////////////////////
    /// v_e = (1/eval^(*)) * (-i/Im(eval) * Meo evec_o)
    /////////////////////////////////////////////////////
    action_.Meooe(src_o_, tmp_);
    ComplexD minusI(0, -1.0);
    ComplexD cc = minusI/eval.imag()/eval_;
    sol_e_ = cc * tmp_;
    
    /////////////////////////////////////////////////////
    /// v_o = (1/eval^(*)) * evec_o
    /////////////////////////////////////////////////////
    cc = 1.0/eval_;
    sol_o_ = cc * src_o_;
    if(sign){sol_o_ = -sol_o_;}
    
    setCheckerboard(vout, sol_e_);
    assert(sol_e_.Checkerboard() == Even);
    setCheckerboard(vout, sol_o_);
    assert(sol_o_.Checkerboard() == Odd);
}

template <typename FImpl>
void A2AVectorsSchurStaggered<FImpl>::makeLowModeV5D(FermionField &vout_4d,
                                                     FermionField &vout_5d,
                                                     const FermionField &evec,
                                                     const ComplexD eval,
                                                     const int sign)
{
    makeLowModeV(vout_5d, evec, eval, sign);
    action_.ExportPhysicalFermionSolution(vout_5d, vout_4d);
}

template <typename FImpl>
void A2AVectorsSchurStaggered<FImpl>::makeLowModeW(FermionField &wout,
                                                   const FermionField &evec,
                                                   const ComplexD eval,
                                                   const int sign)
{
    src_o_ = evec;
    src_o_.Checkerboard() = Odd;
    pickCheckerboard(Even, sol_e_, wout);
    pickCheckerboard(Odd, sol_o_, wout);
    
    /////////////////////////////////////////////////////
    /// v_e = (-i/eval * Meo evec_o)
    /////////////////////////////////////////////////////
    action_.Meooe(src_o_, tmp_);
    ComplexD minusI(0, -1.0);
    ComplexD cc = minusI/eval.imag();
    sol_e_ = cc * tmp_;
    
    /////////////////////////////////////////////////////
    /// v_o = evec_o
    /////////////////////////////////////////////////////
    sol_o_ = src_o_;
    if(sign){sol_o_ = -sol_o_;}
    
    setCheckerboard(wout, sol_e_);
    assert(sol_e_.Checkerboard() == Even);
    setCheckerboard(wout, sol_o_);
    assert(sol_o_.Checkerboard() == Odd);
}

template <typename FImpl>
void A2AVectorsSchurStaggered<FImpl>::makeLowModeW5D(FermionField &wout_4d,
                                                     FermionField &wout_5d,
                                                     const FermionField &evec,
                                                     const ComplexD eval,
                                                     const int sign)
{
    makeLowModeW(tmp5_, evec, eval, sign);
    action_.DminusDag(tmp5_, wout_5d);
    action_.ExportPhysicalFermionSource(wout_5d, wout_4d);
}

template <typename FImpl>
void A2AVectorsSchurStaggered<FImpl>::makeHighModeV(FermionField &vout,
                                                  const FermionField &noise)
{
    solver_(vout, noise);
}

template <typename FImpl>
void A2AVectorsSchurStaggered<FImpl>::makeHighModeV5D(FermionField &vout_4d,
                                                      FermionField &vout_5d,
                                                      const FermionField &noise)
{
    if (noise.Grid()->Dimensions() == fGrid_->Dimensions() - 1)
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
void A2AVectorsSchurStaggered<FImpl>::makeHighModeW(FermionField &wout,
                                                  const FermionField &noise)
{
    wout = noise;
}

template <typename FImpl>
void A2AVectorsSchurStaggered<FImpl>::makeHighModeW5D(FermionField &wout_4d,
                                                    FermionField &wout_5d,
                                                    const FermionField &noise)
{
    if (noise.Grid()->Dimensions() == fGrid_->Dimensions() - 1)
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






/******************************************************************************
 *               all-to-all vectors I/O template implementation               *
 ******************************************************************************/
template <typename Field>
void A2AVectorsIo::write(const std::string fileStem, std::vector<Field> &vec, 
                         const bool multiFile, const int trajectory)
{
    Record       record;
    GridBase     *grid = vec[0].Grid();
    ScidacWriter binWriter(grid->IsBoss());
    std::string  filename = vecFilename(fileStem, trajectory, multiFile);

    if (multiFile)
    {
        std::string fullFilename;

        for (unsigned int i = 0; i < vec.size(); ++i)
        {
            fullFilename = filename + "/elem" + std::to_string(i) + ".bin";

            LOG(Message) << "Writing vector " << i << std::endl;
            makeFileDir(fullFilename, grid);
            binWriter.open(fullFilename);
            record.index = i;
            binWriter.writeScidacFieldRecord(vec[i], record);
            binWriter.close();
        }
    }
    else
    {
        makeFileDir(filename, grid);
        binWriter.open(filename);
        for (unsigned int i = 0; i < vec.size(); ++i)
        {
            LOG(Message) << "Writing vector " << i << std::endl;
            record.index = i;
            binWriter.writeScidacFieldRecord(vec[i], record);
        }
        binWriter.close();
    }
}

template <typename Field>
void A2AVectorsIo::read(std::vector<Field> &vec, const std::string fileStem, 
                        const bool multiFile, const int trajectory)
{
    Record       record;
    ScidacReader binReader;
    std::string  filename = vecFilename(fileStem, trajectory, multiFile);

    if (multiFile)
    {
        std::string fullFilename;

        for (unsigned int i = 0; i < vec.size(); ++i)
        {
            fullFilename = filename + "/elem" + std::to_string(i) + ".bin";

            LOG(Message) << "Reading vector " << i << std::endl;
            binReader.open(fullFilename);
            binReader.readScidacFieldRecord(vec[i], record);
            binReader.close();
            if (record.index != i)
            {
                HADRONS_ERROR(Io, "vector index mismatch");
            }
        }
    }
    else
    {
        binReader.open(filename);
        for (unsigned int i = 0; i < vec.size(); ++i)
        {
            LOG(Message) << "Reading vector " << i << std::endl;
            binReader.readScidacFieldRecord(vec[i], record);
            if (record.index != i)
            {
                HADRONS_ERROR(Io, "vector index mismatch");
            }
        }
        binReader.close();
    }
}

END_HADRONS_NAMESPACE

#endif // A2A_Vectors_hpp_
