#ifndef A2A_Vectors_hpp_
#define A2A_Vectors_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Environment.hpp>
#include <Grid/Hadrons/Solver.hpp>

BEGIN_HADRONS_NAMESPACE

////////////////////////////////
// A2A Modes
////////////////////////////////

template <typename FImpl>
class A2AVectorsSchurDiagTwo
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    SOLVER_TYPE_ALIASES(FImpl,);
public:
    A2AVectorsSchurDiagTwo(FMat &action, Solver &solver);
    virtual ~A2AVectorsSchurDiagTwo(void) = default;
    void makeLowModeV(FermionField &vout, const FermionField &evec, const Real &eval);
    void makeLowModeV5D(FermionField &vout_4d, FermionField &vout_5d, const FermionField &evec, const Real &eval);
    void makeLowModeW(FermionField &wout, const FermionField &evec, const Real &eval);
    void makeLowModeW5D(FermionField &wout_4d, FermionField &wout_5d, const FermionField &evec, const Real &eval);
    void makeHighModeV(FermionField &vout, const FermionField &noise);
    void makeHighModeV5D(FermionField &vout_4d, FermionField &vout_5d, const FermionField &noise_5d);
    void makeHighModeW(FermionField &wout, const FermionField &noise);
    void makeHighModeW5D(FermionField &vout_5d, FermionField &wout_5d, const FermionField &noise_5d);
private:
    FMat                                     &action_;
    Solver                                   &solver_;
    GridBase                                 *fGrid_, *frbGrid_, *gGrid_;
    bool                                     is5d_;
    FermionField                             src_o_, sol_e_, sol_o_, tmp_, tmp5_;
    SchurDiagTwoOperator<FMat, FermionField> op_;
};

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


    // A2AVectorsSchurDiagTwo(const int Nl, const int Nh,
    //                        std::vector<FermionField> &v,
    //                        std::vector<FermionField> &w,
    //                         const bool _return_5d)
    //                     : Nl(_Nl), Nh(_Nh),
    //                         return_5d(_return_5d)
    // {
    //     if (!return_5d)
    //     {
    //         init_resize(1, Nl + Nh);
    //     }
    //     else
    //     {
    //         init_resize(Nl + Nh, Nl + Nh);
    //     }   
    // }

    // void init_resize(const size_t size_5d, const size_t size_4d,
    //                     GridBase *grid_5d, GridBase *grid_4d)
    // {
    //     w_5d.resize(size_5d, grid_5d);
    //     v_5d.resize(size_5d, grid_5d);
    //     w_4d.resize(size_4d, grid_4d);
    //     v_4d.resize(size_4d, grid_4d);
    // }

    // int get_Nh(void) const
    // {
    //     return Nh;
    // }

    // int get_Nl(void) const
    // {
    //     return Nl;
    // }

    // void low_modes(int il, const Field &evec, const Real &eval, Matrix &action)
    // {
    //     int i5d;

    //     i5d = 0;
    //     if (return_5d) i5d = il;
    //     this->low_mode_v(v_5d[i5d], v_4d[il], evec, eval, action);
    //     this->low_mode_w(w_5d[i5d], w_4d[il], evec, eval, action);
    // }

    // void high_modes(int ih, Field &source_5d, Field &w_source_5d, 
    //                 Field &source_4d, Solver &solver)
    // {
    //     int i5d;

    //     i5d = 0;
    //     if (return_5d) i5d = ih + Nl;
    //     this->high_mode_v(source_5d, v_5d[i5d], v_4d[ih + Nl], solver);
    //     this->high_mode_w(w_source_5d, source_4d, w_5d[i5d], w_4d[ih + Nl]);
    // }

    // void return_v(int i, Field &vout_5d, Field &vout_4d)
    // {
    //     vout_4d = v_4d[i];
    //     if (!(return_5d)) i = 0;
    //     vout_5d = v_5d[i];
    // }

    // void return_w(int i, Field &wout_5d, Field &wout_4d)
    // {
    //     wout_4d = w_4d[i];
    //     if (!(return_5d)) i = 0;
    //     wout_5d = w_5d[i];
    // }

    // void low_mode_v(Field &vout_5d, Field &vout_4d, const Field &evec,
    //                 const Real &eval, Matrix &action)
    // {
    //     GridBase *grid = action.RedBlackGrid();
    //     Field src_o(grid);
    //     Field sol_e(grid);
    //     Field sol_o(grid);
    //     Field tmp(grid);

    //     src_o = evec;
    //     src_o.checkerboard = Odd;
    //     pickCheckerboard(Even, sol_e, vout_5d);
    //     pickCheckerboard(Odd, sol_o, vout_5d);

    //     /////////////////////////////////////////////////////
    //     // v_ie = -(1/eval_i) * MeeInv Meo MooInv evec_i
    //     /////////////////////////////////////////////////////
    //     action.MooeeInv(src_o, tmp);
    //     assert(tmp.checkerboard == Odd);
    //     action.Meooe(tmp, sol_e);
    //     assert(sol_e.checkerboard == Even);
    //     action.MooeeInv(sol_e, tmp);
    //     assert(tmp.checkerboard == Even);
    //     sol_e = (-1.0 / eval) * tmp;
    //     assert(sol_e.checkerboard == Even);

    //     /////////////////////////////////////////////////////
    //     // v_io = (1/eval_i) * MooInv evec_i
    //     /////////////////////////////////////////////////////
    //     action.MooeeInv(src_o, tmp);
    //     assert(tmp.checkerboard == Odd);
    //     sol_o = (1.0 / eval) * tmp;
    //     assert(sol_o.checkerboard == Odd);

    //     setCheckerboard(vout_5d, sol_e);
    //     assert(sol_e.checkerboard == Even);
    //     setCheckerboard(vout_5d, sol_o);
    //     assert(sol_o.checkerboard == Odd);

    //     action.ExportPhysicalFermionSolution(vout_5d, vout_4d);
    // }

    // void low_mode_w(Field &wout_5d, Field &wout_4d, const Field &evec,
    //                 const Real &eval, Matrix &action)
    // {
    //     GridBase *grid = action.RedBlackGrid();
    //     SchurDiagTwoOperator<Matrix, Field> _HermOpEO(action);

    //     Field src_o(grid);
    //     Field sol_e(grid);
    //     Field sol_o(grid);
    //     Field tmp(grid);

    //     GridBase *fgrid = action.Grid();
    //     Field tmp_wout(fgrid);

    //     src_o = evec;
    //     src_o.checkerboard = Odd;
    //     pickCheckerboard(Even, sol_e, tmp_wout);
    //     pickCheckerboard(Odd, sol_o, tmp_wout);

    //     /////////////////////////////////////////////////////
    //     // w_ie = - MeeInvDag MoeDag Doo evec_i
    //     /////////////////////////////////////////////////////
    //     _HermOpEO.Mpc(src_o, tmp);
    //     assert(tmp.checkerboard == Odd);
    //     action.MeooeDag(tmp, sol_e);
    //     assert(sol_e.checkerboard == Even);
    //     action.MooeeInvDag(sol_e, tmp);
    //     assert(tmp.checkerboard == Even);
    //     sol_e = (-1.0) * tmp;

    //     /////////////////////////////////////////////////////
    //     // w_io = Doo evec_i
    //     /////////////////////////////////////////////////////
    //     _HermOpEO.Mpc(src_o, sol_o);
    //     assert(sol_o.checkerboard == Odd);

    //     setCheckerboard(tmp_wout, sol_e);
    //     assert(sol_e.checkerboard == Even);
    //     setCheckerboard(tmp_wout, sol_o);
    //     assert(sol_o.checkerboard == Odd);

    //     action.DminusDag(tmp_wout, wout_5d);

    //     action.ExportPhysicalFermionSource(wout_5d, wout_4d);
    // }

    // void high_mode_v(const Field &source, Field &vout_5d, Field &vout_4d,
    //                     Matrix &action, Solver &solver)
    // {
    //     GridBase *fgrid = action.Grid();
    //     solver(vout_5d, source); // Note: solver is solver(out, in)
    //     action.ExportPhysicalFermionSolution(vout_5d, vout_4d);
    // }

    // void high_mode_w(const Field &w_source_5d, const Field &source_4d, 
    //                     Field &wout_5d, Field &wout_4d)
    // {
    //     wout_5d = w_source_5d;
    //     wout_4d = source_4d;
    // }

// TODO: A2A for coarse eigenvectors

// template <class FineField, class CoarseField, class Matrix, class Solver>
// class A2ALMSchurDiagTwoCoarse : public A2AModesSchurDiagTwo<FineField, Matrix, Solver>
// {
//   private:
//     const std::vector<FineField> &subspace;
//     const std::vector<CoarseField> &evec_coarse;
//     const std::vector<RealD> &eval_coarse;
//     Matrix &action;

//   public:
//     A2ALMSchurDiagTwoCoarse(const std::vector<FineField> &_subspace, const std::vector<CoarseField> &_evec_coarse, const std::vector<RealD> &_eval_coarse, Matrix &_action)
//         : subspace(_subspace), evec_coarse(_evec_coarse), eval_coarse(_eval_coarse), action(_action){};

//     void operator()(int i, FineField &vout, FineField &wout)
//     {
//         FineField prom_evec(subspace[0]._grid);
//         blockPromote(evec_coarse[i], prom_evec, subspace);
//         this->low_mode_v(action, prom_evec, eval_coarse[i], vout);
//         this->low_mode_w(action, prom_evec, eval_coarse[i], wout);
//     }
// };

END_HADRONS_NAMESPACE

#endif // A2A_Vectors_hpp_