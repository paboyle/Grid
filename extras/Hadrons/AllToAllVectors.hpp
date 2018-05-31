#ifndef A2A_Vectors_hpp_
#define A2A_Vectors_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Environment.hpp>

BEGIN_HADRONS_NAMESPACE

////////////////////////////////
// A2A Modes
////////////////////////////////

template <class Field, class Matrix>
class A2AModesSchurDiagTwo
{
  public:
    A2AModesSchurDiagTwo(void) = default;
    virtual ~A2AModesSchurDiagTwo(void) = default;

    void low_mode_v(Matrix &_Matrix, const Field &evec, const RealD &eval, Field &vout)
    {

        GridBase *grid = _Matrix.RedBlackGrid();

        Field src_o(grid);
        Field sol_e(grid);
        Field sol_o(grid);
        Field tmp(grid);

        pickCheckerboard(Odd, src_o, evec);
        pickCheckerboard(Even, sol_e, vout);
        pickCheckerboard(Odd, sol_o, vout);

        /////////////////////////////////////////////////////
        // v_ie = -(1/eval_i)* Mee Meo MooInv evec_i
        /////////////////////////////////////////////////////
        _Matrix.MooeeInv(src_o, tmp);
        assert(tmp.checkerboard == Odd);
        _Matrix.Meooe(tmp, sol_e);
        assert(sol_e.checkerboard == Even);
        _Matrix.Mooee(sol_e, tmp);
        assert(tmp.checkerboard == Even);
        sol_e = -(1.0 / eval) * tmp;
        assert(sol_e.checkerboard == Even);

        /////////////////////////////////////////////////////
        // v_io = -(1/eval_i)* MooInv evec_i
        /////////////////////////////////////////////////////
        _Matrix.MooeeInv(src_o, tmp);
        assert(tmp.checkerboard == Odd);
        sol_o = -(1.0 / eval) * tmp;
        assert(sol_o.checkerboard == Odd);

        setCheckerboard(vout, sol_e);
        assert(sol_e.checkerboard == Even);
        setCheckerboard(vout, sol_o);
        assert(sol_o.checkerboard == Odd);
    }

    void low_mode_w(Matrix &_Matrix, const Field &evec, const RealD &eval, Field &wout)
    {
        GridBase *grid = _Matrix.RedBlackGrid();
        SchurDiagTwoOperator<Matrix, Field> _HermOpEO(_Matrix);

        Field src_o(grid);
        Field sol_e(grid);
        Field sol_o(grid);
        Field tmp(grid);

        GridBase *fgrid = _Matrix.Grid();
        Field tmp_out(fgrid);

        pickCheckerboard(Odd, src_o, evec);
        pickCheckerboard(Even, sol_e, wout);
        pickCheckerboard(Odd, sol_o, wout);

        /////////////////////////////////////////////////////
        // w_ie = MeeInvDag MoeDag Doo evec_i
        /////////////////////////////////////////////////////
        _HermOpEO.Mpc(src_o, sol_e);
        assert(sol_e.checkerboard == Odd);
        _Matrix.MeooeDag(sol_e, tmp);
        assert(tmp.checkerboard == Even);
        _Matrix.MooeeInvDag(tmp, sol_e);
        assert(sol_e.checkerboard == Even);

        /////////////////////////////////////////////////////
        // w_io = Doo evec_i
        /////////////////////////////////////////////////////
        _HermOpEO.Mpc(src_o, sol_o);
        assert(sol_o.checkerboard == Odd);

        setCheckerboard(tmp_out, sol_e);
        assert(sol_e.checkerboard == Even);
        setCheckerboard(tmp_out, sol_o);
        assert(sol_o.checkerboard == Odd);

        _Matrix.Dminus(tmp_out, wout);
    }

    void high_mode_v(Matrix &_Matrix, std::function<void(Field &, const Field &)> &Solver, const Field &source, Field &vout)
    {
        GridBase *fgrid = _Matrix.Grid();
        Field tmp(fgrid);

        _Matrix.Dminus(source, tmp);
        Solver(vout, tmp);
    }

    void high_mode_w(Matrix &_Matrix, const Field &source, Field &wout)
    {
        wout = source;
    }
};

////////////////////////////////
// Low Modes
////////////////////////////////

template <class Field, class Matrix>
class A2ALMSchurDiagTwo : public A2AModesSchurDiagTwo<Field, Matrix>
{
  private:
    const std::vector<Field> &evec;
    const std::vector<RealD> &eval;
    Matrix                   &action;

  public:
    A2ALMSchurDiagTwo(const std::vector<Field> &_evec, const std::vector<RealD> &_eval, Matrix &_action) : evec(_evec), eval(_eval), action(_action){};
    void operator()(int i, Field &vout, Field &wout)
    {
        this->low_mode_v(action, evec[i], eval[i], vout);
        this->low_mode_w(action, evec[i], eval[i], wout);
    }
};

template <class FineField, class CoarseField, class Matrix>
class A2ALMSchurDiagTwoCoarse : public A2AModesSchurDiagTwo<FineField, Matrix>
{
  private:
    const std::vector<FineField>   &subspace;
    const std::vector<CoarseField> &evec_coarse;
    const std::vector<RealD>       &eval_coarse;
    Matrix                         &action;

  public:
    A2ALMSchurDiagTwoCoarse(const std::vector<FineField> &_subspace, const std::vector<CoarseField> &_evec_coarse, const std::vector<RealD> &_eval_coarse, Matrix &_action)
        : subspace(_subspace), evec_coarse(_evec_coarse), eval_coarse(_eval_coarse), action(_action){};

    void operator()(int i, FineField &vout, FineField &wout)
    {
        FineField prom_evec(subspace[0]._grid);
        blockPromote(evec_coarse[i], prom_evec, subspace);
        this->low_mode_v(action, prom_evec, eval_coarse[i], vout);
        this->low_mode_w(action, prom_evec, eval_coarse[i], wout);
    }
};

////////////////////////////////
// High Modes
////////////////////////////////

template <class Field, class Matrix>
class A2AHMSchurDiagTwo : virtual public A2AModesSchurDiagTwo<Field, Matrix>
{
  public:
    void operator()(Matrix &_Matrix, std::function<void(Field &, const Field &)> &Solver, const Field &source, Field &vout, Field &wout)
    {
        this->high_mode_v(_Matrix, Solver, source, vout);
        this->high_mode_w(_Matrix, source, wout);
    }
};

template <class Field, class Matrix>
class A2AVectorsReturnHigh : public A2AModesSchurDiagTwo<Field, Matrix>
{
  private:
    Matrix &action;
    std::function<void(Field &, const Field &)> &Solver;
    const int Nh, Ls;

  public:
    std::vector<Field> w_high, v_high;
    A2AVectorsReturnHigh(Matrix &_action,
                         std::function<void(Field &, const Field &)> &_Solver,
                         const int _Nh, const int _Ls)
        : action(_action),
          Solver(_Solver),
          Nh(_Nh), Ls(_Ls)
    {
        GridBase *fgrid = action.Grid();
        resize(Nh, fgrid);
    };

    void resize(const size_t size, GridBase *grid)
    {
        w_high.resize(size, grid);
        v_high.resize(size, grid);
    }

    void high_modes(Field &source, int i)
    {
        this->high_mode_v(action, Solver, source, v_high[i]);
        this->high_mode_w(action, source, w_high[i]);
    }

    void operator()(int i, Field &vout, Field &wout)
    {
        if (Ls > 1)
        {
            vout = v_high[i];
            wout = w_high[i];
        }
        else
        {
            action.ExportPhysicalFermionSolution(v_high[i], vout);
            action.ExportPhysicalFermionSolution(w_high[i], wout);
        }
    }
};

////////////////////////////////
// Both Modes
////////////////////////////////

template <class Field, class Matrix>
class A2AVectorsReturn : public A2AModesSchurDiagTwo<Field, Matrix>
{
  private:
    const std::vector<Field> &evec;
    const std::vector<RealD> &eval;
    Matrix &action;
    std::function<void(Field &, const Field &)> &Solver;
    const int Nl, Nh, Ls;

  public:
    std::vector<Field> w_high, v_high;

    A2AVectorsReturn(const std::vector<Field> &_evec, const std::vector<RealD> &_eval,
                     Matrix &_action,
                     std::function<void(Field &, const Field &)> &_Solver,
                     const int _Nl, const int _Nh, const int _Ls)
        : evec(_evec), eval(_eval),
          action(_action),
          Solver(_Solver),
          Nl(_Nl), Nh(_Nh), Ls(_Ls)
    {
        GridBase *fgrid = action.Grid();
        resize(Nh, fgrid);
    };

    void resize(const size_t size, GridBase *grid)
    {
        w_high.resize(size, grid);
        v_high.resize(size, grid);
    }

    void high_modes(Field &source, int i)
    {
        this->high_mode_v(action, Solver, source, v_high[i]);
        this->high_mode_w(action, source, w_high[i]);
    }

    void operator()(int i, Field &vout, Field &wout)
    {

        GridBase *fgrid = action.Grid();
        Field vtmp(fgrid);
        Field wtmp(fgrid);
        if (i < Nl)
        {   
            this->low_mode_v(action, evec[i], eval[i], vtmp);
            this->low_mode_w(action, evec[i], eval[i], wtmp);
        }
        else
        {
            vtmp = v_high[i-Nl];
            wtmp = w_high[i-Nl];
        }

        if (Ls > 1)
        {
            vout = vtmp;
            wout = wtmp;
        }
        else
        {
            action.ExportPhysicalFermionSolution(vtmp, vout);
            action.ExportPhysicalFermionSolution(wtmp, wout);
        }
    }
};

END_HADRONS_NAMESPACE

#endif // A2A_Vectors_hpp_