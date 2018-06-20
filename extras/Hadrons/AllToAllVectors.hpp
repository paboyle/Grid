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

    void Doo(Matrix &action, const Field &in, Field &out)
    {
        Field tmp(in._grid);

        action.MooeeInv(in, out);
        action.Meooe(out, tmp);
        action.MooeeInv(tmp, out);
        action.Meooe(out, tmp);

        axpy(out, -1.0, tmp, in);
    }

    void low_mode_v(Matrix &action, const Field &evec, const RealD &eval, Field &vout, bool return_5d = true)
    {

        GridBase *grid = action.RedBlackGrid();
        Field src_o(grid);
        Field sol_e(grid);
        Field sol_o(grid);
        Field tmp(grid);

        GridBase *fgrid = action.Grid();
        Field tmp_out(fgrid);

        src_o = evec;
        src_o.checkerboard = Odd;
        pickCheckerboard(Even, sol_e, tmp_out);
        pickCheckerboard(Odd, sol_o, tmp_out);

        /////////////////////////////////////////////////////
        // v_ie = -(1/eval_i) * MeeInv Meo MooInv evec_i
        /////////////////////////////////////////////////////
        action.MooeeInv(src_o, tmp);
        assert(tmp.checkerboard == Odd);
        action.Meooe(tmp, sol_e);
        assert(sol_e.checkerboard == Even);
        action.MooeeInv(sol_e, tmp);
        assert(tmp.checkerboard == Even);
        sol_e = (-1.0 / eval) * tmp;
        assert(sol_e.checkerboard == Even);

        /////////////////////////////////////////////////////
        // v_io = (1/eval_i) * MooInv evec_i
        /////////////////////////////////////////////////////
        action.MooeeInv(src_o, tmp);
        assert(tmp.checkerboard == Odd);
        sol_o = (1.0 / eval) * tmp;
        assert(sol_o.checkerboard == Odd);

        setCheckerboard(tmp_out, sol_e);
        assert(sol_e.checkerboard == Even);
        setCheckerboard(tmp_out, sol_o);
        assert(sol_o.checkerboard == Odd);

        this->return_dim(action, tmp_out, vout, return_5d);
    }

    void low_mode_w(Matrix &action, const Field &evec, const RealD &eval, Field &wout, bool return_5d = true)
    {
        GridBase *grid = action.RedBlackGrid();
        SchurDiagTwoOperator<Matrix, Field> _HermOpEO(action);

        Field src_o(grid);
        Field sol_e(grid);
        Field sol_o(grid);
        Field tmp(grid);

        GridBase *fgrid = action.Grid();
        Field tmp_out(fgrid);
        Field tmp_wout(fgrid);

        src_o = evec;
        src_o.checkerboard = Odd;
        pickCheckerboard(Even, sol_e, tmp_wout);
        pickCheckerboard(Odd, sol_o, tmp_wout);

        /////////////////////////////////////////////////////
        // w_ie = - MeeInvDag MoeDag Doo evec_i
        /////////////////////////////////////////////////////
        Doo(action, src_o, tmp);
        assert(tmp.checkerboard == Odd);
        action.MeooeDag(tmp, sol_e);
        assert(sol_e.checkerboard == Even);
        action.MooeeInvDag(sol_e, tmp);
        assert(tmp.checkerboard == Even);
        sol_e = (-1.0) * tmp;

        /////////////////////////////////////////////////////
        // w_io = Doo evec_i
        /////////////////////////////////////////////////////
        Doo(action, src_o, sol_o);
        assert(sol_o.checkerboard == Odd);

        setCheckerboard(tmp_wout, sol_e);
        assert(sol_e.checkerboard == Even);
        setCheckerboard(tmp_wout, sol_o);
        assert(sol_o.checkerboard == Odd);

        action.DminusDag(tmp_wout, tmp_out);
        this->return_dim(action, tmp_out, wout, return_5d);
    }

    void high_mode_v(Matrix &action, std::function<void(Field &, const Field &)> &Solver, const Field &source, Field &vout, bool return_5d = true)
    {
        GridBase *fgrid = action.Grid();
        Field tmp(fgrid);
        Field tmp_out(fgrid);

        action.Dminus(source, tmp);
        Solver(tmp_out, source); // Note: Solver is Solver(out, in)
        this->return_dim(action, tmp_out, vout, return_5d);
    }

    void high_mode_w(Matrix &action, const Field &source4d, Field &wout, bool return_5d = true)
    {
        // GridBase *fgrid = action.Grid();
        // Field tmp_out(fgrid);

        // tmp_out = source;
        // this->return_dim(action, tmp_out, wout, return_5d);
        wout = source4d;
    }

    void return_dim(Matrix &action, const Field &in, Field &out, bool return_5d)
    {
        if (return_5d)
        {
            out = in;
        }
        else
        {
            action.ExportPhysicalFermionSolution(in, out);
        }
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
    Matrix &action;

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
    const std::vector<FineField> &subspace;
    const std::vector<CoarseField> &evec_coarse;
    const std::vector<RealD> &eval_coarse;
    Matrix &action;

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
    void operator()(Matrix &action, std::function<void(Field &, const Field &)> &Solver, const Field &source, Field &vout, Field &wout)
    {
        this->high_mode_v(action, Solver, source, vout);
        this->high_mode_w(action, source, wout);
    }
};

////////////////////////////////
// Both Modes
////////////////////////////////

template <class Field, class Matrix>
class A2AVectorsReturn : public A2AModesSchurDiagTwo<Field, Matrix>
{
  private:
    const std::vector<Field> *evec;
    const std::vector<RealD> *eval;
    Matrix &action;
    std::function<void(Field &, const Field &)> &Solver;
    const int Nl, Nh;
    const bool return_5d;
    std::vector<Field> w_high, v_high;

  public:
    A2AVectorsReturn(const std::vector<Field> *_evec, const std::vector<RealD> *_eval,
                     Matrix &_action,
                     std::function<void(Field &, const Field &)> &_Solver,
                     const int _Nl, const int _Nh,
                     const bool _return_5d)
        : evec(_evec), eval(_eval),
          action(_action),
          Solver(_Solver),
          Nl(_Nl), Nh(_Nh),
          return_5d(_return_5d)
    {
        GridBase *grid;
        if (return_5d)
        {
            grid = action.Grid();
        }
        else
        {
            grid = action.GaugeGrid();
        }
        resize(Nh, grid);
    };

    void resize(const size_t size, GridBase *grid)
    {
        w_high.resize(size, grid);
        v_high.resize(size, grid);
    }

    void high_modes(Field &source5d, Field &source4d, int i)
    {
        LOG(Message) << "A2A high modes for i = " << i << std::endl;
        this->high_mode_v(action, Solver, source5d, v_high[i], return_5d);
        this->high_mode_w(action, source4d, w_high[i], return_5d);
    }

    void operator()(int i, Field &vout, Field &wout)
    {
        if (i < Nl)
        {
            LOG(Message) << "A2A low modes for i = " << i << std::endl;
            this->low_mode_v(action, evec->at(i), eval->at(i), vout, return_5d);
            this->low_mode_w(action, evec->at(i), eval->at(i), wout, return_5d);
        }
        else
        {
            vout = v_high[i - Nl];
            wout = w_high[i - Nl];
        }
    }
};

END_HADRONS_NAMESPACE

#endif // A2A_Vectors_hpp_