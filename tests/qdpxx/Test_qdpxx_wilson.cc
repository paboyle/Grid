/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/qdpxx/Test_qdpxx_wilson.cc

    Copyright (C) 2017

    Author: Guido Cossu <guido.cossu@ed.ac.uk>

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

#include <Grid/Grid.h>
#include <chroma.h>
#include <actions/ferm/invert/syssolver_linop_cg_array.h>
#include <actions/ferm/invert/syssolver_linop_aggregate.h>

// Mass
double mq = 0.1;

// Define Wilson Types
typedef Grid::QCD::WilsonImplR::FermionField FermionField;
typedef Grid::QCD::LatticeGaugeField GaugeField;

enum ChromaAction
{
  Wilson,      // Wilson
  WilsonClover // CloverFermions
};

namespace Chroma
{

class ChromaWrapper
{
public:
  typedef multi1d<LatticeColorMatrix> U;
  typedef LatticeFermion T4;

  static void ImportGauge(GaugeField &gr,
                          QDP::multi1d<QDP::LatticeColorMatrix> &ch)
  {
    Grid::QCD::LorentzColourMatrix LCM;
    Grid::Complex cc;
    QDP::ColorMatrix cm;
    QDP::Complex c;

    std::vector<int> x(4);
    QDP::multi1d<int> cx(4);
    std::vector<int> gd = gr._grid->GlobalDimensions();

    for (x[0] = 0; x[0] < gd[0]; x[0]++)
    {
      for (x[1] = 0; x[1] < gd[1]; x[1]++)
      {
        for (x[2] = 0; x[2] < gd[2]; x[2]++)
        {
          for (x[3] = 0; x[3] < gd[3]; x[3]++)
          {
            cx[0] = x[0];
            cx[1] = x[1];
            cx[2] = x[2];
            cx[3] = x[3];
            Grid::peekSite(LCM, gr, x);

            for (int mu = 0; mu < 4; mu++)
            {
              for (int i = 0; i < 3; i++)
              {
                for (int j = 0; j < 3; j++)
                {
                  cc = LCM(mu)()(i, j);
                  c = QDP::cmplx(QDP::Real(real(cc)), QDP::Real(imag(cc)));
                  QDP::pokeColor(cm, c, i, j);
                }
              }
              QDP::pokeSite(ch[mu], cm, cx);
            }
          }
        }
      }
    }
  }

  static void ExportGauge(GaugeField &gr,
                          QDP::multi1d<QDP::LatticeColorMatrix> &ch)
  {
    Grid::QCD::LorentzColourMatrix LCM;
    Grid::Complex cc;
    QDP::ColorMatrix cm;
    QDP::Complex c;

    std::vector<int> x(4);
    QDP::multi1d<int> cx(4);
    std::vector<int> gd = gr._grid->GlobalDimensions();

    for (x[0] = 0; x[0] < gd[0]; x[0]++)
    {
      for (x[1] = 0; x[1] < gd[1]; x[1]++)
      {
        for (x[2] = 0; x[2] < gd[2]; x[2]++)
        {
          for (x[3] = 0; x[3] < gd[3]; x[3]++)
          {
            cx[0] = x[0];
            cx[1] = x[1];
            cx[2] = x[2];
            cx[3] = x[3];

            for (int mu = 0; mu < 4; mu++)
            {
              for (int i = 0; i < 3; i++)
              {
                for (int j = 0; j < 3; j++)
                {
                  cm = QDP::peekSite(ch[mu], cx);
                  c = QDP::peekColor(cm, i, j);
                  cc = Grid::Complex(toDouble(real(c)), toDouble(imag(c)));
                  LCM(mu)
                  ()(i, j) = cc;
                }
              }
            }
            Grid::pokeSite(LCM, gr, x);
          }
        }
      }
    }
  }

  // Specific for Wilson Fermions
  static void ImportFermion(Grid::QCD::LatticeFermion &gr,
                            QDP::LatticeFermion &ch)
  {
    Grid::QCD::SpinColourVector F;
    Grid::Complex c;

    QDP::Fermion cF;
    QDP::SpinVector cS;
    QDP::Complex cc;

    std::vector<int> x(4); // explicit 4d fermions in Grid
    QDP::multi1d<int> cx(4);
    std::vector<int> gd = gr._grid->GlobalDimensions();

    for (x[0] = 0; x[0] < gd[0]; x[0]++)
    {
      for (x[1] = 0; x[1] < gd[1]; x[1]++)
      {
        for (x[2] = 0; x[2] < gd[2]; x[2]++)
        {
          for (x[3] = 0; x[3] < gd[3]; x[3]++)
          {
            cx[0] = x[0];
            cx[1] = x[1];
            cx[2] = x[2];
            cx[3] = x[3];

            Grid::peekSite(F, gr, x);

            for (int j = 0; j < 3; j++)
            {
              for (int sp = 0; sp < 4; sp++)
              {

                c = F()(sp)(j);

                cc = QDP::cmplx(QDP::Real(real(c)), QDP::Real(imag(c)));

                QDP::pokeSpin(cS, cc, sp);
              }
              QDP::pokeColor(cF, cS, j);
            }
            QDP::pokeSite(ch, cF, cx);
          }
        }
      }
    }
  }

  // Specific for 4d Wilson fermions
  static void ExportFermion(Grid::QCD::LatticeFermion &gr,
                            QDP::LatticeFermion &ch)
  {
    Grid::QCD::SpinColourVector F;
    Grid::Complex c;

    QDP::Fermion cF;
    QDP::SpinVector cS;
    QDP::Complex cc;

    std::vector<int> x(4); // 4d fermions
    QDP::multi1d<int> cx(4);
    std::vector<int> gd = gr._grid->GlobalDimensions();

    for (x[0] = 0; x[0] < gd[0]; x[0]++)
    {
      for (x[1] = 0; x[1] < gd[1]; x[1]++)
      {
        for (x[2] = 0; x[2] < gd[2]; x[2]++)
        {
          for (x[3] = 0; x[3] < gd[3]; x[3]++)
          {
            cx[0] = x[0];
            cx[1] = x[1];
            cx[2] = x[2];
            cx[3] = x[3];

            cF = QDP::peekSite(ch, cx);
            for (int sp = 0; sp < 4; sp++)
            {
              for (int j = 0; j < 3; j++)
              {
                cS = QDP::peekColor(cF, j);
                cc = QDP::peekSpin(cS, sp);
                c = Grid::Complex(QDP::toDouble(QDP::real(cc)),
                                  QDP::toDouble(QDP::imag(cc)));
                F()
                (sp)(j) = c;
              }
            }
            Grid::pokeSite(F, gr, x);
          }
        }
      }
    }
  }

  static Handle<Chroma::UnprecLinearOperator<T4, U, U>> GetLinOp(U &u, ChromaAction params)
  {
    QDP::Real _mq(mq);
    QDP::multi1d<int> bcs(QDP::Nd);

    // Boundary conditions
    bcs[0] = bcs[1] = bcs[2] = bcs[3] = 1;

    if (params == Wilson)
    {

      Chroma::WilsonFermActParams p;
      p.Mass = _mq;
      AnisoParam_t _apar;
      _apar.anisoP = true;
      _apar.t_dir = 3; // in 4d
      _apar.xi_0 = 2.0;
      _apar.nu = 1.0;
      p.anisoParam = _apar;

      Chroma::Handle<Chroma::FermBC<T4, U, U>> fbc(new Chroma::SimpleFermBC<T4, U, U>(bcs));
      Chroma::Handle<Chroma::CreateFermState<T4, U, U>> cfs(new Chroma::CreateSimpleFermState<T4, U, U>(fbc));
      Chroma::UnprecWilsonFermAct S_f(cfs, p);
      Chroma::Handle<Chroma::FermState<T4, U, U>> ffs(S_f.createState(u));
      return S_f.linOp(ffs);
    }

    if (params == WilsonClover)
    {
      Chroma::CloverFermActParams p;
      p.Mass = _mq;
      p.clovCoeffR = QDP::Real(1.0);
      p.clovCoeffT = QDP::Real(2.0);
      p.u0 = QDP::Real(1.0);
      AnisoParam_t _apar;
      _apar.anisoP = true;
      _apar.t_dir = 3; // in 4d
      _apar.xi_0 = 2.0;
      _apar.nu = 1.0;
      p.anisoParam = _apar;

      Chroma::Handle<Chroma::FermBC<T4, U, U>> fbc(new Chroma::SimpleFermBC<T4, U, U>(bcs));
      Chroma::Handle<Chroma::CreateFermState<T4, U, U>> cfs(new Chroma::CreateSimpleFermState<T4, U, U>(fbc));
      Chroma::UnprecCloverFermAct S_f(cfs, p);
      Chroma::Handle<Chroma::FermState<T4, U, U>> ffs(S_f.createState(u));
      return S_f.linOp(ffs);
    }
  }
};
} // namespace Chroma

void calc_chroma(ChromaAction action, GaugeField &lat, FermionField &src, FermionField &res, int dag)
{
  QDP::multi1d<QDP::LatticeColorMatrix> u(4);
  Chroma::ChromaWrapper::ImportGauge(lat, u);

  QDP::LatticeFermion check;
  QDP::LatticeFermion result;
  QDP::LatticeFermion psi;

  Chroma::ChromaWrapper::ImportFermion(src, psi);

  for (int mu = 0; mu < 4; mu++)
  {
    std::cout << "Imported Gauge norm [" << mu << "] " << QDP::norm2(u[mu]) << std::endl;
  }
  std::cout << "Imported Fermion norm " << QDP::norm2(psi) << std::endl;

  typedef QDP::LatticeFermion T;
  typedef QDP::multi1d<QDP::LatticeColorMatrix> U;

  auto linop = Chroma::ChromaWrapper::GetLinOp(u, action);

  printf("Calling Chroma Linop\n");
  fflush(stdout);

  if (dag)
    (*linop)(check, psi, Chroma::MINUS);
  else
    (*linop)(check, psi, Chroma::PLUS);

  printf("Called Chroma Linop\n");
  fflush(stdout);

  // std::cout << "Calling Chroma Linop " << std::endl;
  // linop->evenEvenLinOp(tmp, psi, isign);
  // check[rb[0]] = tmp;
  // linop->oddOddLinOp(tmp, psi, isign);
  // check[rb[1]] = tmp;
  // linop->evenOddLinOp(tmp, psi, isign);
  // check[rb[0]] += tmp;
  // linop->oddEvenLinOp(tmp, psi, isign);
  // check[rb[1]] += tmp;

  Chroma::ChromaWrapper::ExportFermion(res, check);
}

void make_gauge(GaugeField &Umu, FermionField &src)
{
  using namespace Grid;
  using namespace Grid::QCD;

  std::vector<int> seeds4({1, 2, 3, 4});

  Grid::GridCartesian *UGrid = (Grid::GridCartesian *)Umu._grid;
  Grid::GridParallelRNG RNG4(UGrid);
  RNG4.SeedFixedIntegers(seeds4);
  Grid::QCD::SU3::HotConfiguration(RNG4, Umu);

  // Fermion field
  Grid::gaussian(RNG4, src);
  /*
  Grid::QCD::SpinColourVector F;
  Grid::Complex c;

  

  std::vector<int> x(4); // 4d fermions
  std::vector<int> gd = src._grid->GlobalDimensions();

  for (x[0] = 0; x[0] < gd[0]; x[0]++)
  {
    for (x[1] = 0; x[1] < gd[1]; x[1]++)
    {
      for (x[2] = 0; x[2] < gd[2]; x[2]++)
      {
        for (x[3] = 0; x[3] < gd[3]; x[3]++)
        {
          for (int sp = 0; sp < 4; sp++)
          {
            for (int j = 0; j < 3; j++) // colours
            {
              F()(sp)(j) = Grid::Complex(0.0,0.0);
              if (((sp == 0)|| (sp==3)) && (j==2))
              {
                c = Grid::Complex(1.0, 0.0);
                F()(sp)(j) = c;
              }
            }
          }
          Grid::pokeSite(F, src, x);
          
        }
      }
    }
  }
  */
}

void calc_grid(ChromaAction action, Grid::QCD::LatticeGaugeField &Umu, Grid::QCD::LatticeFermion &src, Grid::QCD::LatticeFermion &res, int dag)
{
  using namespace Grid;
  using namespace Grid::QCD;

  Grid::GridCartesian *UGrid = (Grid::GridCartesian *)Umu._grid;
  Grid::GridRedBlackCartesian *UrbGrid = Grid::QCD::SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);

  Grid::RealD _mass = mq;

  if (action == Wilson)
  {
    WilsonAnisotropyCoefficients anis;
    anis.isAnisotropic = true;
    anis.t_direction = 3;
    anis.xi_0 = 2.0;
    anis.nu = 1.0;
    WilsonImplParams iParam;
    Grid::QCD::WilsonFermionR Wf(Umu, *UGrid, *UrbGrid, _mass, iParam, anis);

    std::cout << Grid::GridLogMessage << " Calling Grid Wilson Fermion multiply " << std::endl;

    if (dag)
      Wf.Mdag(src, res);
    else
      Wf.M(src, res);
    return;
  }

  if (action == WilsonClover)
  {
    Grid::RealD _csw_r = 1.0;
    Grid::RealD _csw_t = 2.0;
    WilsonAnisotropyCoefficients anis;
    anis.isAnisotropic = true;
    anis.t_direction = 3;
    anis.xi_0 = 2.0;
    anis.nu = 1.0;
    WilsonImplParams CloverImplParam;
    Grid::QCD::WilsonCloverFermionR Wf(Umu, *UGrid, *UrbGrid, _mass, _csw_r, _csw_t, anis, CloverImplParam);
    Wf.ImportGauge(Umu);

    std::cout << Grid::GridLogMessage << " Calling Grid Wilson Clover Fermion multiply " << std::endl;

    if (dag)
      Wf.Mdag(src, res);
    else
      Wf.M(src, res);
    return;
  }

  assert(0);
}

int main(int argc, char **argv)
{

  /********************************************************
   * Setup QDP
   *********************************************************/
  Chroma::initialize(&argc, &argv);
  Chroma::WilsonTypeFermActs4DEnv::registerAll();

  /********************************************************
   * Setup Grid
   *********************************************************/
  Grid::Grid_init(&argc, &argv);
  Grid::GridCartesian *UGrid = Grid::QCD::SpaceTimeGrid::makeFourDimGrid(Grid::GridDefaultLatt(),
                                                                         Grid::GridDefaultSimd(Grid::QCD::Nd, Grid::vComplex::Nsimd()),
                                                                         Grid::GridDefaultMpi());

  std::vector<int> gd = UGrid->GlobalDimensions();
  QDP::multi1d<int> nrow(QDP::Nd);
  for (int mu = 0; mu < 4; mu++)
    nrow[mu] = gd[mu];

  QDP::Layout::setLattSize(nrow);
  QDP::Layout::create();

  GaugeField Ug(UGrid);
  FermionField src(UGrid);
  FermionField res_chroma(UGrid);
  FermionField res_grid(UGrid);
  FermionField only_wilson(UGrid);
  FermionField difference(UGrid);

  std::vector<ChromaAction> ActionList({Wilson, WilsonClover});
  std::vector<std::string> ActionName({"Wilson", "WilsonClover"});

  {

    for (int i = 0; i < ActionList.size(); i++)
    {
      std::cout << "*****************************" << std::endl;
      std::cout << "Action " << ActionName[i] << std::endl;
      std::cout << "*****************************" << std::endl;
      make_gauge(Ug, src); // fills the gauge field and the fermion field with random numbers

      for (int dag = 0; dag < 2; dag++)
      {

        {

          std::cout << "Dag =  " << dag << std::endl;

          calc_chroma(ActionList[i], Ug, src, res_chroma, dag);

          // Remove the normalisation of Chroma Gauge links ????????
          std::cout << "Norm of Chroma " << ActionName[i] << " multiply " << Grid::norm2(res_chroma) << std::endl;
          calc_grid(ActionList[i], Ug, src, res_grid, dag);

          std::cout << "Norm of gauge " << Grid::norm2(Ug) << std::endl;

          std::cout << "Norm of Grid " << ActionName[i] << " multiply " << Grid::norm2(res_grid) << std::endl;

          difference = res_chroma - res_grid;
          std::cout << "Norm of difference " << Grid::norm2(difference) << std::endl;
        }
      }

      std::cout << "Finished test " << std::endl;

      Chroma::finalize();
    }
  }
}
