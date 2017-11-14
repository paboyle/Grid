/*************************************************************************************

  Grid physics library, www.github.com/paboyle/Grid

  Source file: ./lib/qcd/action/gauge/WilsonGaugeAction.h

  Copyright (C) 2015

  Author: Guido Cossu <guido,cossu@ed.ac.uk>

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

  See the full license in the file "LICENSE" in the top level distribution
directory
  *************************************************************************************/
/*  END LEGAL */

#ifndef SCALAR_INT_ACTION_H
#define SCALAR_INT_ACTION_H

// Note: this action can completely absorb the ScalarAction for real float fields
// use the scalarObjs to generalise the structure

namespace Grid
{
// FIXME drop the QCD namespace everywhere here

template <class Impl, int Ndim>
class ScalarInteractionAction : public QCD::Action<typename Impl::Field>
{
public:
  INHERIT_FIELD_TYPES(Impl);

private:
  RealD mass_square;
  RealD lambda;
  RealD g;
  const unsigned int N = Impl::Group::Dimension;

  typedef typename Field::vector_object vobj;
  typedef CartesianStencil<vobj, vobj> Stencil;

  SimpleCompressor<vobj> compressor;
  int npoint = 2 * Ndim;
  std::vector<int> directions;    //
  std::vector<int> displacements; //

public:
  ScalarInteractionAction(RealD ms, RealD l, RealD gval) : mass_square(ms), lambda(l), g(gval), displacements(2 * Ndim, 0), directions(2 * Ndim, 0)
  {
    for (int mu = 0; mu < Ndim; mu++)
    {
      directions[mu] = mu;
      directions[mu + Ndim] = mu;
      displacements[mu] = 1;
      displacements[mu + Ndim] = -1;
    }
  }

  virtual std::string LogParameters()
  {
    std::stringstream sstream;
    sstream << GridLogMessage << "[ScalarAction] lambda      : " << lambda << std::endl;
    sstream << GridLogMessage << "[ScalarAction] mass_square : " << mass_square << std::endl;
    sstream << GridLogMessage << "[ScalarAction] g           : " << g << std::endl;
    return sstream.str();
  }

  virtual std::string action_name() { return "ScalarAction"; }

  virtual void refresh(const Field &U, GridParallelRNG &pRNG) {}

  virtual RealD S(const Field &p)
  {
    assert(p._grid->Nd() == Ndim);
    static Stencil phiStencil(p._grid, npoint, 0, directions, displacements);
    phiStencil.HaloExchange(p, compressor);
    Field action(p._grid), pshift(p._grid), phisquared(p._grid);
    phisquared = p * p;
    action = (2.0 * Ndim + mass_square) * phisquared - lambda * phisquared * phisquared;
    for (int mu = 0; mu < Ndim; mu++)
    {
      //  pshift = Cshift(p, mu, +1);  // not efficient, implement with stencils
      parallel_for(int i = 0; i < p._grid->oSites(); i++)
      {
        int permute_type;
        StencilEntry *SE;
        vobj temp2;
        const vobj *temp, *t_p;

        SE = phiStencil.GetEntry(permute_type, mu, i);
        t_p = &p._odata[i];
        if (SE->_is_local)
        {
          temp = &p._odata[SE->_offset];
          if (SE->_permute)
          {
            permute(temp2, *temp, permute_type);
            action._odata[i] -= temp2 * (*t_p) + (*t_p) * temp2;
          }
          else
          {
            action._odata[i] -= (*temp) * (*t_p) + (*t_p) * (*temp);
          }
        }
        else
        {
          action._odata[i] -= phiStencil.CommBuf()[SE->_offset] * (*t_p) + (*t_p) * phiStencil.CommBuf()[SE->_offset];
        }
      }
      //  action -= pshift*p + p*pshift;
    }
    // NB the trace in the algebra is normalised to 1/2
    // minus sign coming from the antihermitian fields
    return -(TensorRemove(sum(trace(action)))).real() * N / g;
  };

  virtual void deriv(const Field &p, Field &force)
  {
    double t0 = usecond();
    assert(p._grid->Nd() == Ndim);
    force = (2. * Ndim + mass_square) * p - 2. * lambda * p * p * p;
    double interm_t = usecond();

    // move this outside
    static Stencil phiStencil(p._grid, npoint, 0, directions, displacements);

    phiStencil.HaloExchange(p, compressor);
    double halo_t = usecond();
    int chunk = 128;
    //for (int mu = 0; mu < QCD::Nd; mu++) force -= Cshift(p, mu, -1) + Cshift(p, mu, 1);

    // inverting the order of the loops slows down the code(! g++ 7)
    // cannot try to reduce the number of  force writes by factor npoint...
    // use cache blocking
    for (int point = 0; point < npoint; point++)
    {

#pragma omp parallel 
{
        int permute_type;
        StencilEntry *SE;
        const vobj *temp;

#pragma omp for schedule(static, chunk)
      for (int i = 0; i < p._grid->oSites(); i++)
      {
        SE = phiStencil.GetEntry(permute_type, point, i);
        // prefetch next p?

        if (SE->_is_local)
        {
          temp = &p._odata[SE->_offset];
      
          if (SE->_permute)
          {
            vobj temp2;
            permute(temp2, *temp, permute_type);
            force._odata[i] -= temp2;
          }
          else
          {
            force._odata[i] -= *temp; // slow part. Dominated by this read/write (BW)
          }
        }
        else
        {
          force._odata[i] -= phiStencil.CommBuf()[SE->_offset];
        }
      }

    }
  }
  force *= N / g;

  double t1 = usecond();
  double total_time = (t1 - t0) / 1e6;
  double interm_time = (interm_t - t0) / 1e6;
  double halo_time = (halo_t - interm_t) / 1e6;
  double stencil_time = (t1 - halo_t) / 1e6;
  std::cout << GridLogIntegrator << "Total time for force computation (s)       : " << total_time << std::endl;
  std::cout << GridLogIntegrator << "Intermediate time for force computation (s): " << interm_time << std::endl;
  std::cout << GridLogIntegrator << "Halo time in force computation (s)         : " << halo_time << std::endl;
  std::cout << GridLogIntegrator << "Stencil time in force computation (s)      : " << stencil_time << std::endl;
  double flops = p._grid->gSites() * (14 * N * N * N + 18 * N * N + 2);
  double flops_no_stencil = p._grid->gSites() * (14 * N * N * N + 6 * N * N + 2);
  double Gflops = flops / (total_time * 1e9);
  double Gflops_no_stencil = flops_no_stencil / (interm_time * 1e9);
  std::cout << GridLogIntegrator << "Flops: " << flops << "  - Gflop/s : " << Gflops << std::endl;
  std::cout << GridLogIntegrator << "Flops NS: " << flops_no_stencil << "  - Gflop/s NS: " << Gflops_no_stencil << std::endl;
}
};

} // namespace Grid

#endif // SCALAR_INT_ACTION_H
