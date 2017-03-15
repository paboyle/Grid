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

namespace Grid {
  // FIXME drop the QCD namespace everywhere here

template <class Impl>
class ScalarInteractionAction : public QCD::Action<typename Impl::Field> {
public:
    INHERIT_FIELD_TYPES(Impl);
private:
    RealD mass_square;
    RealD lambda;


    typedef typename Field::vector_object vobj;
    typedef CartesianStencil<vobj,vobj> Stencil;

    SimpleCompressor<vobj> compressor;
    int npoint = 8;
    std::vector<int> directions    = {0,1,2,3,0,1,2,3};  // forcing 4 dimensions
    std::vector<int> displacements = {1,1,1,1, -1,-1,-1,-1};


 public:

    ScalarInteractionAction(RealD ms, RealD l) : mass_square(ms), lambda(l){}

    virtual std::string LogParameters() {
      std::stringstream sstream;
      sstream << GridLogMessage << "[ScalarAction] lambda      : " << lambda      << std::endl;
      sstream << GridLogMessage << "[ScalarAction] mass_square : " << mass_square << std::endl;
      return sstream.str();
    }

    virtual std::string action_name() {return "ScalarAction";}

    virtual void refresh(const Field &U, GridParallelRNG &pRNG) {}

    virtual RealD S(const Field &p) {
        static Stencil phiStencil(p._grid, npoint, 0, directions, displacements);
        phiStencil.HaloExchange(p, compressor);

        Field action(p._grid), pshift(p._grid), phisquared(p._grid);
        phisquared = p*p;
        action = (2.0*QCD::Nd + mass_square)*phisquared + lambda*phisquared*phisquared;
        for (int mu = 0; mu < QCD::Nd; mu++) {
            //  pshift = Cshift(p, mu, +1);  // not efficient, implement with stencils
            PARALLEL_FOR_LOOP
            for (int i = 0; i < p._grid->oSites(); i++) {
                int permute_type;
                StencilEntry *SE;
                vobj temp2;
                vobj *temp;
                vobj *t_p;

                SE = phiStencil.GetEntry(permute_type, mu, i);
                t_p  = &p._odata[i];
                if ( SE->_is_local ) {
                    temp = &p._odata[SE->_offset];
                    if ( SE->_permute ) {
                        permute(temp2, *temp, permute_type);
                        action._odata[i] -= temp2*(*t_p) + (*t_p)*temp2;
                    } else {
                  action._odata[i] -= *temp*(*t_p) + (*t_p)*(*temp);
                    }
                } else {
                  action._odata[i] -= phiStencil.CommBuf()[SE->_offset]*(*t_p) + (*t_p)*phiStencil.CommBuf()[SE->_offset];
                }
            }
            //  action -= pshift*p + p*pshift;
        }
        // NB the trace in the algebra is normalised to 1/2
        // minus sign coming from the antihermitian fields
        return -(TensorRemove(sum(trace(action)))).real();
    };

    virtual void deriv(const Field &p, Field &force) {
        force = (2.0*QCD::Nd + mass_square)*p + 2.0*lambda*p*p*p;
        // move this outside
        static Stencil phiStencil(p._grid, npoint, 0, directions, displacements);
        phiStencil.HaloExchange(p, compressor);

        //for (int mu = 0; mu < QCD::Nd; mu++) force -= Cshift(p, mu, -1) + Cshift(p, mu, 1);
        for (int point = 0; point < npoint; point++) {
            PARALLEL_FOR_LOOP
            for (int i = 0; i < p._grid->oSites(); i++) {
                vobj *temp;
                vobj temp2;
                int permute_type;
                StencilEntry *SE;
                SE = phiStencil.GetEntry(permute_type, point, i);

                if ( SE->_is_local ) {
                    temp = &p._odata[SE->_offset];
                    if ( SE->_permute ) {
                        permute(temp2, *temp, permute_type);
                        force._odata[i] -= temp2;
                    } else {
                        force._odata[i] -= *temp;
                    }
                } else {
                    force._odata[i] -= phiStencil.CommBuf()[SE->_offset];
                }
            }
        }
    }
};

}  // namespace Grid

#endif  // SCALAR_INT_ACTION_H
