    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/action/fermion/ContinuedFractionFermion5D.h

    Copyright (C) 2015

Author: Peter Boyle <pabobyle@ph.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>

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
#ifndef  GRID_QCD_CONTINUED_FRACTION_H
#define  GRID_QCD_CONTINUED_FRACTION_H

namespace Grid {

  namespace QCD {

    template<class Impl>
    class ContinuedFractionFermion5D : public WilsonFermion5D<Impl>
    {
    public:
     INHERIT_IMPL_TYPES(Impl);
    public:

      // override multiply
      virtual RealD  M    (const FermionField &in, FermionField &out);
      virtual RealD  Mdag (const FermionField &in, FermionField &out);

      // half checkerboard operaions
      virtual void   Meooe       (const FermionField &in, FermionField &out);
      virtual void   MeooeDag    (const FermionField &in, FermionField &out);
      virtual void   Mooee       (const FermionField &in, FermionField &out);
      virtual void   MooeeDag    (const FermionField &in, FermionField &out);
      virtual void   MooeeInv    (const FermionField &in, FermionField &out);
      virtual void   MooeeInvDag (const FermionField &in, FermionField &out);

      // force terms; five routines; default to Dhop on diagonal
      virtual void MDeriv  (GaugeField &mat,const FermionField &U,const FermionField &V,int dag);
      virtual void MoeDeriv(GaugeField &mat,const FermionField &U,const FermionField &V,int dag);
      virtual void MeoDeriv(GaugeField &mat,const FermionField &U,const FermionField &V,int dag);

      //      virtual void   Instantiatable(void)=0;
      virtual void   Instantiatable(void) =0;

      // Efficient support for multigrid coarsening
      virtual void  Mdir (const FermionField &in, FermionField &out,int dir,int disp);

      // Constructors
      ContinuedFractionFermion5D(GaugeField &_Umu,
				 GridCartesian         &FiveDimGrid,
				 GridRedBlackCartesian &FiveDimRedBlackGrid,
				 GridCartesian         &FourDimGrid,
				 GridRedBlackCartesian &FourDimRedBlackGrid,
				 RealD _mass,RealD M5,const ImplParams &p= ImplParams());

    protected:

      void SetCoefficientsTanh(Approx::zolotarev_data *zdata,RealD scale);
      void SetCoefficientsZolotarev(RealD zolo_hi,Approx::zolotarev_data *zdata);;

      // Cont frac
      RealD dw_diag;
      RealD mass;
      RealD R;
      RealD ZoloHiInv;
      std::vector<double> Beta;
      std::vector<double> cc;;
      std::vector<double> cc_d;;
      std::vector<double> sqrt_cc;
      std::vector<double> See;
      std::vector<double> Aee;

    };


  }
}

#endif
