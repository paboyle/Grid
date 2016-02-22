    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/action/fermion/FermionOperator.h

    Copyright (C) 2015

Author: Peter Boyle <pabobyle@ph.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: Peter Boyle <peterboyle@Peters-MacBook-Pro-2.local>

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
#ifndef  GRID_QCD_FERMION_OPERATOR_H
#define  GRID_QCD_FERMION_OPERATOR_H

namespace Grid {

  namespace QCD {

    ////////////////////////////////////////////////////////////////
    // Allow to select  between gauge representation rank bc's, flavours etc.
    // and single/double precision.
    ////////////////////////////////////////////////////////////////
    
    template<class Impl>
    class FermionOperator : public CheckerBoardedSparseMatrixBase<typename Impl::FermionField>, public Impl
    {
    public:

      INHERIT_IMPL_TYPES(Impl);

      FermionOperator(const ImplParams &p= ImplParams()) : Impl(p) {};

      GridBase * Grid(void)   { return FermionGrid(); };   // this is all the linalg routines need to know
      GridBase * RedBlackGrid(void) { return FermionRedBlackGrid(); };

      virtual GridBase *FermionGrid(void)         =0;
      virtual GridBase *FermionRedBlackGrid(void) =0;
      virtual GridBase *GaugeGrid(void)           =0;
      virtual GridBase *GaugeRedBlackGrid(void)   =0;

      // override multiply
      virtual RealD  M    (const FermionField &in, FermionField &out)=0;
      virtual RealD  Mdag (const FermionField &in, FermionField &out)=0;

      // half checkerboard operaions
      virtual int    ConstEE(void) { return 1; }; // clover returns zero as EE depends on gauge field

      virtual void   Meooe       (const FermionField &in, FermionField &out)=0;
      virtual void   MeooeDag    (const FermionField &in, FermionField &out)=0;
      virtual void   Mooee       (const FermionField &in, FermionField &out)=0;
      virtual void   MooeeDag    (const FermionField &in, FermionField &out)=0;
      virtual void   MooeeInv    (const FermionField &in, FermionField &out)=0;
      virtual void   MooeeInvDag (const FermionField &in, FermionField &out)=0;

      // non-hermitian hopping term; half cb or both
      virtual void Dhop  (const FermionField &in, FermionField &out,int dag)=0;
      virtual void DhopOE(const FermionField &in, FermionField &out,int dag)=0;
      virtual void DhopEO(const FermionField &in, FermionField &out,int dag)=0;
      virtual void DhopDir(const FermionField &in, FermionField &out,int dir,int disp)=0; // implemented by WilsonFermion and WilsonFermion5D

      // force terms; five routines; default to Dhop on diagonal
      virtual void MDeriv  (GaugeField &mat,const FermionField &U,const FermionField &V,int dag){DhopDeriv(mat,U,V,dag);};
      virtual void MoeDeriv(GaugeField &mat,const FermionField &U,const FermionField &V,int dag){DhopDerivOE(mat,U,V,dag);};
      virtual void MeoDeriv(GaugeField &mat,const FermionField &U,const FermionField &V,int dag){DhopDerivEO(mat,U,V,dag);};
      virtual void MooDeriv(GaugeField &mat,const FermionField &U,const FermionField &V,int dag){mat=zero;}; // Clover can override these
      virtual void MeeDeriv(GaugeField &mat,const FermionField &U,const FermionField &V,int dag){mat=zero;};

      virtual void DhopDeriv  (GaugeField &mat,const FermionField &U,const FermionField &V,int dag)=0;
      virtual void DhopDerivEO(GaugeField &mat,const FermionField &U,const FermionField &V,int dag)=0;
      virtual void DhopDerivOE(GaugeField &mat,const FermionField &U,const FermionField &V,int dag)=0;


      virtual void  Mdiag  (const FermionField &in, FermionField &out) { Mooee(in,out);};   // Same as Mooee applied to both CB's
      virtual void  Mdir   (const FermionField &in, FermionField &out,int dir,int disp)=0;   // case by case Wilson, Clover, Cayley, ContFrac, PartFrac

      ///////////////////////////////////////////////
      // Updates gauge field during HMC
      ///////////////////////////////////////////////
      virtual void ImportGauge(const GaugeField & _U)=0;

    };

  }
}

#endif
