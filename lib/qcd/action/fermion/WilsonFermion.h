    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/action/fermion/WilsonFermion.h

    Copyright (C) 2015

Author: Peter Boyle <pabobyle@ph.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: paboyle <paboyle@ph.ed.ac.uk>

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
#ifndef  GRID_QCD_WILSON_FERMION_H
#define  GRID_QCD_WILSON_FERMION_H

namespace Grid {

  namespace QCD {

    class WilsonFermionStatic {
    public:
      static int HandOptDslash; // these are a temporary hack
      static int MortonOrder;
      static const std::vector<int> directions   ;
      static const std::vector<int> displacements;
      static const int npoint=8;
    };

    template<class Impl>
    class WilsonFermion : public WilsonKernels<Impl>, public WilsonFermionStatic
    {
    public:
    INHERIT_IMPL_TYPES(Impl);
    typedef WilsonKernels<Impl> Kernels;

      ///////////////////////////////////////////////////////////////
      // Implement the abstract base
      ///////////////////////////////////////////////////////////////
      GridBase *GaugeGrid(void)              { return _grid ;}
      GridBase *GaugeRedBlackGrid(void)      { return _cbgrid ;}
      GridBase *FermionGrid(void)            { return _grid;}
      GridBase *FermionRedBlackGrid(void)    { return _cbgrid;}

      //////////////////////////////////////////////////////////////////
      // override multiply; cut number routines if pass dagger argument
      // and also make interface more uniformly consistent
      //////////////////////////////////////////////////////////////////
      RealD M(const FermionField &in, FermionField &out);
      RealD Mdag(const FermionField &in, FermionField &out);

      /////////////////////////////////////////////////////////
      // half checkerboard operations
      // could remain virtual so we  can derive Clover from Wilson base
      /////////////////////////////////////////////////////////
      void Meooe(const FermionField &in, FermionField &out) ;
      void MeooeDag(const FermionField &in, FermionField &out) ;

      // allow override for twisted mass and clover
      virtual void Mooee(const FermionField &in, FermionField &out) ;
      virtual void MooeeDag(const FermionField &in, FermionField &out) ;
      virtual void MooeeInv(const FermionField &in, FermionField &out) ;
      virtual void MooeeInvDag(const FermionField &in, FermionField &out) ;

      ////////////////////////
      // Derivative interface
      ////////////////////////
      // Interface calls an internal routine
      void DhopDeriv(GaugeField &mat,const FermionField &U,const FermionField &V,int dag);
      void DhopDerivOE(GaugeField &mat,const FermionField &U,const FermionField &V,int dag);
      void DhopDerivEO(GaugeField &mat,const FermionField &U,const FermionField &V,int dag);


      ///////////////////////////////////////////////////////////////
      // non-hermitian hopping term; half cb or both
      ///////////////////////////////////////////////////////////////
      void Dhop(const FermionField &in, FermionField &out,int dag) ;
      void DhopOE(const FermionField &in, FermionField &out,int dag) ;
      void DhopEO(const FermionField &in, FermionField &out,int dag) ;

      ///////////////////////////////////////////////////////////////
      // Multigrid assistance; force term uses too
      ///////////////////////////////////////////////////////////////
      void Mdir (const FermionField &in, FermionField &out,int dir,int disp) ;
      void DhopDir(const FermionField &in, FermionField &out,int dir,int disp);
      void DhopDirDisp(const FermionField &in, FermionField &out,int dirdisp,int gamma,int dag) ;

      ///////////////////////////////////////////////////////////////
      // Extra methods added by derived
      ///////////////////////////////////////////////////////////////
      void DerivInternal(StencilImpl & st,
			 DoubledGaugeField & U,
			 GaugeField &mat,
			 const FermionField &A,
			 const FermionField &B,
			 int dag);

      void DhopInternal(StencilImpl & st,DoubledGaugeField & U,
			const FermionField &in, FermionField &out,int dag) ;

      void DhopInternalCommsThenCompute(StencilImpl & st,DoubledGaugeField & U,
				    const FermionField &in, FermionField &out,int dag) ;
      void DhopInternalCommsOverlapCompute(StencilImpl & st,DoubledGaugeField & U,
				    const FermionField &in, FermionField &out,int dag) ;


      // Constructor
      WilsonFermion(GaugeField &_Umu,
		    GridCartesian         &Fgrid,
		    GridRedBlackCartesian &Hgrid, 
		    RealD _mass,
		    const ImplParams &p= ImplParams()
		    ) ;

      // DoubleStore impl dependent
      void ImportGauge(const GaugeField &_Umu);

      ///////////////////////////////////////////////////////////////
      // Data members require to support the functionality
      ///////////////////////////////////////////////////////////////

      //    protected:
    public:

      RealD                        mass;

      GridBase                     *    _grid; 
      GridBase                     *  _cbgrid;

      //Defines the stencils for even and odd
      StencilImpl Stencil; 
      StencilImpl StencilEven; 
      StencilImpl StencilOdd; 

      // Copy of the gauge field , with even and odd subsets
      DoubledGaugeField Umu;
      DoubledGaugeField UmuEven;
      DoubledGaugeField UmuOdd;
      
    };

    typedef WilsonFermion<WilsonImplF> WilsonFermionF;
    typedef WilsonFermion<WilsonImplD> WilsonFermionD;

  }
}
#endif
