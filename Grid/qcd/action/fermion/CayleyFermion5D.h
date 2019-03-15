    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./lib/qcd/action/fermion/CayleyFermion5D.h

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
#ifndef  GRID_QCD_CAYLEY_FERMION_H
#define  GRID_QCD_CAYLEY_FERMION_H

#include <Grid/qcd/action/fermion/WilsonFermion5D.h>

namespace Grid {

  namespace QCD {

     template<typename T> struct switcheroo   {
       static inline int iscomplex()  { return 0; }

       template<class vec>
       static inline vec mult(vec a, vec b) {
	 return real_mult(a,b);
       }
     };
     template<> struct switcheroo<ComplexD> {
       static inline int iscomplex()  { return 1; }

       template<class vec>
       static inline vec mult(vec a, vec b) {
	 return a*b;
       }
     };
     template<> struct switcheroo<ComplexF> {
       static inline int iscomplex()  { return 1; }
       template<class vec>
       static inline vec mult(vec a, vec b) {
	 return a*b;
       }
     };


    template<class Impl>
    class CayleyFermion5D : public WilsonFermion5D<Impl>
    {
    public:
     INHERIT_IMPL_TYPES(Impl);
    public:

      // override multiply
      virtual RealD  M    (const FermionField &in, FermionField &out);
      virtual RealD  Mdag (const FermionField &in, FermionField &out);

      // half checkerboard operations
      virtual void   Meooe       (const FermionField &in, FermionField &out);
      virtual void   MeooeDag    (const FermionField &in, FermionField &out);
      virtual void   Mooee       (const FermionField &in, FermionField &out);
      virtual void   MooeeDag    (const FermionField &in, FermionField &out);
      virtual void   MooeeInv    (const FermionField &in, FermionField &out);
      virtual void   MooeeInvDag (const FermionField &in, FermionField &out);
      virtual void   Meo5D (const FermionField &psi, FermionField &chi);

      virtual void   M5D   (const FermionField &psi, FermionField &chi);
      virtual void   M5Ddag(const FermionField &psi, FermionField &chi);

      ///////////////////////////////////////////////////////////////
      // Physical surface field utilities
      ///////////////////////////////////////////////////////////////
      virtual void   Dminus(const FermionField &psi, FermionField &chi);
      virtual void   DminusDag(const FermionField &psi, FermionField &chi);
      virtual void ExportPhysicalFermionSolution(const FermionField &solution5d,FermionField &exported4d);
      virtual void ExportPhysicalFermionSource(const FermionField &solution5d, FermionField &exported4d);
      virtual void ImportPhysicalFermionSource(const FermionField &input4d,FermionField &imported5d);
      virtual void ImportUnphysicalFermion(const FermionField &solution5d, FermionField &exported4d);

      ///////////////////////////////////////////////////////////////
      // Support for MADWF tricks
      ///////////////////////////////////////////////////////////////
      RealD Mass(void) { return mass; };
      void  SetMass(RealD _mass) { 
	mass=_mass; 
	SetCoefficientsInternal(_zolo_hi,_gamma,_b,_c);  // Reset coeffs
      } ;
      void  P(const FermionField &psi, FermionField &chi);
      void  Pdag(const FermionField &psi, FermionField &chi);

      /////////////////////////////////////////////////////
      // Instantiate different versions depending on Impl
      /////////////////////////////////////////////////////
      void M5D(const FermionField &psi,
	       const FermionField &phi,
	       FermionField &chi,
	       std::vector<Coeff_t> &lower,
	       std::vector<Coeff_t> &diag,
	       std::vector<Coeff_t> &upper);

      void M5Ddag(const FermionField &psi,
		  const FermionField &phi,
		  FermionField &chi,
		  std::vector<Coeff_t> &lower,
		  std::vector<Coeff_t> &diag,
		  std::vector<Coeff_t> &upper);

      void MooeeInternal(const FermionField &in, FermionField &out,int dag,int inv);
      void MooeeInternalCompute(int dag, int inv, Vector<iSinglet<Simd> > & Matp, Vector<iSinglet<Simd> > & Matm);

      void MooeeInternalAsm(const FermionField &in, FermionField &out,
			    int LLs, int site,
			    Vector<iSinglet<Simd> > &Matp,
			    Vector<iSinglet<Simd> > &Matm);
      void MooeeInternalZAsm(const FermionField &in, FermionField &out,
			    int LLs, int site,
			    Vector<iSinglet<Simd> > &Matp,
			    Vector<iSinglet<Simd> > &Matm);


      virtual void   Instantiatable(void)=0;

      // force terms; five routines; default to Dhop on diagonal
      virtual void MDeriv  (GaugeField &mat,const FermionField &U,const FermionField &V,int dag);
      virtual void MoeDeriv(GaugeField &mat,const FermionField &U,const FermionField &V,int dag);
      virtual void MeoDeriv(GaugeField &mat,const FermionField &U,const FermionField &V,int dag);

      // Efficient support for multigrid coarsening
      virtual void  Mdir (const FermionField &in, FermionField &out,int dir,int disp);

      void   Meooe5D       (const FermionField &in, FermionField &out);
      void   MeooeDag5D    (const FermionField &in, FermionField &out);

      //    protected:
      RealD mass;

      // Save arguments to SetCoefficientsInternal
      std::vector<Coeff_t> _gamma;
      RealD                _zolo_hi;
      RealD                _b;
      RealD                _c;

      // Cayley form Moebius (tanh and zolotarev)
      std::vector<Coeff_t> omega;
      std::vector<Coeff_t> bs;    // S dependent coeffs
      std::vector<Coeff_t> cs;
      std::vector<Coeff_t> as;
      // For preconditioning Cayley form
      std::vector<Coeff_t> bee;
      std::vector<Coeff_t> cee;
      std::vector<Coeff_t> aee;
      std::vector<Coeff_t> beo;
      std::vector<Coeff_t> ceo;
      std::vector<Coeff_t> aeo;
      // LDU factorisation of the eeoo matrix
      std::vector<Coeff_t> lee;
      std::vector<Coeff_t> leem;
      std::vector<Coeff_t> uee;
      std::vector<Coeff_t> ueem;
      std::vector<Coeff_t> dee;

      // Matrices of 5d ee inverse params
      Vector<iSinglet<Simd> >  MatpInv;
      Vector<iSinglet<Simd> >  MatmInv;
      Vector<iSinglet<Simd> >  MatpInvDag;
      Vector<iSinglet<Simd> >  MatmInvDag;

      // Constructors
      CayleyFermion5D(GaugeField &_Umu,
		      GridCartesian         &FiveDimGrid,
		      GridRedBlackCartesian &FiveDimRedBlackGrid,
		      GridCartesian         &FourDimGrid,
		      GridRedBlackCartesian &FourDimRedBlackGrid,
		      RealD _mass,RealD _M5,const ImplParams &p= ImplParams());



     void CayleyReport(void);
     void CayleyZeroCounters(void);

     double M5Dflops;
     double M5Dcalls;
     double M5Dtime;

     double MooeeInvFlops;
     double MooeeInvCalls;
     double MooeeInvTime;

    protected:
      virtual void SetCoefficientsZolotarev(RealD zolohi,Approx::zolotarev_data *zdata,RealD b,RealD c);
      virtual void SetCoefficientsTanh(Approx::zolotarev_data *zdata,RealD b,RealD c);
      virtual void SetCoefficientsInternal(RealD zolo_hi,std::vector<Coeff_t> & gamma,RealD b,RealD c);
    };

  }
}
#define INSTANTIATE_DPERP(A)\
template void CayleyFermion5D< A >::M5D(const FermionField &psi,const FermionField &phi,FermionField &chi,\
					std::vector<Coeff_t> &lower,std::vector<Coeff_t> &diag,std::vector<Coeff_t> &upper); \
template void CayleyFermion5D< A >::M5Ddag(const FermionField &psi,const FermionField &phi,FermionField &chi,\
					   std::vector<Coeff_t> &lower,std::vector<Coeff_t> &diag,std::vector<Coeff_t> &upper); \
template void CayleyFermion5D< A >::MooeeInv    (const FermionField &psi, FermionField &chi); \
template void CayleyFermion5D< A >::MooeeInvDag (const FermionField &psi, FermionField &chi);

#undef  CAYLEY_DPERP_DENSE
#define  CAYLEY_DPERP_CACHE
#undef  CAYLEY_DPERP_LINALG
#define CAYLEY_DPERP_VEC

#endif
