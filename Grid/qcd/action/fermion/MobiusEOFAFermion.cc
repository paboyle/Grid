/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/action/fermion/MobiusEOFAFermion.cc

Copyright (C) 2017

Author: Peter Boyle <pabobyle@ph.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: Peter Boyle <peterboyle@Peters-MacBook-Pro-2.local>
Author: paboyle <paboyle@ph.ed.ac.uk>
Author: David Murphy <dmurphy@phys.columbia.edu>

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

#include <Grid/Grid_Eigen_Dense.h>
#include <Grid/qcd/action/fermion/FermionCore.h>
#include <Grid/qcd/action/fermion/MobiusEOFAFermion.h>

namespace Grid {
namespace QCD {

  template<class Impl>
    MobiusEOFAFermion<Impl>::MobiusEOFAFermion(
      GaugeField            &_Umu,
      GridCartesian         &FiveDimGrid,
      GridRedBlackCartesian &FiveDimRedBlackGrid,
      GridCartesian         &FourDimGrid,
      GridRedBlackCartesian &FourDimRedBlackGrid,
      RealD _mq1, RealD _mq2, RealD _mq3,
      RealD _shift, int _pm, RealD _M5,
      RealD _b, RealD _c, const ImplParams &p) :
    AbstractEOFAFermion<Impl>(_Umu, FiveDimGrid, FiveDimRedBlackGrid,
        FourDimGrid, FourDimRedBlackGrid, _mq1, _mq2, _mq3,
        _shift, _pm, _M5, _b, _c, p)
    {
      int Ls = this->Ls;

      RealD eps = 1.0;
      Approx::zolotarev_data *zdata = Approx::higham(eps, this->Ls);
      assert(zdata->n == this->Ls);

      std::cout << GridLogMessage << "MobiusEOFAFermion (b=" << _b <<
        ",c=" << _c << ") with Ls=" << Ls << std::endl;
      this->SetCoefficientsTanh(zdata, _b, _c);
      std::cout << GridLogMessage << "EOFA parameters: (mq1=" << _mq1 <<
        ",mq2=" << _mq2 << ",mq3=" << _mq3 << ",shift=" << _shift <<
        ",pm=" << _pm << ")" << std::endl;

      Approx::zolotarev_free(zdata);

      if(_shift != 0.0){
        SetCoefficientsPrecondShiftOps();
      } else {
        Mooee_shift.resize(Ls, 0.0);
        MooeeInv_shift_lc.resize(Ls, 0.0);
        MooeeInv_shift_norm.resize(Ls, 0.0);
        MooeeInvDag_shift_lc.resize(Ls, 0.0);
        MooeeInvDag_shift_norm.resize(Ls, 0.0);
      }
    }

    /****************************************************************
     * Additional EOFA operators only called outside the inverter.  
     * Since speed is not essential, simple axpby-style
     * implementations should be fine.
     ***************************************************************/
    template<class Impl>
    void MobiusEOFAFermion<Impl>::Omega(const FermionField& psi, FermionField& Din, int sign, int dag)
    {
      int Ls = this->Ls;
      RealD alpha = this->alpha;

      Din = zero;
      if((sign == 1) && (dag == 0)) { // \Omega_{+}
        for(int s=0; s<Ls; ++s){
          axpby_ssp(Din, 0.0, psi, 2.0*std::pow(1.0-alpha,Ls-s-1)/std::pow(1.0+alpha,Ls-s), psi, s, 0);
        }
      } else if((sign == -1) && (dag == 0)) { // \Omega_{-}
        for(int s=0; s<Ls; ++s){
          axpby_ssp(Din, 0.0, psi, 2.0*std::pow(1.0-alpha,s)/std::pow(1.0+alpha,s+1), psi, s, 0);
        }
      } else if((sign == 1 ) && (dag == 1)) { // \Omega_{+}^{\dagger}
        for(int sp=0; sp<Ls; ++sp){
          axpby_ssp(Din, 1.0, Din, 2.0*std::pow(1.0-alpha,Ls-sp-1)/std::pow(1.0+alpha,Ls-sp), psi, 0, sp);
        }
      } else if((sign == -1) && (dag == 1)) { // \Omega_{-}^{\dagger}
        for(int sp=0; sp<Ls; ++sp){
          axpby_ssp(Din, 1.0, Din, 2.0*std::pow(1.0-alpha,sp)/std::pow(1.0+alpha,sp+1), psi, 0, sp);
        }
      }
    }

    // This is the operator relating the usual Ddwf to TWQCD's EOFA Dirac operator (arXiv:1706.05843, Eqn. 6).
    // It also relates the preconditioned and unpreconditioned systems described in Appendix B.2.
    template<class Impl>
    void MobiusEOFAFermion<Impl>::Dtilde(const FermionField& psi, FermionField& chi)
    {
      int Ls    = this->Ls;
      RealD b   = 0.5 * ( 1.0 + this->alpha );
      RealD c   = 0.5 * ( 1.0 - this->alpha );
      RealD mq1 = this->mq1;

      for(int s=0; s<Ls; ++s){
        if(s == 0) {
          axpby_ssp_pminus(chi, b, psi, -c, psi, s, s+1);
          axpby_ssp_pplus (chi, 1.0, chi, mq1*c, psi, s, Ls-1);
        } else if(s == (Ls-1)) {
          axpby_ssp_pminus(chi, b, psi, mq1*c, psi, s, 0);
          axpby_ssp_pplus (chi, 1.0, chi, -c, psi, s, s-1);
        } else {
          axpby_ssp_pminus(chi, b, psi, -c, psi, s, s+1);
          axpby_ssp_pplus (chi, 1.0, chi, -c, psi, s, s-1);
        }
      }
    }

    template<class Impl>
    void MobiusEOFAFermion<Impl>::DtildeInv(const FermionField& psi, FermionField& chi)
    {
      int Ls = this->Ls;
      RealD m = this->mq1;
      RealD c = 0.5 * this->alpha;
      RealD d = 0.5;

      RealD DtInv_p(0.0), DtInv_m(0.0);
      RealD N = std::pow(c+d,Ls) + m*std::pow(c-d,Ls);
      FermionField tmp(this->FermionGrid());

      for(int s=0; s<Ls; ++s){
      for(int sp=0; sp<Ls; ++sp){

        DtInv_p = m * std::pow(-1.0,s-sp+1) * std::pow(c-d,Ls+s-sp) / std::pow(c+d,s-sp+1) / N;
        DtInv_p += (s < sp) ? 0.0 : std::pow(-1.0,s-sp) * std::pow(c-d,s-sp) / std::pow(c+d,s-sp+1);
        DtInv_m = m * std::pow(-1.0,sp-s+1) * std::pow(c-d,Ls+sp-s) / std::pow(c+d,sp-s+1) / N;
        DtInv_m += (s > sp) ? 0.0 : std::pow(-1.0,sp-s) * std::pow(c-d,sp-s) / std::pow(c+d,sp-s+1);

        if(sp == 0){
          axpby_ssp_pplus (tmp, 0.0, tmp, DtInv_p, psi, s, sp);
          axpby_ssp_pminus(tmp, 0.0, tmp, DtInv_m, psi, s, sp);
        } else {
          axpby_ssp_pplus (tmp, 1.0, tmp, DtInv_p, psi, s, sp);
          axpby_ssp_pminus(tmp, 1.0, tmp, DtInv_m, psi, s, sp);
        }

      }}
    }

    /*****************************************************************************************************/

    template<class Impl>
    RealD MobiusEOFAFermion<Impl>::M(const FermionField& psi, FermionField& chi)
    {
      int Ls = this->Ls;

      FermionField Din(psi._grid);

      this->Meooe5D(psi, Din);
      this->DW(Din, chi, DaggerNo);
      axpby(chi, 1.0, 1.0, chi, psi);
      this->M5D(psi, chi);
      return(norm2(chi));
    }

    template<class Impl>
    RealD MobiusEOFAFermion<Impl>::Mdag(const FermionField& psi, FermionField& chi)
    {
      int Ls = this->Ls;

      FermionField Din(psi._grid);

      this->DW(psi, Din, DaggerYes);
      this->MeooeDag5D(Din, chi);
      this->M5Ddag(psi, chi);
      axpby(chi, 1.0, 1.0, chi, psi);
      return(norm2(chi));
    }

    /********************************************************************
     * Performance critical fermion operators called inside the inverter
     ********************************************************************/

    template<class Impl>
    void MobiusEOFAFermion<Impl>::M5D(const FermionField& psi, FermionField& chi)
    {
      int Ls = this->Ls;

      std::vector<Coeff_t> diag(Ls,1.0);
      std::vector<Coeff_t> upper(Ls,-1.0);  upper[Ls-1] = this->mq1;
      std::vector<Coeff_t> lower(Ls,-1.0);  lower[0]    = this->mq1;

      // no shift term
      if(this->shift == 0.0){ this->M5D(psi, chi, chi, lower, diag, upper); }

      // fused M + shift operation
      else{ this->M5D_shift(psi, chi, chi, lower, diag, upper, Mooee_shift); }
    }

    template<class Impl>
    void MobiusEOFAFermion<Impl>::M5Ddag(const FermionField& psi, FermionField& chi)
    {
      int Ls = this->Ls;

      std::vector<Coeff_t> diag(Ls,1.0);
      std::vector<Coeff_t> upper(Ls,-1.0);  upper[Ls-1] = this->mq1;
      std::vector<Coeff_t> lower(Ls,-1.0);  lower[0]    = this->mq1;

      // no shift term
      if(this->shift == 0.0){ this->M5Ddag(psi, chi, chi, lower, diag, upper); }

      // fused M + shift operation
      else{ this->M5Ddag_shift(psi, chi, chi, lower, diag, upper, Mooee_shift); }
    }

    // half checkerboard operations
    template<class Impl>
    void MobiusEOFAFermion<Impl>::Mooee(const FermionField& psi, FermionField& chi)
    {
      int Ls = this->Ls;

      // coefficients of Mooee
      std::vector<Coeff_t> diag = this->bee;
      std::vector<Coeff_t> upper(Ls);
      std::vector<Coeff_t> lower(Ls);
      for(int s=0; s<Ls; s++){
        upper[s] = -this->cee[s];
        lower[s] = -this->cee[s];
      }
      upper[Ls-1] *= -this->mq1;
      lower[0]    *= -this->mq1;

      // no shift term
      if(this->shift == 0.0){ this->M5D(psi, psi, chi, lower, diag, upper); }

      // fused M + shift operation
      else { this->M5D_shift(psi, psi, chi, lower, diag, upper, Mooee_shift); }
    }

    template<class Impl>
    void MobiusEOFAFermion<Impl>::MooeeDag(const FermionField& psi, FermionField& chi)
    {
      int Ls = this->Ls;

      // coefficients of MooeeDag
      std::vector<Coeff_t> diag = this->bee;
      std::vector<Coeff_t> upper(Ls);
      std::vector<Coeff_t> lower(Ls);
      for(int s=0; s<Ls; s++){
        if(s==0) {
          upper[s] = -this->cee[s+1];
          lower[s] = this->mq1*this->cee[Ls-1];
        } else if(s==(Ls-1)) {
          upper[s] = this->mq1*this->cee[0];
          lower[s] = -this->cee[s-1];
        } else {
          upper[s] = -this->cee[s+1];
          lower[s] = -this->cee[s-1];
        }
      }

      // no shift term
      if(this->shift == 0.0){ this->M5Ddag(psi, psi, chi, lower, diag, upper); }

      // fused M + shift operation
      else{ this->M5Ddag_shift(psi, psi, chi, lower, diag, upper, Mooee_shift); }
    }

    /****************************************************************************************/

    // Computes coefficients for applying Cayley preconditioned shift operators
    //  (Mooee + \Delta) --> Mooee_shift
    //  (Mooee + \Delta)^{-1} --> MooeeInv_shift_lc, MooeeInv_shift_norm
    //  (Mooee + \Delta)^{-dag} --> MooeeInvDag_shift_lc, MooeeInvDag_shift_norm
    // For the latter two cases, the operation takes the form
    //  [ (Mooee + \Delta)^{-1} \psi ]_{i} = Mooee_{ij} \psi_{j} +
    //      ( MooeeInv_shift_norm )_{i} ( \sum_{j} [ MooeeInv_shift_lc ]_{j} P_{pm} \psi_{j} )
    template<class Impl>
    void MobiusEOFAFermion<Impl>::SetCoefficientsPrecondShiftOps()
    {
      int   Ls    = this->Ls;
      int   pm    = this->pm;
      RealD alpha = this->alpha;
      RealD k     = this->k;
      RealD mq1   = this->mq1;
      RealD shift = this->shift;

      // Initialize
      Mooee_shift.resize(Ls);
      MooeeInv_shift_lc.resize(Ls);
      MooeeInv_shift_norm.resize(Ls);
      MooeeInvDag_shift_lc.resize(Ls);
      MooeeInvDag_shift_norm.resize(Ls);

      // Construct Mooee_shift
      int idx(0);
      Coeff_t N = ( (pm == 1) ? 1.0 : -1.0 ) * (2.0*shift*k) *
                  ( std::pow(alpha+1.0,Ls) + mq1*std::pow(alpha-1.0,Ls) );
      for(int s=0; s<Ls; ++s){
        idx = (pm == 1) ? (s) : (Ls-1-s);
        Mooee_shift[idx] = N * std::pow(-1.0,s) * std::pow(alpha-1.0,s) / std::pow(alpha+1.0,Ls+s+1);
      }

      // Tridiagonal solve for MooeeInvDag_shift_lc
      {
        Coeff_t m(0.0);
        std::vector<Coeff_t> d = Mooee_shift;
        std::vector<Coeff_t> u(Ls,0.0);
        std::vector<Coeff_t> y(Ls,0.0);
        std::vector<Coeff_t> q(Ls,0.0);
        if(pm == 1){ u[0] = 1.0; }
        else{ u[Ls-1] = 1.0; }

        // Tridiagonal matrix algorithm + Sherman-Morrison formula
        //
        // We solve
        //  ( Mooee' + u \otimes v ) MooeeInvDag_shift_lc = Mooee_shift
        // where Mooee' is the tridiagonal part of Mooee_{+}, and
        // u = (1,0,...,0) and v = (0,...,0,mq1*cee[0]) are chosen
        // so that the outer-product u \otimes v gives the (0,Ls-1)
        // entry of Mooee_{+}.
        //
        // We do this as two solves: Mooee'*y = d and Mooee'*q = u,
        // and then construct the solution to the original system
        //  MooeeInvDag_shift_lc = y - <v,y> / ( 1 + <v,q> ) q
        if(pm == 1){
          for(int s=1; s<Ls; ++s){
            m = -this->cee[s] / this->bee[s-1];
            d[s] -= m*d[s-1];
            u[s] -= m*u[s-1];
          }
        }
        y[Ls-1] = d[Ls-1] / this->bee[Ls-1];
        q[Ls-1] = u[Ls-1] / this->bee[Ls-1];
        for(int s=Ls-2; s>=0; --s){
          if(pm == 1){
            y[s] = d[s] / this->bee[s];
            q[s] = u[s] / this->bee[s];
          } else {
            y[s] = ( d[s] + this->cee[s]*y[s+1] ) / this->bee[s];
            q[s] = ( u[s] + this->cee[s]*q[s+1] ) / this->bee[s];
          }
        }

        // Construct MooeeInvDag_shift_lc
        for(int s=0; s<Ls; ++s){
          if(pm == 1){
            MooeeInvDag_shift_lc[s] = y[s] - mq1*this->cee[0]*y[Ls-1] /
              (1.0+mq1*this->cee[0]*q[Ls-1]) * q[s];
          } else {
            MooeeInvDag_shift_lc[s] = y[s] - mq1*this->cee[Ls-1]*y[0] /
              (1.0+mq1*this->cee[Ls-1]*q[0]) * q[s];
          }
        }

        // Compute remaining coefficients
        N = (pm == 1) ? (1.0 + MooeeInvDag_shift_lc[Ls-1]) : (1.0 + MooeeInvDag_shift_lc[0]);
        for(int s=0; s<Ls; ++s){

          // MooeeInv_shift_lc
          if(pm == 1){ MooeeInv_shift_lc[s] = std::pow(this->bee[s],s) * std::pow(this->cee[s],Ls-1-s); }
          else{ MooeeInv_shift_lc[s] = std::pow(this->bee[s],Ls-1-s) * std::pow(this->cee[s],s); }

          // MooeeInv_shift_norm
          MooeeInv_shift_norm[s] = -MooeeInvDag_shift_lc[s] /
            ( std::pow(this->bee[s],Ls) + mq1*std::pow(this->cee[s],Ls) ) / N;

          // MooeeInvDag_shift_norm
          if(pm == 1){ MooeeInvDag_shift_norm[s] = -std::pow(this->bee[s],s) * std::pow(this->cee[s],Ls-1-s) /
            ( std::pow(this->bee[s],Ls) + mq1*std::pow(this->cee[s],Ls) ) / N; }
          else{ MooeeInvDag_shift_norm[s] = -std::pow(this->bee[s],Ls-1-s) * std::pow(this->cee[s],s) /
            ( std::pow(this->bee[s],Ls) + mq1*std::pow(this->cee[s],Ls) ) / N; }
        }
      }
    }

    // Recompute coefficients for a different value of shift constant
    template<class Impl>
    void MobiusEOFAFermion<Impl>::RefreshShiftCoefficients(RealD new_shift)
    {
      this->shift = new_shift;
      if(new_shift != 0.0){
        SetCoefficientsPrecondShiftOps();
      } else {
        int Ls = this->Ls;
        Mooee_shift.resize(Ls,0.0);
        MooeeInv_shift_lc.resize(Ls,0.0);
        MooeeInv_shift_norm.resize(Ls,0.0);
        MooeeInvDag_shift_lc.resize(Ls,0.0);
        MooeeInvDag_shift_norm.resize(Ls,0.0);
      }
    }

    template<class Impl>
    void MobiusEOFAFermion<Impl>::MooeeInternalCompute(int dag, int inv,
      Vector<iSinglet<Simd> >& Matp, Vector<iSinglet<Simd> >& Matm)
    {
      int Ls = this->Ls;

      GridBase* grid = this->FermionRedBlackGrid();
      int LLs = grid->_rdimensions[0];

      if(LLs == Ls){ return; } // Not vectorised in 5th direction

      Eigen::MatrixXcd Pplus  = Eigen::MatrixXcd::Zero(Ls,Ls);
      Eigen::MatrixXcd Pminus = Eigen::MatrixXcd::Zero(Ls,Ls);

      for(int s=0; s<Ls; s++){
        Pplus(s,s)  = this->bee[s];
        Pminus(s,s) = this->bee[s];
      }

      for(int s=0; s<Ls-1; s++){
        Pminus(s,s+1) = -this->cee[s];
        Pplus(s+1,s) = -this->cee[s+1];
      }

      Pplus (0,Ls-1) = this->mq1*this->cee[0];
      Pminus(Ls-1,0) = this->mq1*this->cee[Ls-1];

      if(this->shift != 0.0){
        RealD c = 0.5 * this->alpha;
        RealD d = 0.5;
        RealD N = this->shift * this->k * ( std::pow(c+d,Ls) + this->mq1*std::pow(c-d,Ls) );
        if(this->pm == 1) {
          for(int s=0; s<Ls; ++s){
            Pplus(s,Ls-1) += N * std::pow(-1.0,s) * std::pow(c-d,s) / std::pow(c+d,Ls+s+1);
          }
        } else {
          for(int s=0; s<Ls; ++s){
            Pminus(s,0) += N * std::pow(-1.0,s+1) * std::pow(c-d,Ls-1-s) / std::pow(c+d,2*Ls-s);
          }
        }
      }

      Eigen::MatrixXcd PplusMat ;
      Eigen::MatrixXcd PminusMat;

      if(inv) {
        PplusMat  = Pplus.inverse();
        PminusMat = Pminus.inverse();
      } else {
        PplusMat  = Pplus;
        PminusMat = Pminus;
      }

      if(dag){
        PplusMat.adjointInPlace();
        PminusMat.adjointInPlace();
      }

      typedef typename SiteHalfSpinor::scalar_type scalar_type;
      const int Nsimd = Simd::Nsimd();
      Matp.resize(Ls*LLs);
      Matm.resize(Ls*LLs);

      for(int s2=0; s2<Ls; s2++){
      for(int s1=0; s1<LLs; s1++){
        int istride = LLs;
        int ostride = 1;
        Simd Vp;
        Simd Vm;
        scalar_type *sp = (scalar_type*) &Vp;
        scalar_type *sm = (scalar_type*) &Vm;
        for(int l=0; l<Nsimd; l++){
          if(switcheroo<Coeff_t>::iscomplex()) {
            sp[l] = PplusMat (l*istride+s1*ostride,s2);
            sm[l] = PminusMat(l*istride+s1*ostride,s2);
          } else {
            // if real
            scalar_type tmp;
            tmp = PplusMat (l*istride+s1*ostride,s2);
            sp[l] = scalar_type(tmp.real(),tmp.real());
            tmp = PminusMat(l*istride+s1*ostride,s2);
            sm[l] = scalar_type(tmp.real(),tmp.real());
          }
        }
        Matp[LLs*s2+s1] = Vp;
        Matm[LLs*s2+s1] = Vm;
      }}
  }

  FermOpTemplateInstantiate(MobiusEOFAFermion);
  GparityFermOpTemplateInstantiate(MobiusEOFAFermion);

}}
