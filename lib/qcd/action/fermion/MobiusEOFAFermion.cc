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

    /***************************************************************
    /* Additional EOFA operators only called outside the inverter.
    /* Since speed is not essential, simple axpby-style
    /* implementations should be fine.
    /***************************************************************/
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
    void MobiusEOFAFermion<Impl>::DtildeInv(const FermionField& psi, FermionField& chi){ }

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
    /* Performance critical fermion operators called inside the inverter
    /********************************************************************/

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
        int   Ls    = this->Ls;
        int   pm    = this->pm;
        RealD shift = this->shift;
        RealD mq1   = this->mq1;
        RealD mq2   = this->mq2;
        RealD mq3   = this->mq3;

        // coefficients for shift operator ( = shift*\gamma_{5}*R_{5}*\Delta_{\pm}(mq2,mq3)*P_{\pm} )
        Coeff_t shiftp(0.0), shiftm(0.0);
        if(shift != 0.0){
          if(pm == 1){ shiftp = shift*(mq3-mq2); }
          else{ shiftm = -shift*(mq3-mq2); }
        }

        std::vector<Coeff_t> diag(Ls,1.0);
        std::vector<Coeff_t> upper(Ls,-1.0); upper[Ls-1] = mq1 + shiftp;
        std::vector<Coeff_t> lower(Ls,-1.0); lower[0]    = mq1 + shiftm;

        #if(0)
            std::cout << GridLogMessage << "MobiusEOFAFermion::M5Ddag(FF&,FF&):" << std::endl;
            for(int i=0; i<diag.size(); ++i){
                std::cout << GridLogMessage << "diag[" << i << "] =" << diag[i] << std::endl;
            }
            for(int i=0; i<upper.size(); ++i){
                std::cout << GridLogMessage << "upper[" << i << "] =" << upper[i] << std::endl;
            }
            for(int i=0; i<lower.size(); ++i){
                std::cout << GridLogMessage << "lower[" << i << "] =" << lower[i] << std::endl;
            }
        #endif

        this->M5Ddag(psi, chi, chi, lower, diag, upper);
    }

    // half checkerboard operations
    template<class Impl>
    void MobiusEOFAFermion<Impl>::Mooee(const FermionField& psi, FermionField& chi)
    {
        int Ls = this->Ls;

        std::vector<Coeff_t> diag = this->bee;
        std::vector<Coeff_t> upper(Ls);
        std::vector<Coeff_t> lower(Ls);

        for(int s=0; s<Ls; s++){
          upper[s] = -this->cee[s];
          lower[s] = -this->cee[s];
        }
        upper[Ls-1] = this->dm;
        lower[0]    = this->dp;

        this->M5D(psi, psi, chi, lower, diag, upper);
    }

    template<class Impl>
    void MobiusEOFAFermion<Impl>::MooeeDag(const FermionField& psi, FermionField& chi)
    {
        int Ls = this->Ls;

        std::vector<Coeff_t> diag = this->bee;
        std::vector<Coeff_t> upper(Ls);
        std::vector<Coeff_t> lower(Ls);

        for(int s=0; s<Ls; s++){
          upper[s] = -this->cee[s];
          lower[s] = -this->cee[s];
        }
        upper[Ls-1] = this->dp;
        lower[0]    = this->dm;

        this->M5Ddag(psi, psi, chi, lower, diag, upper);
    }

    /****************************************************************************************/

    //Zolo
    template<class Impl>
    void MobiusEOFAFermion<Impl>::SetCoefficientsInternal(RealD zolo_hi, std::vector<Coeff_t>& gamma, RealD b, RealD c)
    {
        int   Ls    = this->Ls;
        int   pm    = this->pm;
        RealD mq1   = this->mq1;
        RealD mq2   = this->mq2;
        RealD mq3   = this->mq3;
        RealD shift = this->shift;

        ////////////////////////////////////////////////////////
        // Constants for the preconditioned matrix Cayley form
        ////////////////////////////////////////////////////////
        this->bs.resize(Ls);
        this->cs.resize(Ls);
        this->aee.resize(Ls);
        this->aeo.resize(Ls);
        this->bee.resize(Ls);
        this->beo.resize(Ls);
        this->cee.resize(Ls);
        this->ceo.resize(Ls);

        for(int i=0; i<Ls; ++i){
          this->bee[i] = 4.0 - this->M5 + 1.0;
          this->cee[i] = 1.0;
        }

        for(int i=0; i<Ls; ++i){
          this->aee[i] = this->cee[i];
          this->bs[i] = this->beo[i] = 1.0;
          this->cs[i] = this->ceo[i] = 0.0;
        }

        //////////////////////////////////////////
        // EOFA shift terms
        //////////////////////////////////////////
        if(pm == 1){
          this->dp = mq1*this->cee[0] + shift*(mq3-mq2);
          this->dm = mq1*this->cee[Ls-1];
        } else if(this->pm == -1) {
          this->dp = mq1*this->cee[0];
          this->dm = mq1*this->cee[Ls-1] - shift*(mq3-mq2);
        } else {
          this->dp = mq1*this->cee[0];
          this->dm = mq1*this->cee[Ls-1];
        }

        //////////////////////////////////////////
        // LDU decomposition of eeoo
        //////////////////////////////////////////
        this->dee.resize(Ls+1);
        this->lee.resize(Ls);
        this->leem.resize(Ls);
        this->uee.resize(Ls);
        this->ueem.resize(Ls);

        for(int i=0; i<Ls; ++i){

          if(i < Ls-1){

            this->lee[i] = -this->cee[i+1]/this->bee[i]; // sub-diag entry on the ith column

            this->leem[i] = this->dm/this->bee[i];
            for(int j=0; j<i; j++){ this->leem[i] *= this->aee[j]/this->bee[j]; }

            this->dee[i] = this->bee[i];

            this->uee[i] = -this->aee[i]/this->bee[i];   // up-diag entry on the ith row

            this->ueem[i] = this->dp / this->bee[0];
            for(int j=1; j<=i; j++){ this->ueem[i] *= this->cee[j]/this->bee[j]; }

          } else {

            this->lee[i]  = 0.0;
            this->leem[i] = 0.0;
            this->uee[i]  = 0.0;
            this->ueem[i] = 0.0;

          }
        }

        {
          Coeff_t delta_d = 1.0 / this->bee[0];
          for(int j=1; j<Ls-1; j++){ delta_d *= this->cee[j] / this->bee[j]; }
          this->dee[Ls-1] = this->bee[Ls-1] + this->cee[0] * this->dm * delta_d;
          this->dee[Ls] = this->bee[Ls-1] + this->cee[Ls-1] * this->dp * delta_d;
        }

        int inv = 1;
        this->MooeeInternalCompute(0, inv, this->MatpInv, this->MatmInv);
        this->MooeeInternalCompute(1, inv, this->MatpInvDag, this->MatmInvDag);
    }

    // Recompute Cayley-form coefficients for different shift
    template<class Impl>
    void MobiusEOFAFermion<Impl>::RefreshShiftCoefficients(RealD new_shift)
    {
        this->shift = new_shift;
        Approx::zolotarev_data *zdata = Approx::higham(1.0, this->Ls);
        this->SetCoefficientsTanh(zdata, 1.0, 0.0);
    }

    template<class Impl>
    void MobiusEOFAFermion<Impl>::MooeeInternalCompute(int dag, int inv,
        Vector<iSinglet<Simd> >& Matp, Vector<iSinglet<Simd> >& Matm)
    {
        int Ls = this->Ls;

        GridBase* grid = this->FermionRedBlackGrid();
        int LLs = grid->_rdimensions[0];

        if(LLs == Ls){
            return; // Not vectorised in 5th direction
        }

        Eigen::MatrixXcd Pplus  = Eigen::MatrixXcd::Zero(Ls,Ls);
        Eigen::MatrixXcd Pminus = Eigen::MatrixXcd::Zero(Ls,Ls);

        for(int s=0; s<Ls; s++){
            Pplus(s,s)  = this->bee[s];
            Pminus(s,s) = this->bee[s];
        }

        for(int s=0; s<Ls-1; s++){
            Pminus(s,s+1) = -this->cee[s];
        }

        for(int s=0; s<Ls-1; s++){
            Pplus(s+1,s) = -this->cee[s+1];
        }

        Pplus (0,Ls-1) = this->dp;
        Pminus(Ls-1,0) = this->dm;

        Eigen::MatrixXcd PplusMat ;
        Eigen::MatrixXcd PminusMat;

        #if(0)
            std::cout << GridLogMessage << "Pplus:" << std::endl;
            for(int s=0; s<Ls; ++s){
                for(int ss=0; ss<Ls; ++ss){
                    std::cout << Pplus(s,ss) << "\t";
                }
                std::cout << std::endl;
            }
            std::cout << GridLogMessage << "Pminus:" << std::endl;
            for(int s=0; s<Ls; ++s){
                for(int ss=0; ss<Ls; ++ss){
                    std::cout << Pminus(s,ss) << "\t";
                }
                std::cout << std::endl;
            }
        #endif

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
