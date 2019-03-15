/*************************************************************************************
 
 Grid physics library, www.github.com/paboyle/Grid
 
 Source file: ./lib/qcd/action/gauge/Photon.h
 
Copyright (C) 2015-2018
 
 Author: Peter Boyle <paboyle@ph.ed.ac.uk>
 Author: Antonin Portelli <antonin.portelli@me.com>
 Author: James Harrison <J.Harrison@soton.ac.uk>
 
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
#ifndef QCD_PHOTON_ACTION_H
#define QCD_PHOTON_ACTION_H

namespace Grid{
namespace QCD{

  template <class S>
  class QedGImpl
  {
  public:
    typedef S Simd;
    
    template <typename vtype>
    using iImplGaugeLink  = iScalar<iScalar<iScalar<vtype>>>;
    template <typename vtype>
    using iImplGaugeField = iVector<iScalar<iScalar<vtype>>, Nd>;
    
    typedef iImplGaugeLink<Simd>  SiteLink;
    typedef iImplGaugeField<Simd> SiteField;
    typedef SiteLink              SiteComplex;
    
    typedef Lattice<SiteLink>  LinkField;
    typedef Lattice<SiteField> Field;
    typedef Field              ComplexField;
  };
  
  typedef QedGImpl<vComplex> QedGImplR;
  
  template <class GImpl>
  class Photon
  {
  public:
    INHERIT_GIMPL_TYPES(GImpl);
    typedef typename SiteGaugeLink::scalar_object ScalarSite;
    typedef typename ScalarSite::scalar_type      ScalarComplex;
    GRID_SERIALIZABLE_ENUM(Gauge, undef, feynman, 1, coulomb, 2, landau, 3);
    GRID_SERIALIZABLE_ENUM(ZmScheme, undef, qedL, 1, qedTL, 2);
  public:
    Photon(GridBase *grid, Gauge gauge, ZmScheme zmScheme, std::vector<Real> improvement);
    Photon(GridBase *grid, Gauge gauge, ZmScheme zmScheme);
    virtual ~Photon(void) = default;
    void FreePropagator(const GaugeField &in, GaugeField &out);
    void MomentumSpacePropagator(const GaugeField &in, GaugeField &out);
    void StochasticWeight(GaugeLinkField &weight);
    void StochasticField(GaugeField &out, GridParallelRNG &rng);
    void StochasticField(GaugeField &out, GridParallelRNG &rng,
                         const GaugeLinkField &weight);
    void UnitField(GaugeField &out);
  private:
    void makeSpatialNorm(LatticeInteger &spNrm);
    void makeKHat(std::vector<GaugeLinkField> &khat);
    void makeInvKHatSquared(GaugeLinkField &out);
    void zmSub(GaugeLinkField &out);
    void transverseProjectSpatial(GaugeField &out);
    void gaugeTransform(GaugeField &out);
  private:
    GridBase          *grid_;
    Gauge             gauge_;
    ZmScheme          zmScheme_;
    std::vector<Real> improvement_;
  };

  typedef Photon<QedGImplR>  PhotonR;
  
  template<class GImpl>
  Photon<GImpl>::Photon(GridBase *grid, Gauge gauge, ZmScheme zmScheme,
                        std::vector<Real> improvements)
  : grid_(grid), gauge_(gauge), zmScheme_(zmScheme), improvement_(improvements)
  {}

  template<class GImpl>
  Photon<GImpl>::Photon(GridBase *grid, Gauge gauge, ZmScheme zmScheme)
  : Photon(grid, gauge, zmScheme, std::vector<Real>())
  {}

  template<class GImpl>
  void Photon<GImpl>::FreePropagator(const GaugeField &in, GaugeField &out)
  {
    FFT        theFFT(dynamic_cast<GridCartesian *>(grid_));
    GaugeField in_k(grid_);
    GaugeField prop_k(grid_);
    
    theFFT.FFT_all_dim(in_k, in, FFT::forward);
    MomentumSpacePropagator(prop_k, in_k);
    theFFT.FFT_all_dim(out, prop_k, FFT::backward);
  }

  template<class GImpl>
  void Photon<GImpl>::makeSpatialNorm(LatticeInteger &spNrm)
  {
    LatticeInteger   coor(grid_);
    std::vector<int> l = grid_->FullDimensions();

    spNrm = zero;
    for(int mu = 0; mu < grid_->Nd() - 1; mu++)
    {
      LatticeCoordinate(coor, mu);
      coor  = where(coor < Integer(l[mu]/2), coor, coor - Integer(l[mu]));
      spNrm = spNrm + coor*coor;
    }
  }

  template<class GImpl>
  void Photon<GImpl>::makeKHat(std::vector<GaugeLinkField> &khat)
  {
    const unsigned int nd = grid_->Nd();
    std::vector<int>   l  = grid_->FullDimensions();
    Complex            ci(0., 1.);

    khat.resize(nd, grid_);
    for (unsigned int mu = 0; mu < nd; ++mu)
    {
      Real piL = M_PI/l[mu];

      LatticeCoordinate(khat[mu], mu);
      khat[mu] = exp(piL*ci*khat[mu])*2.*sin(piL*khat[mu]);
    }
  }

  template<class GImpl>
  void Photon<GImpl>::makeInvKHatSquared(GaugeLinkField &out)
  {
    std::vector<GaugeLinkField> khat;
    GaugeLinkField              lone(grid_);
    const unsigned int          nd = grid_->Nd();
    std::vector<int>            zm(nd, 0);
    ScalarSite                  one = ScalarComplex(1., 0.), z = ScalarComplex(0., 0.);
    
    out = zero;
    makeKHat(khat);
    for(int mu = 0; mu < nd; mu++)
    {
      out = out + khat[mu]*conjugate(khat[mu]);
    }
    lone = ScalarComplex(1., 0.);
    pokeSite(one, out, zm);
    out = lone/out;
    pokeSite(z, out, zm);
  }
  
  template<class GImpl>
  void Photon<GImpl>::zmSub(GaugeLinkField &out)
  {
    switch (zmScheme_)
    {
      case ZmScheme::qedTL:
      {
        std::vector<int> zm(grid_->Nd(), 0);
        ScalarSite       z = ScalarComplex(0., 0.);
        
        pokeSite(z, out, zm);
        break;
      }
      case ZmScheme::qedL:
      {
        LatticeInteger spNrm(grid_);

        makeSpatialNorm(spNrm);
        out = where(spNrm == Integer(0), 0.*out, out);
        for(int i = 0; i < improvement_.size(); i++)
        {
          Real f = sqrt(improvement_[i] + 1);
          out = where(spNrm == Integer(i + 1), f*out, out);
        }
        break;
      }
      default:
        assert(0);
        break;
    }
  }

  template<class GImpl>
  void Photon<GImpl>::transverseProjectSpatial(GaugeField &out)
  {
    const unsigned int          nd = grid_->Nd();
    GaugeLinkField              invKHat(grid_), cst(grid_), spdiv(grid_);
    LatticeInteger              spNrm(grid_);
    std::vector<GaugeLinkField> khat, a(nd, grid_), aProj(nd, grid_);

    invKHat = zero;
    makeSpatialNorm(spNrm);
    makeKHat(khat);
    for (unsigned int mu = 0; mu < nd; ++mu)
    {
      a[mu] = peekLorentz(out, mu);
      if (mu < nd - 1)
      {
        invKHat += khat[mu]*conjugate(khat[mu]);
      }
    }
    cst     = ScalarComplex(1., 0.);
    invKHat = where(spNrm == Integer(0), cst, invKHat);
    invKHat = cst/invKHat;
    cst     = zero;
    invKHat = where(spNrm == Integer(0), cst, invKHat);
    spdiv   = zero;
    for (unsigned int nu = 0; nu < nd - 1; ++nu)
    {
      spdiv += conjugate(khat[nu])*a[nu];
    }
    spdiv *= invKHat;
    for (unsigned int mu = 0; mu < nd; ++mu)
    {
      aProj[mu] = a[mu] - khat[mu]*spdiv;
      pokeLorentz(out, aProj[mu], mu);
    }
  }

  template<class GImpl>
  void Photon<GImpl>::gaugeTransform(GaugeField &out)
  {
    switch (gauge_)
    {
      case Gauge::feynman:
        break;
      case Gauge::coulomb:
        transverseProjectSpatial(out);
        break;
      case Gauge::landau:
        assert(0);
        break;
      default:
        assert(0);
        break;
    }
  }

  template<class GImpl>
  void Photon<GImpl>::MomentumSpacePropagator(const GaugeField &in,
                                              GaugeField &out)
  {
    LatticeComplex momProp(grid_);
    
    makeInvKHatSquared(momProp);
    zmSub(momProp);
    
    out = in*momProp;
  }
  
  template<class GImpl>
  void Photon<GImpl>::StochasticWeight(GaugeLinkField &weight)
  {
    const unsigned int nd  = grid_->Nd();
    std::vector<int>   l   = grid_->FullDimensions();
    Integer            vol = 1;

    for(unsigned int mu = 0; mu < nd; mu++)
    {
      vol = vol*l[mu];
    }
    makeInvKHatSquared(weight);
    weight = sqrt(vol)*sqrt(weight);
    zmSub(weight);
  }
  
  template<class GImpl>
  void Photon<GImpl>::StochasticField(GaugeField &out, GridParallelRNG &rng)
  {
    GaugeLinkField weight(grid_);
    
    StochasticWeight(weight);
    StochasticField(out, rng, weight);
  }
  
  template<class GImpl>
  void Photon<GImpl>::StochasticField(GaugeField &out, GridParallelRNG &rng,
                                      const GaugeLinkField &weight)
  {
    const unsigned int nd = grid_->Nd();
    GaugeLinkField     r(grid_);
    GaugeField         aTilde(grid_);
    FFT                fft(dynamic_cast<GridCartesian *>(grid_));
    
    for(unsigned int mu = 0; mu < nd; mu++)
    {
      gaussian(rng, r);
      r = weight*r;
      pokeLorentz(aTilde, r, mu);
    }
    gaugeTransform(aTilde);
    fft.FFT_all_dim(out, aTilde, FFT::backward);
    out = real(out);
  }

  template<class GImpl>
  void Photon<GImpl>::UnitField(GaugeField &out)
  {
    const unsigned int nd = grid_->Nd();
    GaugeLinkField     r(grid_);
    
    r = ScalarComplex(1., 0.);
    for(unsigned int mu = 0; mu < nd; mu++)
    {
      pokeLorentz(out, r, mu);
    }
    out = real(out);
  }
  
}}
#endif
