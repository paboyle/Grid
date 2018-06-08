/*************************************************************************************
 
 Grid physics library, www.github.com/paboyle/Grid
 
 Source file: ./lib/qcd/action/gauge/Photon.h
 
 Copyright (C) 2015
 
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
#ifndef QCD_PHOTON_ACTION_H
#define QCD_PHOTON_ACTION_H

namespace Grid{
namespace QCD{
  template <class S>
  class QedGimpl
  {
  public:
    typedef S Simd;
    
    template <typename vtype>
    using iImplGaugeLink  = iScalar<iScalar<iScalar<vtype>>>;
    template <typename vtype>
    using iImplGaugeField = iVector<iScalar<iScalar<vtype>>, Nd>;
    
    typedef iImplGaugeLink<Simd>  SiteLink;
    typedef iImplGaugeField<Simd> SiteField;
    typedef SiteField             SiteComplex;
    
    typedef Lattice<SiteLink>  LinkField;
    typedef Lattice<SiteField> Field;
    typedef Field              ComplexField;
  };
  
  typedef QedGimpl<vComplex> QedGimplR;
  
  template<class Gimpl>
  class Photon
  {
  public:
    INHERIT_GIMPL_TYPES(Gimpl);
    GRID_SERIALIZABLE_ENUM(Gauge, undef, feynman, 1, coulomb, 2, landau, 3);
    GRID_SERIALIZABLE_ENUM(ZmScheme, undef, qedL, 1, qedTL, 2, qedInf, 3);
  public:
    Photon(Gauge gauge, ZmScheme zmScheme);
    Photon(Gauge gauge, ZmScheme zmScheme, std::vector<Real> improvements);
    Photon(Gauge gauge, ZmScheme zmScheme, Real G0);
    Photon(Gauge gauge, ZmScheme zmScheme, std::vector<Real> improvements, Real G0);
    virtual ~Photon(void) = default;
    void FreePropagator(const GaugeField &in, GaugeField &out);
    void MomentumSpacePropagator(const GaugeField &in, GaugeField &out);
    void StochasticWeight(GaugeLinkField &weight);
    void StochasticField(GaugeField &out, GridParallelRNG &rng);
    void StochasticField(GaugeField &out, GridParallelRNG &rng,
                         const GaugeLinkField &weight);
    void UnitField(GaugeField &out);
  private:
    void infVolPropagator(GaugeLinkField &out);
    void invKHatSquared(GaugeLinkField &out);
    void zmSub(GaugeLinkField &out);
  private:
    Gauge    gauge_;
    ZmScheme zmScheme_;
    std::vector<Real>  improvement_;
    Real     G0_;
  };

  typedef Photon<QedGimplR>  PhotonR;
  
  template<class Gimpl>
  Photon<Gimpl>::Photon(Gauge gauge, ZmScheme zmScheme)
  : gauge_(gauge), zmScheme_(zmScheme), improvement_(std::vector<Real>()),
    G0_(0.15493339023106021408483720810737508876916113364521)
  {}

  template<class Gimpl>
  Photon<Gimpl>::Photon(Gauge gauge, ZmScheme zmScheme,
                        std::vector<Real> improvements)
  : gauge_(gauge), zmScheme_(zmScheme), improvement_(improvements),
    G0_(0.15493339023106021408483720810737508876916113364521)
  {}

  template<class Gimpl>
  Photon<Gimpl>::Photon(Gauge gauge, ZmScheme zmScheme, Real G0)
  : gauge_(gauge), zmScheme_(zmScheme), improvement_(std::vector<Real>()), G0_(G0)
  {}

  template<class Gimpl>
  Photon<Gimpl>::Photon(Gauge gauge, ZmScheme zmScheme,
                        std::vector<Real> improvements, Real G0)
  : gauge_(gauge), zmScheme_(zmScheme), improvement_(improvements), G0_(G0)
  {}

  template<class Gimpl>
  void Photon<Gimpl>::FreePropagator (const GaugeField &in,GaugeField &out)
  {
    FFT theFFT(in._grid);
    
    GaugeField in_k(in._grid);
    GaugeField prop_k(in._grid);
    
    theFFT.FFT_all_dim(in_k,in,FFT::forward);
    MomentumSpacePropagator(prop_k,in_k);
    theFFT.FFT_all_dim(out,prop_k,FFT::backward);
  }

  template<class Gimpl>
  void Photon<Gimpl>::infVolPropagator(GaugeLinkField &out)
  {
    auto               *grid = dynamic_cast<GridCartesian *>(out._grid);
    LatticeReal        xmu(grid);
    GaugeLinkField     one(grid);
    const unsigned int nd    = grid->_ndimension;
    std::vector<int>   &l    = grid->_fdimensions;
    std::vector<int>   x0(nd,0);
    TComplex           Tone  = Complex(1.0,0.0);
    TComplex           Tzero = Complex(G0_,0.0);
    FFT                fft(grid);
    
    one = Complex(1.0,0.0);
    out = zero;
    for(int mu = 0; mu < nd; mu++)
    {
      LatticeCoordinate(xmu,mu);
      Real lo2 = l[mu]/2.0;
      xmu = where(xmu < lo2, xmu, xmu-double(l[mu]));
      out = out + toComplex(4*M_PI*M_PI*xmu*xmu);
    }
    pokeSite(Tone, out, x0);
    out = one/out;
    pokeSite(Tzero, out, x0);
    fft.FFT_all_dim(out, out, FFT::forward);
  }
  
  template<class Gimpl>
  void Photon<Gimpl>::invKHatSquared(GaugeLinkField &out)
  {
    GridBase           *grid = out._grid;
    GaugeLinkField     kmu(grid), one(grid);
    const unsigned int nd    = grid->_ndimension;
    std::vector<int>   &l    = grid->_fdimensions;
    std::vector<int>   zm(nd,0);
    TComplex           Tone = Complex(1.0,0.0);
    TComplex           Tzero= Complex(0.0,0.0);
    
    one = Complex(1.0,0.0);
    out = zero;
    for(int mu = 0; mu < nd; mu++)
    {
      Real twoPiL = M_PI*2./l[mu];
      
      LatticeCoordinate(kmu,mu);
      kmu = 2.*sin(.5*twoPiL*kmu);
      out = out + kmu*kmu;
    }
    pokeSite(Tone, out, zm);
    out = one/out;
    pokeSite(Tzero, out, zm);
  }
  
  template<class Gimpl>
  void Photon<Gimpl>::zmSub(GaugeLinkField &out)
  {
    GridBase           *grid = out._grid;
    const unsigned int nd    = grid->_ndimension;
    std::vector<int>   &l    = grid->_fdimensions;
    
    switch (zmScheme_)
    {
      case ZmScheme::qedTL:
      {
        std::vector<int> zm(nd,0);
        TComplex         Tzero = Complex(0.0,0.0);
        
        pokeSite(Tzero, out, zm);
        
        break;
      }
      case ZmScheme::qedL:
      {
        LatticeInteger spNrm(grid), coor(grid);
        GaugeLinkField z(grid);
        
        spNrm = zero;
        for(int d = 0; d < grid->_ndimension - 1; d++)
        {
          LatticeCoordinate(coor,d);
          coor = where(coor < Integer(l[d]/2), coor, coor-Integer(l[d]));
          spNrm = spNrm + coor*coor;
        }
        out = where(spNrm == Integer(0), 0.*out, out);

        // IR improvement
        for(int i = 0; i < improvement_.size(); i++)
        {
          Real f = sqrt(improvement_[i]+1);
          out = where(spNrm == Integer(i+1), f*out, out);
        }
      }
      default:
        break;
    }
  }

  template<class Gimpl>
  void Photon<Gimpl>::MomentumSpacePropagator(const GaugeField &in,
                                               GaugeField &out)
  {
  GridBase           *grid = out._grid;
    LatticeComplex     momProp(grid);
    
    switch (zmScheme_)
    {
      case ZmScheme::qedTL:
      case ZmScheme::qedL:
      {
        invKHatSquared(momProp);
        zmSub(momProp);
        break;
      }
      case ZmScheme::qedInf:
      {
        infVolPropagator(momProp);
        break;
      }
      default:
        break;
    }
    
    out = in*momProp;
  }
  
  template<class Gimpl>
  void Photon<Gimpl>::StochasticWeight(GaugeLinkField &weight)
  {
    auto               *grid     = dynamic_cast<GridCartesian *>(weight._grid);
    const unsigned int nd        = grid->_ndimension;
    std::vector<int>   latt_size = grid->_fdimensions;
    
    switch (zmScheme_)
    {
      case ZmScheme::qedTL:
      case ZmScheme::qedL:
      {
        Integer vol = 1;
        for(int d = 0; d < nd; d++)
        {
          vol = vol * latt_size[d];
        }
        invKHatSquared(weight);
        weight = sqrt(vol)*sqrt(weight);
        zmSub(weight);
        break;
      }
      case ZmScheme::qedInf:
      {
        infVolPropagator(weight);
        weight = sqrt(real(weight));
        break;
      }
      default:
        break;
    }
  }
  
  template<class Gimpl>
  void Photon<Gimpl>::StochasticField(GaugeField &out, GridParallelRNG &rng)
  {
    auto           *grid = dynamic_cast<GridCartesian *>(out._grid);
    GaugeLinkField weight(grid);
    
    StochasticWeight(weight);
    StochasticField(out, rng, weight);
  }
  
  template<class Gimpl>
  void Photon<Gimpl>::StochasticField(GaugeField &out, GridParallelRNG &rng,
                                      const GaugeLinkField &weight)
  {
    auto               *grid = dynamic_cast<GridCartesian *>(out._grid);
    const unsigned int nd = grid->_ndimension;
    GaugeLinkField     r(grid);
    GaugeField         aTilde(grid);
    FFT                fft(grid);
    
    switch (zmScheme_)
    {
      case ZmScheme::qedTL:
      case ZmScheme::qedL:
      {
        for(int mu = 0; mu < nd; mu++)
        {
          gaussian(rng, r);
          r = weight*r;
          pokeLorentz(aTilde, r, mu);
        }
        break;
      }
      case ZmScheme::qedInf:
      {
        Complex                    shift(1., 1.); // This needs to be a GaugeLink element?
        for(int mu = 0; mu < nd; mu++)
        {
          bernoulli(rng, r);
          r = weight*(2.*r - shift);
          pokeLorentz(aTilde, r, mu);
        }
        break;
      }
      default:
        break;
    }

    fft.FFT_all_dim(out, aTilde, FFT::backward);
    
    out = real(out);
  }

  template<class Gimpl>
  void Photon<Gimpl>::UnitField(GaugeField &out)
  {
    auto               *grid = dynamic_cast<GridCartesian *>(out._grid);
    const unsigned int nd = grid->_ndimension;
    GaugeLinkField     r(grid);
    
    r = Complex(1.0,0.0);

    for(int mu = 0; mu < nd; mu++)
    {
      pokeLorentz(out, r, mu);
    }
    
    out = real(out);
  }
//  template<class Gimpl>
//  void Photon<Gimpl>::FeynmanGaugeMomentumSpacePropagator_L(GaugeField &out,
//                                                            const GaugeField &in)
//  {
//    
//    FeynmanGaugeMomentumSpacePropagator_TL(out,in);
//    
//    GridBase *grid = out._grid;
//    LatticeInteger     coor(grid);
//    GaugeField zz(grid); zz=zero;
//    
//    // xyzt
//    for(int d = 0; d < grid->_ndimension-1;d++){
//      LatticeCoordinate(coor,d);
//      out = where(coor==Integer(0),zz,out);
//    }
//  }
//  
//  template<class Gimpl>
//  void Photon<Gimpl>::FeynmanGaugeMomentumSpacePropagator_TL(GaugeField &out,
//                                                             const GaugeField &in)
//  {
//    
//    // what type LatticeComplex
//    GridBase *grid = out._grid;
//    int nd = grid->_ndimension;
//    
//    typedef typename GaugeField::vector_type vector_type;
//    typedef typename GaugeField::scalar_type ScalComplex;
//    typedef Lattice<iSinglet<vector_type> > LatComplex;
//    
//    std::vector<int> latt_size   = grid->_fdimensions;
//    
//    LatComplex denom(grid); denom= zero;
//    LatComplex   one(grid); one = ScalComplex(1.0,0.0);
//    LatComplex   kmu(grid);
//    
//    ScalComplex ci(0.0,1.0);
//    // momphase = n * 2pi / L
//    for(int mu=0;mu<Nd;mu++) {
//      
//      LatticeCoordinate(kmu,mu);
//      
//      RealD TwoPiL =  M_PI * 2.0/ latt_size[mu];
//      
//      kmu = TwoPiL * kmu ;
//      
//      denom = denom + 4.0*sin(kmu*0.5)*sin(kmu*0.5); // Wilson term
//    }
//    std::vector<int> zero_mode(nd,0);
//    TComplexD Tone = ComplexD(1.0,0.0);
//    TComplexD Tzero= ComplexD(0.0,0.0);
//    
//    pokeSite(Tone,denom,zero_mode);
//    
//    denom= one/denom;
//    
//    pokeSite(Tzero,denom,zero_mode);
//    
//    out = zero;
//    out = in*denom;
//  };
  
}}
#endif
