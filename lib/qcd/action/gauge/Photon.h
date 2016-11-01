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
  
  template<class Gimpl>
  class Photon
  {
  public:
    INHERIT_GIMPL_TYPES(Gimpl);
    enum class Gauge    {Feynman, Coulomb, Landau};
    enum class ZmScheme {QedL, QedTL};
  public:
    Photon(Gauge gauge, ZmScheme zmScheme);
    virtual ~Photon(void) = default;
    void FreePropagator(const GaugeField &in, GaugeField &out);
    void MomentumSpacePropagator(const GaugeField &in, GaugeField &out);
    void StochasticField(GaugeField &out, GridParallelRNG &rng);
  private:
    void invKHatSquared(GaugeLinkField &out);
    void zmSub(GaugeLinkField &out);
  private:
    Gauge    gauge_;
    ZmScheme zmScheme_;
  };

  template<class Gimpl>
  Photon<Gimpl>::Photon(Gauge gauge, ZmScheme zmScheme)
  : gauge_(gauge), zmScheme_(zmScheme)
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
    
    switch (zmScheme_)
    {
      case ZmScheme::QedTL:
      {
        std::vector<int> zm(nd,0);
        TComplex         Tzero = Complex(0.0,0.0);
        
        pokeSite(Tzero, out, zm);
        
        break;
      }
      case ZmScheme::QedL:
      {
        LatticeInteger spNrm(grid), coor(grid);
        GaugeLinkField z(grid);
        
        spNrm = zero;
        for(int d = 0; d < grid->_ndimension - 1; d++)
        {
          LatticeCoordinate(coor,d);
          spNrm = spNrm + coor*coor;
        }
        out = where(spNrm == Integer(0), 0.*out, out);
        
        break;
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
    LatticeComplex     k2Inv(grid);
    
    invKHatSquared(k2Inv);
    zmSub(k2Inv);
    
    out = in*k2Inv;
  }
  
  template<class Gimpl>
  void Photon<Gimpl>::StochasticField(GaugeField &out, GridParallelRNG &rng)
  {
    auto               *grid = dynamic_cast<GridCartesian *>(out._grid);
    const unsigned int nd = grid->_ndimension;
    std::vector<int> latt_size   = grid->_fdimensions;
    GaugeLinkField     sqrtK2Inv(grid), r(grid);
    GaugeField         aTilde(grid);
    FFT                fft(grid);
    
    Integer vol = 1;
    for(int d = 0; d < nd; d++)
    {
      vol = vol * latt_size[d];
    }

    invKHatSquared(sqrtK2Inv);
    sqrtK2Inv = sqrt(vol*real(sqrtK2Inv));
    zmSub(sqrtK2Inv);
    for(int mu = 0; mu < nd; mu++)
    {
      gaussian(rng, r);
      r = sqrtK2Inv*r;
      pokeLorentz(aTilde, r, mu);
    }
    fft.FFT_all_dim(out, aTilde, FFT::backward);
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
