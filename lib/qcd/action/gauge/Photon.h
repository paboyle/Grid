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

NAMESPACE_BEGIN(Grid);

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
  GRID_SERIALIZABLE_ENUM(ZmScheme, undef, qedL, 1, qedTL, 2);
public:
  Photon(Gauge gauge, ZmScheme zmScheme);
  virtual ~Photon(void) = default;
  void FreePropagator(const GaugeField &in, GaugeField &out);
  void MomentumSpacePropagator(const GaugeField &in, GaugeField &out);
  void StochasticWeight(GaugeLinkField &weight);
  void StochasticField(GaugeField &out, GridParallelRNG &rng);
  void StochasticField(GaugeField &out, GridParallelRNG &rng,
		       const GaugeLinkField &weight);
private:
  void invKHatSquared(GaugeLinkField &out);
  void zmSub(GaugeLinkField &out);
private:
  Gauge    gauge_;
  ZmScheme zmScheme_;
};

typedef Photon<QedGimplR>  PhotonR;
  
template<class Gimpl>
Photon<Gimpl>::Photon(Gauge gauge, ZmScheme zmScheme)
  : gauge_(gauge), zmScheme_(zmScheme)
{}
  
template<class Gimpl>
void Photon<Gimpl>::FreePropagator (const GaugeField &in,GaugeField &out)
{
  FFT theFFT(in.Grid());
    
  GaugeField in_k(in.Grid());
  GaugeField prop_k(in.Grid());
    
  theFFT.FFT_all_dim(in_k,in,FFT::forward);
  MomentumSpacePropagator(prop_k,in_k);
  theFFT.FFT_all_dim(out,prop_k,FFT::backward);
}
  
template<class Gimpl>
void Photon<Gimpl>::invKHatSquared(GaugeLinkField &out)
{
  GridBase           *grid = out.Grid();
  GaugeLinkField     kmu(grid), one(grid);
  const unsigned int nd    = grid->_ndimension;
  Coordinate   &l    = grid->_fdimensions;
  Coordinate   zm(nd,0);
  TComplex           Tone = Complex(1.0,0.0);
  TComplex           Tzero= Complex(0.0,0.0);
    
  one = Complex(1.0,0.0);
  out = Zero();
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
  GridBase           *grid = out.Grid();
  const unsigned int nd    = grid->_ndimension;
    
  switch (zmScheme_)
    {
    case ZmScheme::qedTL:
      {
        Coordinate zm(nd,0);
        TComplex         Tzero = Complex(0.0,0.0);
        
        pokeSite(Tzero, out, zm);
        
        break;
      }
    case ZmScheme::qedL:
      {
        LatticeInteger spNrm(grid), coor(grid);
        GaugeLinkField z(grid);
        
        spNrm = Zero();
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
  GridBase           *grid = out.Grid();
  LatticeComplex     k2Inv(grid);
    
  invKHatSquared(k2Inv);
  zmSub(k2Inv);
    
  out = in*k2Inv;
}
  
template<class Gimpl>
void Photon<Gimpl>::StochasticWeight(GaugeLinkField &weight)
{
  auto               *grid     = dynamic_cast<GridCartesian *>(weight.Grid());
  const unsigned int nd        = grid->_ndimension;
  Coordinate   latt_size = grid->_fdimensions;
    
  Integer vol = 1;
  for(int d = 0; d < nd; d++)
    {
      vol = vol * latt_size[d];
    }
  invKHatSquared(weight);
  weight = sqrt(vol*real(weight));
  zmSub(weight);
}
  
template<class Gimpl>
void Photon<Gimpl>::StochasticField(GaugeField &out, GridParallelRNG &rng)
{
  auto           *grid = dynamic_cast<GridCartesian *>(out.Grid());
  GaugeLinkField weight(grid);
    
  StochasticWeight(weight);
  StochasticField(out, rng, weight);
}
  
template<class Gimpl>
void Photon<Gimpl>::StochasticField(GaugeField &out, GridParallelRNG &rng,
				    const GaugeLinkField &weight)
{
  auto               *grid = dynamic_cast<GridCartesian *>(out.Grid());
  const unsigned int nd = grid->_ndimension;
  GaugeLinkField     r(grid);
  GaugeField         aTilde(grid);
  FFT                fft(grid);
    
  for(int mu = 0; mu < nd; mu++)
    {
      gaussian(rng, r);
      r = weight*r;
      pokeLorentz(aTilde, r, mu);
    }
  fft.FFT_all_dim(out, aTilde, FFT::backward);
    
  out = real(out);
}
//  template<class Gimpl>
//  void Photon<Gimpl>::FeynmanGaugeMomentumSpacePropagator_L(GaugeField &out,
//                                                            const GaugeField &in)
//  {
//    
//    FeynmanGaugeMomentumSpacePropagator_TL(out,in);
//    
//    GridBase *grid = out.Grid();
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
//    GridBase *grid = out.Grid();
//    int nd = grid->_ndimension;
//    
//    typedef typename GaugeField::vector_type vector_type;
//    typedef typename GaugeField::scalar_type ScalComplex;
//    typedef Lattice<iSinglet<vector_type> > LatComplex;
//    
//    Coordinate latt_size   = grid->_fdimensions;
//    
//    LatComplex denom(grid); denom= Zero();
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
//    Coordinate zero_mode(nd,0);
//    TComplexD Tone = ComplexD(1.0,0.0);
//    TComplexD Tzero= ComplexD(0.0,0.0);
//    
//    pokeSite(Tone,denom,zero_mode);
//    
//    denom= one/denom;
//    
//    pokeSite(Tzero,denom,zero_mode);
//    
//    out = Zero();
//    out = in*denom;
//  };
  
NAMESPACE_END(Grid);
#endif
