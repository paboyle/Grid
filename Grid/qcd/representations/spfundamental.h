

#ifndef SPFUNDAMENTAL_H
#define SPFUNDAMENTAL_H

NAMESPACE_BEGIN(Grid);

/*
 * This is an helper class for the HMC
 * Empty since HMC updates already the fundamental representation 
 */

template <int ncolour>
class SpFundamentalRep {
public:
  static const int Dimension = ncolour;
  static const int nSp = ncolour/2;
  static const bool isFundamental = true;

  // typdef to be used by the Representations class in HMC to get the
  // types for the higher representation fields
  typedef typename Sp<nSp>::LatticeMatrix LatticeMatrix;
  typedef LatticeGaugeField LatticeField;
  
  explicit SpFundamentalRep(GridBase* grid) {} //do nothing
  void update_representation(const LatticeGaugeField& Uin) {} // do nothing

  LatticeField RtoFundamentalProject(const LatticeField& in, Real scale = 1.0) const{
    return (scale * in);
  }

};


    

  
typedef	 SpFundamentalRep<Nc> SpFundamentalRepresentation;

NAMESPACE_END(Grid);  

#endif
