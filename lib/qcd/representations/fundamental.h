/*
 *	Policy classes for the HMC
 *	Author: Guido Cossu
*/	

#ifndef FUNDAMENTAL_H
#define FUNDAMENTAL_H


namespace Grid {
namespace QCD {

/*
* This is an helper class for the HMC
* Empty since HMC updates already the fundamental representation 
*/

template <int ncolour>
class FundamentalRep {
 public:
  static const int Dimension = ncolour;
  static const bool isFundamental = true;

  // typdef to be used by the Representations class in HMC to get the
  // types for the higher representation fields
  typedef typename SU<ncolour>::LatticeMatrix LatticeMatrix;
  typedef LatticeGaugeField LatticeField;
  
  explicit FundamentalRep(GridBase* grid) {} //do nothing
  void update_representation(const LatticeGaugeField& Uin) {} // do nothing

  LatticeField RtoFundamentalProject(const LatticeField& in, Real scale = 1.0) const{
    return (scale * in);
  }

};

  template<class Field> 
  class EmptyRep {
  public:
    typedef Field LatticeField;

    explicit EmptyRep(GridBase* grid) {} //do nothing
    void update_representation(const LatticeField& Uin) {} // do nothing
    LatticeField RtoFundamentalProject(const LatticeField& in, Real scale = 1.0) const{}// do nothing
  };
    

  
typedef	 FundamentalRep<Nc> FundamentalRepresentation;

  
}
}




#endif
