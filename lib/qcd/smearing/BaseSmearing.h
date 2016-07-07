/*
  @brief Declares base smearing class Smear
 */
#ifndef BASE_SMEAR_
#define BASE_SMEAR_

template <class Gimpl> 
class Smear{
public:
  INHERIT_GIMPL_TYPES(Gimpl) // inherits the types for the gauge fields

  virtual ~Smear(){}
  virtual void smear     (GaugeField&,const GaugeField&)const = 0;
  virtual void derivative(GaugeField&,
			  const GaugeField&,const GaugeField&) const = 0;
};
#endif
