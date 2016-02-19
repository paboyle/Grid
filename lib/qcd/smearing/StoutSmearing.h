/*
  @file stoutSmear.hpp
  @brief Declares Stout smearing class
*/
#ifndef STOUT_SMEAR_
#define STOUT_SMEAR_

/*!  @brief Stout smearing of link variable. */
template <class Gimpl> 
class Smear_Stout: public Smear<Gimpl> {
private:
  const std::valarray<double> d_rho;
  const Smear* SmearBase;

  double func_xi0(double w) const;
public:
  INHERIT_GIMPL_TYPES(Gimpl)
  Smear_Stout(Smear* base):SmearBase(base){}

  /*! Default constructor */
  Smear_Stout():SmearBase(new Smear_APE()){}

  ~Smear_Stout(){}

  void smear(GaugeField&,const GaugeField&) const;
  void BaseSmear(GaugeField&, const GaugeField&) const;
  void derivative(GaugeField&, const GaugeField&, const GaugeField&) const;
  void exponentiate_iQ(GaugeField&, const GaugeField&) const;

};

#endif  
