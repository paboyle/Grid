#ifndef QCD_UTILS_COVARIANT_CSHIFT_H
#define QCD_UTILS_COVARIANT_CSHIFT_H
namespace Grid {
namespace QCD {
////////////////////////////////////////////////////////////////////////
// Low performance implementation of CovariantCshift API
////////////////////////////////////////////////////////////////////////
// Make these members of an Impl class for BC's.
template<class covariant,class gauge> Lattice<covariant> CovShiftForward(const Lattice<gauge> &Link, 
									    int mu,
									    const Lattice<covariant> &field)
{
  return Link*Cshift(field,mu,1);// moves towards negative mu
}
template<class covariant,class gauge> Lattice<covariant> CovShiftBackward(const Lattice<gauge> &Link, 
									     int mu,
									     const Lattice<covariant> &field)
{
  Lattice<covariant> tmp(field._grid);
  tmp = adj(Link)*field;
  return Cshift(tmp,mu,-1);// moves towards positive mu
}

}}
#endif
