#ifndef QCD_UTILS_COVARIANT_CSHIFT_H
#define QCD_UTILS_COVARIANT_CSHIFT_H

namespace Grid {
namespace QCD {
////////////////////////////////////////////////////////////////////////
// Low performance implementation of CovariantCshift API
////////////////////////////////////////////////////////////////////////
// Make these members of an Impl class for BC's.

namespace PeriodicBC { 

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
}


namespace GparityBC { 

  // Must give right answers across boundary
  //     <----
  //       --    
  //      |  |              
  // xxxxxxxxxxxxxxxxxxxx
  //      |  | 
  //  
  //     Stap= Cshift(GImpl::CovShiftForward(U[nu],nu,
  //	    	      GImpl::CovShiftForward(U[nu],nu,
  //                  GImpl::CovShiftBackward(U[mu],mu,
  //                  GImpl::CovShiftBackward(U[nu],nu,
  //		      GImpl::CovShiftIdentityBackward(U[nu],nu,-1))))) , mu, 1);
  //  
  //                      U  U^* U^* U^T U^adj =  U  (U U U^dag U^T )^*
  //                                           =  U  (U U U^dag)^* ( U^T )^*
  //
  // So covariant shift rule: conjugate inward shifted plane when crossing boundary applies.
  //
  // This conjugate should be applied to BOTH the link and the covariant field on backward shift
  // boundary wrap.
  // 
  //      |  |              
  // xxxxxxxxxxxxxxxxx         
  //      |  | <---- this link is conjugated, and the path leading into it. Segment crossing in and out is double conjugated.
  //       -- 
  //    ------->
  template<class covariant,class gauge> Lattice<covariant> CovShiftForward(const Lattice<gauge> &Link, 
									    int mu,
									    const Lattice<covariant> &field)
  {
    GridBase * grid = Link._grid;

    int Lmu = grid->GlobalDimensions()[mu]-1;

    conformable(field,Link);

    Lattice<iScalar<vInteger> > coor(grid);    LatticeCoordinate(coor,mu);

    Lattice<covariant> field_bc = Cshift(field,mu,1);// moves towards negative mu;

    field_bc = where(coor==Lmu,conjugate(field_bc),field_bc);
    //    std::cout<<"Gparity::CovCshiftForward mu="<<mu<<std::endl;
    return Link*field_bc;
  }

  template<class covariant,class gauge> Lattice<covariant> CovShiftBackward(const Lattice<gauge> &Link, 
									    int mu,
									    const Lattice<covariant> &field)
  {
    GridBase * grid = field._grid;

    int Lmu = grid->GlobalDimensions()[mu]-1;

    conformable(field,Link);

    Lattice<iScalar<vInteger> > coor(grid);    LatticeCoordinate(coor,mu);

    Lattice<covariant> tmp(grid);

    tmp = adj(Link)*field;
    tmp = where(coor==Lmu,conjugate(tmp),tmp);
    //    std::cout<<"Gparity::CovCshiftBackward mu="<<mu<<std::endl;
    return Cshift(tmp,mu,-1);// moves towards positive mu
  }


}


}}
#endif
