#ifndef GRID_LATTICE_REDUCTION_H
#define GRID_LATTICE_REDUCTION_H

namespace Grid {

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Reduction operations
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    template<class vobj>
    inline RealD norm2(const Lattice<vobj> &arg){

      typedef typename vobj::scalar_type scalar;
      typedef typename vobj::vector_type vector;
      decltype(innerProduct(arg._odata[0],arg._odata[0])) vnrm=zero;
      scalar nrm;
      //FIXME make this loop parallelisable
      vnrm=zero;
      for(int ss=0;ss<arg._grid->oSites(); ss++){
	vnrm = vnrm + innerProduct(arg._odata[ss],arg._odata[ss]);
      }
      vector vvnrm =TensorRemove(vnrm) ;
      nrm = Reduce(vvnrm);
      arg._grid->GlobalSum(nrm);
      return real(nrm);
    }

    template<class vobj>
    inline ComplexD innerProduct(const Lattice<vobj> &left,const Lattice<vobj> &right) 
      //    inline auto innerProduct(const Lattice<vobj> &left,const Lattice<vobj> &right) 
      //->decltype(innerProduct(left._odata[0],right._odata[0]))
    {
      typedef typename vobj::scalar_type scalar;
      decltype(innerProduct(left._odata[0],right._odata[0])) vnrm=zero;

      scalar nrm;
      //FIXME make this loop parallelisable
      for(int ss=0;ss<left._grid->oSites(); ss++){
	vnrm = vnrm + innerProduct(left._odata[ss],right._odata[ss]);
      }
      nrm = Reduce(vnrm);
      right._grid->GlobalSum(nrm);
      return nrm;
    }

    template<class vobj>
      inline typename vobj::scalar_object sum(const Lattice<vobj> &arg){

      GridBase *grid=arg._grid;
      int Nsimd = grid->Nsimd();

      typedef typename vobj::scalar_object sobj;
      typedef typename vobj::scalar_type   scalar_type;

      vobj vsum;
      sobj ssum;

      vsum=zero;
      ssum=zero;
      //FIXME make this loop parallelisable
      for(int ss=0;ss<arg._grid->oSites(); ss++){
	vsum = vsum + arg._odata[ss];
      }
      
      std::vector<sobj>               buf(Nsimd);
      std::vector<scalar_type *> pointers(Nsimd);  
      for(int i=0;i<Nsimd;i++) pointers[i] = (scalar_type *)&buf[i];
      extract(vsum,pointers);

      for(int i=0;i<Nsimd;i++) ssum = ssum + buf[i];

      arg._grid->GlobalSum(ssum);

      return ssum;
    }

}
#endif

