#ifndef GRID_LATTICE_REDUCTION_H
#define GRID_LATTICE_REDUCTION_H

namespace Grid {
#ifdef GRID_WARN_SUBOPTIMAL
#warning "Optimisation alert all these reduction loops are NOT threaded "
#endif     

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Reduction operations
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    template<class vobj>
    inline RealD norm2(const Lattice<vobj> &arg){

      typedef typename vobj::scalar_type scalar;
      typedef typename vobj::vector_type vector;
      decltype(innerProduct(arg._odata[0],arg._odata[0])) vnrm; 
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
      decltype(innerProduct(left._odata[0],right._odata[0])) vnrm; 

      scalar nrm;
      //FIXME make this loop parallelisable
      vnrm=zero;
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
      extract(vsum,buf);

      for(int i=0;i<Nsimd;i++) ssum = ssum + buf[i];

      arg._grid->GlobalSum(ssum);

      return ssum;
    }




template<class vobj> inline void sliceSum(const Lattice<vobj> &Data,std::vector<typename vobj::scalar_object> &result,int orthogdim)
{
  typedef typename vobj::scalar_object sobj;

  GridBase  *grid = Data._grid;
  assert(grid!=NULL);
  const int    Nd = grid->_ndimension;
  const int Nsimd = grid->Nsimd();

  assert(orthogdim >= 0);
  assert(orthogdim < Nd);

  int fd=grid->_fdimensions[orthogdim];
  int ld=grid->_ldimensions[orthogdim];
  int rd=grid->_rdimensions[orthogdim];

  sobj szero; szero=zero;

  std::vector<vobj,alignedAllocator<vobj> > lvSum(rd); // will locally sum vectors first
  std::vector<sobj> lsSum(ld,szero); // sum across these down to scalars
  std::vector<sobj> extracted(Nsimd);     // splitting the SIMD

  result.resize(fd); // And then global sum to return the same vector to every node for IO to file
  for(int r=0;r<rd;r++){
    lvSum[r]=zero;
  }

  std::vector<int>  coor(Nd);  

  // sum over reduced dimension planes, breaking out orthog dir
  for(int ss=0;ss<grid->oSites();ss++){
    GridBase::CoorFromIndex(coor,ss,grid->_rdimensions);
    int r = coor[orthogdim];
    lvSum[r]=lvSum[r]+Data._odata[ss];
  }  

  // Sum across simd lanes in the plane, breaking out orthog dir.
  std::vector<int> icoor(Nd);

  for(int rt=0;rt<rd;rt++){

    extract(lvSum[rt],extracted);

    for(int idx=0;idx<Nsimd;idx++){

      grid->iCoorFromIindex(icoor,idx);

      int ldx =rt+icoor[orthogdim]*rd;

      lsSum[ldx]=lsSum[ldx]+extracted[idx];

    }
  }
  
  // sum over nodes.
  sobj gsum;
  for(int t=0;t<fd;t++){
    int pt = t/ld; // processor plane
    int lt = t%ld;
    if ( pt == grid->_processor_coor[orthogdim] ) {
      gsum=lsSum[lt];
    } else {
      gsum=zero;
    }

    grid->GlobalSum(gsum);

    result[t]=gsum;
  }

}


}
#endif

