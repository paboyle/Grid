#ifndef GRID_WHERE_H
#define GRID_WHERE_H
namespace Grid {
// Must implement the predicate gating the 
// Must be able to reduce the predicate down to a single vInteger per site.
// Must be able to require the type be iScalar x iScalar x ....
//                              give a GetVtype method in iScalar
//                              and blow away the tensor structures.
//
template<class vobj,class iobj>
inline void where(Lattice<vobj> &ret,const Lattice<iobj> &predicate,Lattice<vobj> &iftrue,Lattice<vobj> &iffalse)
{
  conformable(iftrue,iffalse);
  conformable(iftrue,predicate);
  conformable(iftrue,ret);

  GridBase *grid=iftrue._grid;
  typedef typename vobj::scalar_object scalar_object;
  typedef typename vobj::scalar_type scalar_type;
  typedef typename vobj::vector_type vector_type;
  typedef typename iobj::vector_type mask_type;

  const int Nsimd = grid->Nsimd();
  const int words = sizeof(vobj)/sizeof(vector_type);

  std::vector<Integer> mask(Nsimd);
  std::vector<scalar_object> truevals (Nsimd);
  std::vector<scalar_object> falsevals(Nsimd);

#pragma omp parallel for
  for(int ss=0;ss<iftrue._grid->oSites(); ss++){

    extract(iftrue._odata[ss]   ,truevals);
    extract(iffalse._odata[ss]  ,falsevals);
    extract<vInteger,Integer>(TensorRemove(predicate._odata[ss]),mask);

    for(int s=0;s<Nsimd;s++){
      if (mask[s]) falsevals[s]=truevals[s];
    }

    merge(ret._odata[ss],falsevals);
  }
}

template<class vobj,class iobj>
inline Lattice<vobj> where(const Lattice<iobj> &predicate,Lattice<vobj> &iftrue,Lattice<vobj> &iffalse)
{
  conformable(iftrue,iffalse);
  conformable(iftrue,predicate);

  Lattice<vobj> ret(iftrue._grid);

  where(ret,predicate,iftrue,iffalse);

  return ret;
}
}
#endif
