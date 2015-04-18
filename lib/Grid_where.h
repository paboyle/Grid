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
  typedef typename vobj::scalar_type scalar_type;
  typedef typename vobj::vector_type vector_type;
  typedef typename iobj::vector_type mask_type;

  const int Nsimd = grid->Nsimd();
  const int words = sizeof(vobj)/sizeof(vector_type);

  std::vector<Integer> mask(Nsimd);
  std::vector<std::vector<scalar_type> > truevals (Nsimd,std::vector<scalar_type>(words) );
  std::vector<std::vector<scalar_type> > falsevals(Nsimd,std::vector<scalar_type>(words) );
  std::vector<scalar_type *> pointers(Nsimd);

#pragma omp parallel for
  for(int ss=0;ss<iftrue._grid->oSites(); ss++){

    for(int s=0;s<Nsimd;s++) pointers[s] = & truevals[s][0];
    extract(iftrue._odata[ss]   ,pointers);

    for(int s=0;s<Nsimd;s++) pointers[s] = & falsevals[s][0];
    extract(iffalse._odata[ss]  ,pointers);

    extract(TensorRemove(predicate._odata[ss]),mask);

    for(int s=0;s<Nsimd;s++){
      if (mask[s]) pointers[s]=&truevals[s][0];
      else         pointers[s]=&falsevals[s][0];
    }

    merge(ret._odata[ss],pointers);
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
