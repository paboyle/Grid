#ifndef GRID_EXTRACT_H
#define GRID_EXTRACT_H
/////////////////////////////////////////////////////////////////
// Generic extract/merge/permute
/////////////////////////////////////////////////////////////////

namespace Grid{

////////////////////////////////////////////////////////////////////////////////////////////////
// Extract/merge a fundamental vector type, to pointer array with offset
////////////////////////////////////////////////////////////////////////////////////////////////

template<class vsimd,class scalar>
inline void extract(typename std::enable_if<isGridTensor<vsimd>::notvalue, const vsimd >::type * y, 
		    std::vector<scalar *> &extracted,int offset){
  // FIXME: bounce off memory is painful
  int Nextr=extracted.size();
  int Nsimd=vsimd::Nsimd();
  int s=Nsimd/Nextr;

  scalar*buf = (scalar *)y;
  for(int i=0;i<Nextr;i++){
    extracted[i][offset] = buf[i*s];
  }
};
////////////////////////////////////////////////////////////////////////
// Merge simd vector from array of scalars to pointer array with offset
////////////////////////////////////////////////////////////////////////
template<class vsimd,class scalar>
inline void merge(typename std::enable_if<isGridTensor<vsimd>::notvalue, vsimd >::type * y, 
		  std::vector<scalar *> &extracted,int offset){
  int Nextr=extracted.size();
  int Nsimd=vsimd::Nsimd();
  int s=Nsimd/Nextr;

  scalar *buf =(scalar *) y;
  for(int i=0;i<Nextr;i++){
    for(int ii=0;ii<s;ii++){
      buf[i*s+ii]=extracted[i][offset];
    }
  }

};


////////////////////////////////////////////////////////////////////////////////////////////////
// Extract a fundamental vector type to scalar array 
////////////////////////////////////////////////////////////////////////////////////////////////
template<class vsimd,class scalar>
inline void extract(typename std::enable_if<isGridTensor<vsimd>::notvalue, const vsimd >::type  &y,std::vector<scalar> &extracted){

  int Nextr=extracted.size();
  int Nsimd=vsimd::Nsimd();
  int s=Nsimd/Nextr;

  scalar *buf = (scalar *)&y;
  for(int i=0;i<Nextr;i++){
    for(int ii=0;ii<s;ii++){
      extracted[i]=buf[i*s+ii];
    }
  }

};

////////////////////////////////////////////////////////////////////////
// Merge simd vector from array of scalars
////////////////////////////////////////////////////////////////////////
template<class vsimd,class scalar>
inline void merge(typename std::enable_if<isGridTensor<vsimd>::notvalue, vsimd >::type  &y,std::vector<scalar> &extracted){
  int Nextr=extracted.size();
  int Nsimd=vsimd::Nsimd();
  int s=Nsimd/Nextr;
  scalar *buf = (scalar *)&y;

  for(int i=0;i<Nextr;i++){
    for(int ii=0;ii<s;ii++){
      buf[i*s+ii]=extracted[i];
    }
  }

};
template<class vsimd,class scalar>
inline void AmergeA(typename std::enable_if<isGridTensor<vsimd>::notvalue, vsimd >::type  &y,std::vector<scalar> &extracted){
  int Nextr=extracted.size();
  int Nsimd=vsimd::Nsimd();
  int s=Nsimd/Nextr;

  scalar *buf = (scalar *)&y;
  for(int i=0;i<Nextr;i++){
    for(int ii=0;ii<s;ii++){
      buf[i*s+ii]=extracted[i];
    }
  }
};

////////////////////////////////////////////////////////////////////////
// Extract to contiguous array scalar object
////////////////////////////////////////////////////////////////////////
template<class vobj> inline void extract(const vobj &vec,std::vector<typename vobj::scalar_object> &extracted)
{
  typedef typename vobj::scalar_type scalar_type ;
  typedef typename vobj::vector_type vector_type ;

  const int Nsimd=vobj::vector_type::Nsimd();
  const int words=sizeof(vobj)/sizeof(vector_type);

  extracted.resize(Nsimd);

  std::vector<scalar_type *> pointers(Nsimd);
  for(int i=0;i<Nsimd;i++) 
    pointers[i] =(scalar_type *)& extracted[i];

  vector_type *vp = (vector_type *)&vec;
  for(int w=0;w<words;w++){
    extract<vector_type,scalar_type>(&vp[w],pointers,w);
  }
}
////////////////////////////////////////////////////////////////////////
// Extract to a bunch of scalar object pointers, with offset
////////////////////////////////////////////////////////////////////////
template<class vobj> inline 
void extract(const vobj &vec,std::vector<typename vobj::scalar_object *> &extracted, int offset)
{

  typedef typename vobj::scalar_type scalar_type ;
  typedef typename vobj::vector_type vector_type ;

  const int words=sizeof(vobj)/sizeof(vector_type);
  const int Nsimd=vobj::vector_type::Nsimd();

  assert(extracted.size()==Nsimd);

  std::vector<scalar_type *> pointers(Nsimd);
  for(int i=0;i<Nsimd;i++) {
    pointers[i] =(scalar_type *)& extracted[i][offset];
  }

  vector_type *vp = (vector_type *)&vec;
  for(int w=0;w<words;w++){
    extract<vector_type,scalar_type>(&vp[w],pointers,w);
  }
}

////////////////////////////////////////////////////////////////////////
// Merge a contiguous array of scalar objects
////////////////////////////////////////////////////////////////////////
template<class vobj> inline 
void merge(vobj &vec,std::vector<typename vobj::scalar_object> &extracted)
{
  typedef typename vobj::scalar_type scalar_type ;
  typedef typename vobj::vector_type vector_type ;
  
  const int Nsimd=vobj::vector_type::Nsimd();
  const int words=sizeof(vobj)/sizeof(vector_type);

  assert(extracted.size()==Nsimd);

  std::vector<scalar_type *> pointers(Nsimd);
  for(int i=0;i<Nsimd;i++) 
    pointers[i] =(scalar_type *)& extracted[i];
  
  vector_type *vp = (vector_type *)&vec;
  for(int w=0;w<words;w++){
    merge<vector_type,scalar_type>(&vp[w],pointers,w);
  }
}

////////////////////////////////////////////////////////////////////////
// Merge a bunch of different scalar object pointers, with offset
////////////////////////////////////////////////////////////////////////
template<class vobj> inline 
void merge(vobj &vec,std::vector<typename vobj::scalar_object *> &extracted,int offset)
{
  typedef typename vobj::scalar_type scalar_type ;
  typedef typename vobj::vector_type vector_type ;
  
  const int Nsimd=vobj::vector_type::Nsimd();
  const int words=sizeof(vobj)/sizeof(vector_type);

  assert(extracted.size()==Nsimd);

  std::vector<scalar_type *> pointers(Nsimd);
  for(int i=0;i<Nsimd;i++) 
    pointers[i] =(scalar_type *)& extracted[i][offset];
  
  vector_type *vp = (vector_type *)&vec;
  assert((void *)vp!=NULL);
  for(int w=0;w<words;w++){
    merge<vector_type,scalar_type>(&vp[w],pointers,w);
  }
}
}
#endif
