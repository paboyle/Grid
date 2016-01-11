    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/tensors/Tensor_extract_merge.h

    Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: neo <cossu@post.kek.jp>
Author: paboyle <paboyle@ph.ed.ac.uk>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

    See the full license in the file "LICENSE" in the top level distribution directory
    *************************************************************************************/
    /*  END LEGAL */
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
inline void extract(typename std::enable_if<!isGridTensor<vsimd>::value, const vsimd >::type * y, 
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
inline void merge(typename std::enable_if<!isGridTensor<vsimd>::value, vsimd >::type * y, 
		  std::vector<scalar *> &extracted,int offset){
  int Nextr=extracted.size();
  int Nsimd=vsimd::Nsimd();
  int s=Nsimd/Nextr; // can have sparse occupation of simd vector if simd_layout does not fill it
                     // replicate n-fold. Use to allow Integer masks to 
                     // predicate floating point of various width assignments and maintain conformable.
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
inline void extract(typename std::enable_if<!isGridTensor<vsimd>::value, const vsimd >::type  &y,std::vector<scalar> &extracted){

  int Nextr=extracted.size();
  int Nsimd=vsimd::Nsimd();
  int s=Nsimd/Nextr;

  scalar *buf = (scalar *)&y;
  for(int i=0;i<Nextr;i++){
    extracted[i]=buf[i*s];
    for(int ii=1;ii<s;ii++){
      if ( buf[i*s]!=buf[i*s+ii] ){
	std::cout<<GridLogMessage << " SIMD extract failure splat = "<<s<<" ii "<<ii<<" " <<Nextr<<" "<< Nsimd<<" "<<std::endl;
	for(int vv=0;vv<Nsimd;vv++) {
	  std::cout<<GridLogMessage<< buf[vv]<<" ";
	}
	std::cout<<GridLogMessage<<std::endl;
	assert(0);
      }
      assert(buf[i*s]==buf[i*s+ii]);
    }
  }

};

////////////////////////////////////////////////////////////////////////
// Merge simd vector from array of scalars
////////////////////////////////////////////////////////////////////////
template<class vsimd,class scalar>
inline void merge(typename std::enable_if<!isGridTensor<vsimd>::value, vsimd >::type  &y,std::vector<scalar> &extracted){
  int Nextr=extracted.size();
  int Nsimd=vsimd::Nsimd();
  int s=Nsimd/Nextr;
  scalar *buf = (scalar *)&y;

  for(int i=0;i<Nextr;i++){
    for(int ii=0;ii<s;ii++){
      buf[i*s+ii]=extracted[i]; // replicates value
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
  int Nextr=extracted.size();
  const int words=sizeof(vobj)/sizeof(vector_type);
  int s=Nsimd/Nextr;

  std::vector<scalar_type *> pointers(Nextr);
  for(int i=0;i<Nextr;i++) 
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

  int Nextr=extracted.size();
  int s = Nsimd/Nextr;
  scalar_type * vp = (scalar_type *)&vec;

  for(int w=0;w<words;w++){
    for(int i=0;i<Nextr;i++){
      scalar_type * pointer = (scalar_type *)& extracted[i][offset];
      pointer[w] = vp[i*s+w*Nsimd];
    }
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

  int Nextr = extracted.size();
  int splat=Nsimd/Nextr;

  std::vector<scalar_type *> pointers(Nextr);
  for(int i=0;i<Nextr;i++) 
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

  int Nextr=extracted.size();
  int s=Nsimd/Nextr;

  scalar_type *pointer;
  scalar_type *vp = (scalar_type *)&vec;

  //  assert( (((uint64_t)vp)&(sizeof(scalar_type)-1)) == 0);

  for(int w=0;w<words;w++){
    for(int i=0;i<Nextr;i++){
      for(int ii=0;ii<s;ii++){
	pointer=(scalar_type *)&extracted[i][offset];
	vp[w*Nsimd+i*s+ii] = pointer[w];
      }
    }
  }
 }

template<class vobj> inline 
void merge1(vobj &vec,std::vector<typename vobj::scalar_object *> &extracted,int offset)
{
  typedef typename vobj::scalar_type scalar_type ;
  typedef typename vobj::vector_type vector_type ;
  
  const int Nsimd=vobj::vector_type::Nsimd();
  const int words=sizeof(vobj)/sizeof(vector_type);

  scalar_type *pointer;
  scalar_type *vp = (scalar_type *)&vec;

  //  assert( (((uint64_t)vp)&(sizeof(scalar_type)-1)) == 0);

  for(int i=0;i<Nsimd;i++){
    pointer=(scalar_type *)&extracted[i][offset];
    for(int w=0;w<words;w++){
      vp[w*Nsimd+i] = pointer[w];
    }
  }
}

template<class vobj> inline 
void merge2(vobj &vec,std::vector<typename vobj::scalar_object *> &extracted,int offset)
{
  typedef typename vobj::scalar_type scalar_type ;
  typedef typename vobj::vector_type vector_type ;
  
  const int Nsimd=vobj::vector_type::Nsimd();
  const int words=sizeof(vobj)/sizeof(vector_type);

  scalar_type *pointer;
  scalar_type *vp = (scalar_type *)&vec;
  //  assert( (((uint64_t)vp)&(sizeof(scalar_type)-1)) == 0);

  for(int w=0;w<words;w++){
    for(int i=0;i<Nsimd;i++){
      pointer=(scalar_type *)&extracted[i][offset];
      vp[w*Nsimd+i] =pointer[w];
    }
  }
}

}

#endif
