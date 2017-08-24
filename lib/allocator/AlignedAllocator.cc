


#include <Grid/GridCore.h>

namespace Grid {

int PointerCache::victim;

  PointerCache::PointerCacheEntry PointerCache::Entries[PointerCache::Ncache];

void *PointerCache::Insert(void *ptr,size_t bytes) {

  if (bytes < 4096 ) return ptr;

#ifdef GRID_OMP
  assert(omp_in_parallel()==0);
#endif 

  void * ret = NULL;
  int v = -1;

  for(int e=0;e<Ncache;e++) {
    if ( Entries[e].valid==0 ) {
      v=e; 
      break;
    }
  }

  if ( v==-1 ) {
    v=victim;
    victim = (victim+1)%Ncache;
  }

  if ( Entries[v].valid ) {
    ret = Entries[v].address;
    Entries[v].valid = 0;
    Entries[v].address = NULL;
    Entries[v].bytes = 0;
  }

  Entries[v].address=ptr;
  Entries[v].bytes  =bytes;
  Entries[v].valid  =1;

  return ret;
}

void *PointerCache::Lookup(size_t bytes) {

 if (bytes < 4096 ) return NULL;

#ifdef _OPENMP
  assert(omp_in_parallel()==0);
#endif 

  for(int e=0;e<Ncache;e++){
    if ( Entries[e].valid && ( Entries[e].bytes == bytes ) ) {
      Entries[e].valid = 0;
      return Entries[e].address;
    }
  }
  return NULL;
}

}
