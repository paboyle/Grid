#include <Grid/GridCore.h>
#ifndef GRID_UVM

#warning "Using explicit device memory copies"
NAMESPACE_BEGIN(Grid);
#define dprintf 

////////////////////////////////////////////////////////////
// For caching copies of data on device
////////////////////////////////////////////////////////////
const int NaccCacheMax=128; 

typedef struct { 
  void *CpuPtr;
  void *AccPtr;
  size_t bytes;
  uint32_t transient;
  uint32_t state;
  uint32_t accLock;
  uint32_t cpuLock;
} AcceleratorViewEntry;

//////////////////////////////////////////////////////////////////////
// Data tables for ViewCache
//////////////////////////////////////////////////////////////////////
static AcceleratorViewEntry AccCache[NaccCacheMax];
static int AccCacheVictim; // Base for round robin search
static int NaccCache = 32;

////////////////////////////////////
// Priority ordering for unlocked entries
//  Empty
//  CpuDirty 
//  Consistent
//  AccDirty
////////////////////////////////////
#define Empty         (0x0)  /*Entry unoccupied  */
#define CpuDirty      (0x1)  /*CPU copy is golden, Acc buffer MAY not be allocated*/
#define Consistent    (0x2)  /*ACC copy AND CPU copy are valid */
#define AccDirty      (0x4)  /*ACC copy is golden */
#define EvictNext     (0x8)  /*Priority for eviction*/

int   AllocationCache::ViewVictim(void)
{
  int prioEmpty            =-1;
  int prioCpuDirty         =-1;
  int prioConsistent       =-1;
  int prioAccDirty         =-1;
  int prioCpuDirtyEN       =-1;
  int prioConsistentEN     =-1;
  int prioAccDirtyEN       =-1;

  int victim=-1;

  // round robin priority search of unlocked entries offset from current victim
  for(int ep=0;ep<NaccCache;ep++){
    int e = (ep+AccCacheVictim)%NaccCache;
    dprintf("AllocationCacheDeviceMem: Inspecting cache entry %d :",e);

    uint32_t locks = AccCache[e].cpuLock+AccCache[e].accLock;
    uint32_t s = AccCache[e].state;
    uint32_t t = AccCache[e].transient;

    assert( (s==Empty)||(s==CpuDirty)||(s==AccDirty)||(s==Consistent));

    if ( locks==0 ) {

      if( s==Empty       ) { prioEmpty = e; dprintf("Empty"); }

      if( t == EvictNext ) {
	if( s==CpuDirty    ) { prioCpuDirtyEN     = e; dprintf("CpuDirty Transient");}
	if( s==Consistent  ) { prioConsistentEN   = e; dprintf("Consistent Transient");}
	if( s==AccDirty    ) { prioAccDirtyEN     = e; dprintf("AccDirty Transient");}
      } else { 
	if( s==CpuDirty    ) { prioCpuDirty     = e; dprintf("CpuDirty");}
	if( s==Consistent  ) { prioConsistent   = e; dprintf("Consistent");}
	if( s==AccDirty    ) { prioAccDirty     = e; dprintf("AccDirty");}
      } 
      
    } else { 
      if ( AccCache[e].cpuLock ) dprintf("Locked in Cpu ");
      if ( AccCache[e].accLock ) dprintf("Locked in Acc ");
    }
    dprintf("\n");
  }
  // This encodes the prioritisation for device residency
  // EvictNext provides a transient mechanism
  if ( prioAccDirty     >= 0 ) victim = prioAccDirty;
  if ( prioConsistent   >= 0 ) victim = prioConsistent;
  if ( prioCpuDirty     >= 0 ) victim = prioCpuDirty;
  if ( prioAccDirtyEN   >= 0 ) victim = prioAccDirtyEN;
  if ( prioConsistentEN >= 0 ) victim = prioConsistentEN;
  if ( prioCpuDirtyEN   >= 0 ) victim = prioCpuDirtyEN;
  if ( prioEmpty        >= 0 ) victim = prioEmpty;       /*Highest prio is winner*/

  assert(victim >= 0); // Must succeed/
  dprintf("AllocationCacheDeviceMem: Selected victim cache entry %d\n",victim); 

  // advance victim pointer
  AccCacheVictim=(AccCacheVictim+1)%NaccCache;
  dprintf("AllocationCacheDeviceMem: victim pointer now %d / %d\n",AccCacheVictim,NaccCache); 

  return victim;
}
/////////////////////////////////////////////////
// Accelerator cache motion
/////////////////////////////////////////////////

void AllocationCache::Discard(int e) // remove from Accelerator, remove entry, without flush
{
  if(AccCache[e].state!=Empty){
    dprintf("AllocationCache: Discard(%d) %llx,%llx\n",e,(uint64_t)AccCache[e].AccPtr,(uint64_t)AccCache[e].CpuPtr); 
    assert(AccCache[e].accLock==0);
    assert(AccCache[e].cpuLock==0);
    assert(AccCache[e].CpuPtr!=NULL);
    if(AccCache[e].AccPtr) {
      dprintf("AllocationCache: Free(%d) %llx\n",e,(uint64_t)AccCache[e].AccPtr);  
      AcceleratorFree(AccCache[e].AccPtr,AccCache[e].bytes);
    }
  }
  AccCache[e].AccPtr=NULL;
  AccCache[e].CpuPtr=NULL;
  AccCache[e].bytes=0;
  AccCache[e].state=Empty;
  AccCache[e].accLock=0;
  AccCache[e].cpuLock=0;
}

void AllocationCache::Evict(int e) // Make CPU consistent, remove from Accelerator, remove entry
{
  if(AccCache[e].state!=Empty){
    dprintf("AllocationCache: Evict(%d) %llx,%llx\n",e,(uint64_t)AccCache[e].AccPtr,(uint64_t)AccCache[e].CpuPtr); 
    assert(AccCache[e].accLock==0);
    assert(AccCache[e].cpuLock==0);
    if(AccCache[e].state==AccDirty) {
      Flush(e);
    }
    assert(AccCache[e].CpuPtr!=NULL);
    if(AccCache[e].AccPtr) {
      dprintf("AllocationCache: Free(%d) %llx\n",e,(uint64_t)AccCache[e].AccPtr);  
      AcceleratorFree(AccCache[e].AccPtr,AccCache[e].bytes);
    }
  }
  AccCache[e].AccPtr=NULL;
  AccCache[e].CpuPtr=NULL;
  AccCache[e].bytes=0;
  AccCache[e].state=Empty;
  AccCache[e].accLock=0;
  AccCache[e].cpuLock=0;
}
void AllocationCache::Flush(int e)// Copy back from a dirty device state and mark consistent. Do not remove
{
  //  printf("AllocationCache: Flush(%d) %llx -> %llx\n",e,(uint64_t)AccCache[e].AccPtr,(uint64_t)AccCache[e].CpuPtr); fflush(stdout);
  assert(AccCache[e].state==AccDirty);
  assert(AccCache[e].cpuLock==0);
  assert(AccCache[e].accLock==0);
  assert(AccCache[e].AccPtr!=NULL);
  assert(AccCache[e].CpuPtr!=NULL);
  acceleratorCopyFromDevice(AccCache[e].AccPtr,AccCache[e].CpuPtr,AccCache[e].bytes);
  AccCache[e].state=Consistent;
}
void AllocationCache::Clone(int e)// Copy from CPU, mark consistent. Allocate if necessary
{
  assert(AccCache[e].state==CpuDirty);
  assert(AccCache[e].cpuLock==0);
  assert(AccCache[e].accLock==0);
  assert(AccCache[e].CpuPtr!=NULL);
  if(AccCache[e].AccPtr==NULL){
    AccCache[e].AccPtr=AcceleratorAllocate(AccCache[e].bytes);
  }
  //  printf("AllocationCache: Clone(%d) %llx <- %llx\n",e,(uint64_t)AccCache[e].AccPtr,(uint64_t)AccCache[e].CpuPtr); fflush(stdout);
  acceleratorCopyToDevice(AccCache[e].CpuPtr,AccCache[e].AccPtr,AccCache[e].bytes);
  AccCache[e].state=Consistent;
}

void AllocationCache::CpuDiscard(int e)// Mark accelerator dirty without copy. Allocate if necessary
{
  assert(AccCache[e].state!=Empty);
  assert(AccCache[e].cpuLock==0);
  assert(AccCache[e].accLock==0);
  assert(AccCache[e].CpuPtr!=NULL);
  if(AccCache[e].AccPtr==NULL){
    AccCache[e].AccPtr=AcceleratorAllocate(AccCache[e].bytes);
  }
  //  printf("AllocationCache: CpuDiscard(%d) %llx <- %llx\n",e,(uint64_t)AccCache[e].AccPtr,(uint64_t)AccCache[e].CpuPtr); fflush(stdout);
  //  acceleratorCopyToDevice(AccCache[e].CpuPtr,AccCache[e].AccPtr,AccCache[e].bytes);
  AccCache[e].state=AccDirty;
}

/////////////////////////////////////////////////////////////////////////////////
// View management
/////////////////////////////////////////////////////////////////////////////////
void AllocationCache::ViewClose(void* Ptr,ViewMode mode)
{
  if( (mode==AcceleratorRead)||(mode==AcceleratorWrite)||(mode==AcceleratorWriteDiscard) ){
    AcceleratorViewClose(Ptr);
  } else if( (mode==CpuRead)||(mode==CpuWrite)){
    CpuViewClose(Ptr);
  } else { 
    assert(0);
  }
}
void *AllocationCache::ViewOpen(void* CpuPtr,size_t bytes,ViewMode mode,ViewAdvise hint)
{
  if( (mode==AcceleratorRead)||(mode==AcceleratorWrite)||(mode==AcceleratorWriteDiscard) ){
    return AcceleratorViewOpen(CpuPtr,bytes,mode,hint);
  } else if( (mode==CpuRead)||(mode==CpuWrite)){
    return CpuViewOpen(CpuPtr,bytes,mode,hint);
  } else { 
    assert(0);
    return nullptr;
  }
}
void *AllocationCache::AcceleratorViewOpen(void* CpuPtr,size_t bytes,ViewMode mode,ViewAdvise hint)
{
  ////////////////////////////////////////////////////////////////////////////
  // Find if present, otherwise get or force an empty
  ////////////////////////////////////////////////////////////////////////////
  int e=CpuViewLookup(CpuPtr);
  if(e==-1) {
    e = ViewVictim();
    dprintf("AcceleratorViewOpen Victim is %d\n",e); 
    Evict(e); // Does copy back if necessary, frees accelerator pointer if not null, sets to empty
  }

  assert((mode==AcceleratorRead)||(mode==AcceleratorWrite)||(mode==AcceleratorWriteDiscard));
  assert(AccCache[e].cpuLock==0);  // Programming error

  if(AccCache[e].state!=Empty) {
    assert(AccCache[e].CpuPtr == CpuPtr);
    assert(AccCache[e].bytes==bytes);
  }
/*
 *  State transitions and actions
 *
 *  Action  State   StateNext         Flush    Clone
 *
 *  AccRead  Empty   Consistent        -        Y
 *  AccWrite Empty   AccDirty          -        Y
 *  AccRead  CpuDirty Consistent       -        Y
 *  AccWrite CpuDirty AccDirty         -        Y
 *  AccRead  Consistent Consistent     -        - 
 *  AccWrite Consistent AccDirty       -        - 
 *  AccRead  AccDirty   AccDirty       -        - 
 *  AccWrite AccDirty   AccDirty       -        - 
 */
  if(AccCache[e].state==Empty) {
    AccCache[e].CpuPtr = CpuPtr;
    AccCache[e].AccPtr = NULL;
    AccCache[e].bytes  = bytes;
    AccCache[e].state  = CpuDirty;   // Cpu starts primary
    if(mode==AcceleratorWriteDiscard){
      CpuDiscard(e);
      AccCache[e].state  = AccDirty;   // Empty + AcceleratorWrite=> AccDirty
    } else if(mode==AcceleratorWrite){
      Clone(e); 
      AccCache[e].state  = AccDirty;   // Empty + AcceleratorWrite=> AccDirty
    } else {
      Clone(e); 
      AccCache[e].state  = Consistent; // Empty + AccRead => Consistent
    }
    AccCache[e].accLock= 1;
    //    printf("Copied Empy entry %d into device accLock %d\n",e,AccCache[e].accLock);
  } else if(AccCache[e].state==CpuDirty ){
    if(mode==AcceleratorWriteDiscard) {
      CpuDiscard(e);
      AccCache[e].state  = AccDirty;   // CpuDirty + AcceleratorWrite=> AccDirty
    } else if(mode==AcceleratorWrite) {
      Clone(e); 
      AccCache[e].state  = AccDirty;   // CpuDirty + AcceleratorWrite=> AccDirty
    } else {
      Clone(e); 
      AccCache[e].state  = Consistent; // CpuDirty + AccRead => Consistent
    }
    AccCache[e].accLock++;
    //    printf("Copied CpuDirty entry %d into device accLock %d\n",e,AccCache[e].accLock);
  } else if(AccCache[e].state==Consistent) {
    if((mode==AcceleratorWrite)||(mode==AcceleratorWriteDiscard))
      AccCache[e].state  = AccDirty;   // Consistent + AcceleratorWrite=> AccDirty
    else
      AccCache[e].state  = Consistent; // Consistent + AccRead => Consistent
    AccCache[e].accLock++;
    //    printf("Consistent entry %d into device accLock %d\n",e,AccCache[e].accLock);
  } else if(AccCache[e].state==AccDirty) {
    if((mode==AcceleratorWrite)||(mode==AcceleratorWriteDiscard))
      AccCache[e].state  = AccDirty; // AccDirty + AcceleratorWrite=> AccDirty
    else
      AccCache[e].state  = AccDirty; // AccDirty + AccRead => AccDirty
    AccCache[e].accLock++;
    //    printf("AccDirty entry %d into device accLock %d\n",e,AccCache[e].accLock);
  } else {
    assert(0);
  }

  int transient =hint;
  AccCache[e].transient= transient? EvictNext : 0;

  return AccCache[e].AccPtr;
}
/*
 *  Action  State   StateNext         Flush    Clone
 *
 *  CpuRead  Empty   CpuDirty          -        -
 *  CpuWrite Empty   CpuDirty          -        -
 *  CpuRead  CpuDirty CpuDirty         -        -
 *  CpuWrite CpuDirty CpuDirty         -        - 
 *  CpuRead  Consistent Consistent     -        - 
 *  CpuWrite Consistent CpuDirty       -        - 
 *  CpuRead  AccDirty   Consistent     Y        -
 *  CpuWrite AccDirty   CpuDirty       Y        -
 */
////////////////////////////////////
// look up & decrement lock count
////////////////////////////////////
void AllocationCache::AcceleratorViewClose(void* AccPtr)
{
  int e=CpuViewLookup(AccPtr);
  //  printf("AccView close %d lock %d \n",e,AccCache[e].accLock);
  if(e==-1) exit(0);
  if(AccCache[e].cpuLock!=0) exit(0);
  if(AccCache[e].accLock==0) exit(0);
  /*
  assert(e!=-1);
  assert(AccCache[e].cpuLock==0);
  assert(AccCache[e].accLock>0);
  */
  AccCache[e].accLock--;
}
void AllocationCache::CpuViewClose(void* CpuPtr)
{
  int e=CpuViewLookup(CpuPtr);
  assert(e!=-1);
  assert(AccCache[e].cpuLock>0);
  assert(AccCache[e].accLock==0);
  AccCache[e].cpuLock--;
}
void *AllocationCache::CpuViewOpen(void* CpuPtr,size_t bytes,ViewMode mode,ViewAdvise transient)
{
  ////////////////////////////////////////////////////////////////////////////
  // Find if present, otherwise get or force an empty
  ////////////////////////////////////////////////////////////////////////////
  int e=CpuViewLookup(CpuPtr);
  if(e==-1) {
    e = ViewVictim();
    dprintf("CpuViewOpen Victim is %d\n",e); 
    Evict(e); // Does copy back if necessary, frees accelerator pointer if not null, sets to empty
  }

  assert((mode==CpuRead)||(mode==CpuWrite));
  assert(AccCache[e].accLock==0);  // Programming error

  if(AccCache[e].state!=Empty) {
    assert(AccCache[e].CpuPtr == CpuPtr);
    assert(AccCache[e].bytes==bytes);
  }

  if(AccCache[e].state==Empty) {
    AccCache[e].CpuPtr = CpuPtr;
    AccCache[e].AccPtr = NULL;
    AccCache[e].bytes  = bytes;
    AccCache[e].state  = CpuDirty; // Empty + CpuRead/CpuWrite => CpuDirty
    AccCache[e].accLock= 0;
    AccCache[e].cpuLock= 1;
  } else if(AccCache[e].state==CpuDirty ){
    // AccPtr dont care, deferred allocate
    AccCache[e].state = CpuDirty; // CpuDirty +CpuRead/CpuWrite => CpuDirty
    AccCache[e].cpuLock++;
  } else if(AccCache[e].state==Consistent) {
    assert(AccCache[e].AccPtr != NULL);
    if(mode==CpuWrite)
      AccCache[e].state = CpuDirty;   // Consistent +CpuWrite => CpuDirty
    else 
      AccCache[e].state = Consistent; // Consistent +CpuRead  => Consistent
    AccCache[e].cpuLock++;
  } else if(AccCache[e].state==AccDirty) {
    assert(AccCache[e].AccPtr != NULL);
    Flush(e);
    if(mode==CpuWrite) AccCache[e].state = CpuDirty;   // AccDirty +CpuWrite => CpuDirty, Flush
    else            AccCache[e].state = Consistent; // AccDirty +CpuRead  => Consistent, Flush
    AccCache[e].cpuLock++;
  } else {
    assert(0); // should be unreachable
  }

  AccCache[e].transient= transient? EvictNext : 0;

  return AccCache[e].CpuPtr;
}

//////////////////////////////////////////////////////////////////////////////
//loop round robin over entries checking acc pointer
//////////////////////////////////////////////////////////////////////////////
int   AllocationCache::CpuViewLookup(void *CpuPtr)
{
  assert(CpuPtr!=NULL);
  for(int e=0;e<NaccCache;e++){
    if ( (AccCache[e].state!=Empty) && (AccCache[e].CpuPtr==CpuPtr) ) {
      return e;
    }
  }
  return -1;
}


NAMESPACE_END(Grid);

#endif
