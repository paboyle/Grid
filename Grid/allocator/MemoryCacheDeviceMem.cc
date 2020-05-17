#include <Grid/GridCore.h>
#ifndef GRID_UNIFIED

#warning "Using explicit device memory copies"
NAMESPACE_BEGIN(Grid);
#define dprintf(...) 

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

#define Write (1)
#define Read  (2)
#define WriteDiscard (3)
//////////////////////////////////////////////////////////////////////
// Data tables for ViewCache
//////////////////////////////////////////////////////////////////////
static AcceleratorViewEntry AccCache[NaccCacheMax];
static int AccCacheVictim; // Base for round robin search
static int NaccCache = 8;

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

      if( s==Empty       ) { prioEmpty = e; dprintf("Empty");}

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
void AllocationCache::Evict(int e) // Make CPU consistent, remove from Accelerator, remove entry
{
  if(AccCache[e].state!=Empty){
    dprintf("AllocationCache: Evict(%d) %llx,%llxn",e,(uint64_t)AccCache[e].AccPtr,(uint64_t)AccCache[e].CpuPtr);
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
  dprintf("AllocationCache: Flush(%d) %llx -> %llx\n",e,(uint64_t)AccCache[e].AccPtr,(uint64_t)AccCache[e].CpuPtr);
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
  dprintf("AllocationCache: Clone(%d) %llx <- %llx\n",e,(uint64_t)AccCache[e].AccPtr,(uint64_t)AccCache[e].CpuPtr);
  acceleratorCopyToDevice(AccCache[e].CpuPtr,AccCache[e].AccPtr,AccCache[e].bytes);
  AccCache[e].state=Consistent;
}
/////////////////////////////////////////////////////////////////////////////////
// View management
/////////////////////////////////////////////////////////////////////////////////
void *AllocationCache::AccViewOpen(void* CpuPtr,size_t bytes,int mode,int transient)
{
  ////////////////////////////////////////////////////////////////////////////
  // Find if present, otherwise get or force an empty
  ////////////////////////////////////////////////////////////////////////////
  int e=CpuViewLookup(CpuPtr);
  if(e==-1) {
    e = ViewVictim();
    Evict(e); // Does copy back if necessary, frees accelerator pointer if not null, sets to empty
  }

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
    Clone(e); 
    if(mode==Write)
      AccCache[e].state  = AccDirty;   // Empty + AccWrite=> AccDirty
    else
      AccCache[e].state  = Consistent; // Empty + AccRead => Consistent
    AccCache[e].accLock= 1;
  } else if(AccCache[e].state&CpuDirty ){
    Clone(e); 
    if(mode==Write)
      AccCache[e].state  = AccDirty;   // CpuDirty + AccWrite=> AccDirty
    else
      AccCache[e].state  = Consistent; // CpuDirty + AccRead => Consistent
    AccCache[e].accLock++;
  } else if(AccCache[e].state&Consistent) {
    if(mode==Write)
      AccCache[e].state  = AccDirty;   // Consistent + AccWrite=> AccDirty
    else
      AccCache[e].state  = Consistent; // Consistent + AccRead => Consistent
    AccCache[e].accLock++;
  } else if(AccCache[e].state&AccDirty) {
    if(mode==Write)
      AccCache[e].state  = AccDirty; // AccDirty + AccWrite=> AccDirty
    else
      AccCache[e].state  = AccDirty; // AccDirty + AccRead => AccDirty
    AccCache[e].accLock++;
  } else {
    assert(0);
  }

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
void AllocationCache::AccViewClose(void* AccPtr)
{
  int e=AccViewLookup(AccPtr);
  assert(e!=-1);
  assert(AccCache[e].cpuLock==0);
  assert(AccCache[e].accLock>0);
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
void *AllocationCache::CpuViewOpen(void* CpuPtr,size_t bytes,int mode,int transient)
{
  ////////////////////////////////////////////////////////////////////////////
  // Find if present, otherwise get or force an empty
  ////////////////////////////////////////////////////////////////////////////
  int e=CpuViewLookup(CpuPtr);
  if(e==-1) {
    e = ViewVictim();
    Evict(e); // Does copy back if necessary, frees accelerator pointer if not null, sets to empty
  }

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
    if(mode==Write)
      AccCache[e].state = CpuDirty;   // Consistent +CpuWrite => CpuDirty
    else 
      AccCache[e].state = Consistent; // Consistent +CpuRead  => Consistent
    AccCache[e].cpuLock++;
  } else if(AccCache[e].state==AccDirty) {
    assert(AccCache[e].AccPtr != NULL);
    Flush(e);
    if(mode==Write) AccCache[e].state = CpuDirty;   // AccDirty +CpuWrite => CpuDirty, Flush
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
int   AllocationCache::AccViewLookup(void *AccPtr)
{
  assert(AccPtr!=NULL);
  for(int e=0;e<NaccCache;e++){
    if ( (AccCache[e].state!=Empty) && (AccCache[e].AccPtr==AccPtr) ) {
      return e;
    }
  }
  return -1;
}


NAMESPACE_END(Grid);

#endif
