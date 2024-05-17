#include <Grid/GridCore.h>
#ifndef GRID_UVM

#warning "Using explicit device memory copies"
NAMESPACE_BEGIN(Grid);

#define MAXLINE 512
static char print_buffer [ MAXLINE ];

#define mprintf(...) snprintf (print_buffer,MAXLINE, __VA_ARGS__ ); std::cout << GridLogMemory << print_buffer;
#define dprintf(...) snprintf (print_buffer,MAXLINE, __VA_ARGS__ ); std::cout << GridLogDebug << print_buffer;
//#define dprintf(...) 


////////////////////////////////////////////////////////////
// For caching copies of data on device
////////////////////////////////////////////////////////////
MemoryManager::AccViewTable_t MemoryManager::AccViewTable;
MemoryManager::LRU_t MemoryManager::LRU;
  
////////////////////////////////////////////////////////
// Footprint tracking
////////////////////////////////////////////////////////
uint64_t  MemoryManager::DeviceBytes;
uint64_t  MemoryManager::DeviceLRUBytes;
uint64_t  MemoryManager::DeviceMaxBytes = 1024*1024*128;
uint64_t  MemoryManager::HostToDeviceBytes;
uint64_t  MemoryManager::DeviceToHostBytes;
uint64_t  MemoryManager::HostToDeviceXfer;
uint64_t  MemoryManager::DeviceToHostXfer;
uint64_t  MemoryManager::DeviceEvictions;
uint64_t  MemoryManager::DeviceDestroy;

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

/////////////////////////////////////////////////
// Mechanics of data table maintenance
/////////////////////////////////////////////////
int   MemoryManager::EntryPresent(uint64_t CpuPtr)
{
  if(AccViewTable.empty()) return 0;

  auto count = AccViewTable.count(CpuPtr);  assert((count==0)||(count==1));
  return count;
}
void  MemoryManager::EntryCreate(uint64_t CpuPtr,size_t bytes,ViewMode mode,ViewAdvise hint)
{
  assert(!EntryPresent(CpuPtr));
  AcceleratorViewEntry AccCache;
  AccCache.CpuPtr = CpuPtr;
  AccCache.AccPtr = (uint64_t)NULL;
  AccCache.bytes  = bytes;
  AccCache.state  = CpuDirty;
  AccCache.LRU_valid=0;
  AccCache.transient=0;
  AccCache.accLock=0;
  AccCache.cpuLock=0;
  AccViewTable[CpuPtr] = AccCache;
}
MemoryManager::AccViewTableIterator MemoryManager::EntryLookup(uint64_t CpuPtr)
{
  assert(EntryPresent(CpuPtr));
  auto AccCacheIterator = AccViewTable.find(CpuPtr);
  assert(AccCacheIterator!=AccViewTable.end());
  return AccCacheIterator;
}
void MemoryManager::EntryErase(uint64_t CpuPtr)
{
  auto AccCache = EntryLookup(CpuPtr);
  AccViewTable.erase(CpuPtr);
}
void  MemoryManager::LRUinsert(AcceleratorViewEntry &AccCache)
{
  assert(AccCache.LRU_valid==0);
  if (AccCache.transient) { 
    LRU.push_back(AccCache.CpuPtr);
    AccCache.LRU_entry = --LRU.end();
  } else {
    LRU.push_front(AccCache.CpuPtr);
    AccCache.LRU_entry = LRU.begin();
  }
  AccCache.LRU_valid = 1;
  DeviceLRUBytes+=AccCache.bytes;
}
void  MemoryManager::LRUremove(AcceleratorViewEntry &AccCache)
{
  assert(AccCache.LRU_valid==1);
  LRU.erase(AccCache.LRU_entry);
  AccCache.LRU_valid = 0;
  DeviceLRUBytes-=AccCache.bytes;
}
/////////////////////////////////////////////////
// Accelerator cache motion & consistency logic
/////////////////////////////////////////////////
void MemoryManager::AccDiscard(AcceleratorViewEntry &AccCache)
{
  ///////////////////////////////////////////////////////////
  // Remove from Accelerator, remove entry, without flush
  // Cannot be locked. If allocated Must be in LRU pool.
  ///////////////////////////////////////////////////////////
  assert(AccCache.state!=Empty);
  
  dprintf("MemoryManager: Discard(%lx) %lx\n",(uint64_t)AccCache.CpuPtr,(uint64_t)AccCache.AccPtr); 
  assert(AccCache.accLock==0);
  assert(AccCache.cpuLock==0);
  assert(AccCache.CpuPtr!=(uint64_t)NULL);
  if(AccCache.AccPtr) {
    AcceleratorFree((void *)AccCache.AccPtr,AccCache.bytes);
    DeviceDestroy++;
    DeviceBytes   -=AccCache.bytes;
    LRUremove(AccCache);
    AccCache.AccPtr=(uint64_t) NULL;
    dprintf("MemoryManager: Free(%lx) LRU %ld Total %ld\n",(uint64_t)AccCache.AccPtr,DeviceLRUBytes,DeviceBytes);  
  }
  uint64_t CpuPtr = AccCache.CpuPtr;
  EntryErase(CpuPtr);
}

void MemoryManager::Evict(AcceleratorViewEntry &AccCache)
{
  ///////////////////////////////////////////////////////////////////////////
  // Make CPU consistent, remove from Accelerator, remove from LRU, LEAVE CPU only entry
  // Cannot be acclocked. If allocated must be in LRU pool.
  //
  // Nov 2022... Felix issue: Allocating two CpuPtrs, can have an entry in LRU-q with CPUlock.
  //                          and require to evict the AccPtr copy. Eviction was a mistake in CpuViewOpen
  //                          but there is a weakness where CpuLock entries are attempted for erase
  //                          Take these OUT LRU queue when CPU locked?
  //                          Cannot take out the table as cpuLock data is important.
  ///////////////////////////////////////////////////////////////////////////
  assert(AccCache.state!=Empty);
  
  mprintf("MemoryManager: Evict CpuPtr %lx AccPtr %lx cpuLock %ld accLock %ld\n",
	  (uint64_t)AccCache.CpuPtr,(uint64_t)AccCache.AccPtr,
	  (uint64_t)AccCache.cpuLock,(uint64_t)AccCache.accLock); 
  if (AccCache.accLock!=0) return;
  if (AccCache.cpuLock!=0) return;
  if(AccCache.state==AccDirty) {
    Flush(AccCache);
  }
  if(AccCache.AccPtr) {
    AcceleratorFree((void *)AccCache.AccPtr,AccCache.bytes);
    LRUremove(AccCache);
    AccCache.AccPtr=(uint64_t)NULL;
    AccCache.state=CpuDirty; // CPU primary now
    DeviceBytes   -=AccCache.bytes;
    dprintf("MemoryManager: Free(AccPtr %lx) footprint now %ld \n",(uint64_t)AccCache.AccPtr,DeviceBytes);  
  }
  //  uint64_t CpuPtr = AccCache.CpuPtr;
  DeviceEvictions++;
  //  EntryErase(CpuPtr);
}
void MemoryManager::Flush(AcceleratorViewEntry &AccCache)
{
  assert(AccCache.state==AccDirty);
  assert(AccCache.cpuLock==0);
  assert(AccCache.accLock==0);
  assert(AccCache.AccPtr!=(uint64_t)NULL);
  assert(AccCache.CpuPtr!=(uint64_t)NULL);
  acceleratorCopyFromDevice((void *)AccCache.AccPtr,(void *)AccCache.CpuPtr,AccCache.bytes);
  mprintf("MemoryManager: acceleratorCopyFromDevice Flush AccPtr %lx -> CpuPtr %lx\n",(uint64_t)AccCache.AccPtr,(uint64_t)AccCache.CpuPtr); fflush(stdout);
  DeviceToHostBytes+=AccCache.bytes;
  DeviceToHostXfer++;
  AccCache.state=Consistent;
}
void MemoryManager::Clone(AcceleratorViewEntry &AccCache)
{
  assert(AccCache.state==CpuDirty);
  assert(AccCache.cpuLock==0);
  assert(AccCache.accLock==0);
  assert(AccCache.CpuPtr!=(uint64_t)NULL);
  if(AccCache.AccPtr==(uint64_t)NULL){
    AccCache.AccPtr=(uint64_t)AcceleratorAllocate(AccCache.bytes);
    DeviceBytes+=AccCache.bytes;
  }
  mprintf("MemoryManager: acceleratorCopyToDevice   Clone AccPtr %lx <- CpuPtr %lx\n",(uint64_t)AccCache.AccPtr,(uint64_t)AccCache.CpuPtr); fflush(stdout);
  acceleratorCopyToDevice((void *)AccCache.CpuPtr,(void *)AccCache.AccPtr,AccCache.bytes);
  HostToDeviceBytes+=AccCache.bytes;
  HostToDeviceXfer++;
  AccCache.state=Consistent;
}

void MemoryManager::CpuDiscard(AcceleratorViewEntry &AccCache)
{
  assert(AccCache.state!=Empty);
  assert(AccCache.cpuLock==0);
  assert(AccCache.accLock==0);
  assert(AccCache.CpuPtr!=(uint64_t)NULL);
  if(AccCache.AccPtr==(uint64_t)NULL){
    AccCache.AccPtr=(uint64_t)AcceleratorAllocate(AccCache.bytes);
    DeviceBytes+=AccCache.bytes;
  }
  AccCache.state=AccDirty;
}

/////////////////////////////////////////////////////////////////////////////////
// View management
/////////////////////////////////////////////////////////////////////////////////
void MemoryManager::ViewClose(void* Ptr,ViewMode mode)
{
  if( (mode==AcceleratorRead)||(mode==AcceleratorWrite)||(mode==AcceleratorWriteDiscard) ){
    dprintf("AcceleratorViewClose %lx\n",(uint64_t)Ptr);
    AcceleratorViewClose((uint64_t)Ptr);
  } else if( (mode==CpuRead)||(mode==CpuWrite)){
    CpuViewClose((uint64_t)Ptr);
  } else { 
    assert(0);
  }
}
void *MemoryManager::ViewOpen(void* _CpuPtr,size_t bytes,ViewMode mode,ViewAdvise hint)
{
  uint64_t CpuPtr = (uint64_t)_CpuPtr;
  if( (mode==AcceleratorRead)||(mode==AcceleratorWrite)||(mode==AcceleratorWriteDiscard) ){
    dprintf("AcceleratorViewOpen %lx\n",(uint64_t)CpuPtr);
    return (void *) AcceleratorViewOpen(CpuPtr,bytes,mode,hint);
  } else if( (mode==CpuRead)||(mode==CpuWrite)){
    return (void *)CpuViewOpen(CpuPtr,bytes,mode,hint);
  } else { 
    assert(0);
    return NULL;
  }
}
void  MemoryManager::EvictVictims(uint64_t bytes)
{
  assert(bytes<DeviceMaxBytes);
  while(bytes+DeviceLRUBytes > DeviceMaxBytes){
    if ( DeviceLRUBytes > 0){
      assert(LRU.size()>0);
      uint64_t victim = LRU.back(); // From the LRU
      auto AccCacheIterator = EntryLookup(victim);
      auto & AccCache = AccCacheIterator->second;
      Evict(AccCache);
    } else {
      return;
    }
  }
}
uint64_t MemoryManager::AcceleratorViewOpen(uint64_t CpuPtr,size_t bytes,ViewMode mode,ViewAdvise hint)
{
  ////////////////////////////////////////////////////////////////////////////
  // Find if present, otherwise get or force an empty
  ////////////////////////////////////////////////////////////////////////////
  if ( EntryPresent(CpuPtr)==0 ){
    EntryCreate(CpuPtr,bytes,mode,hint);
  }

  auto AccCacheIterator = EntryLookup(CpuPtr);
  auto & AccCache = AccCacheIterator->second;
  if (!AccCache.AccPtr) {
    EvictVictims(bytes); 
  } 
  assert((mode==AcceleratorRead)||(mode==AcceleratorWrite)||(mode==AcceleratorWriteDiscard));

  assert(AccCache.cpuLock==0);  // Programming error

  if(AccCache.state!=Empty) {
    dprintf("ViewOpen found entry %lx %lx : %ld %ld accLock %ld\n",
		    (uint64_t)AccCache.CpuPtr,
		    (uint64_t)CpuPtr,
		    (uint64_t)AccCache.bytes,
	            (uint64_t)bytes,
		    (uint64_t)AccCache.accLock);
    assert(AccCache.CpuPtr == CpuPtr);
    assert(AccCache.bytes  ==bytes);
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
  if(AccCache.state==Empty) {
    assert(AccCache.LRU_valid==0);
    AccCache.CpuPtr = CpuPtr;
    AccCache.AccPtr = (uint64_t)NULL;
    AccCache.bytes  = bytes;
    AccCache.state  = CpuDirty;   // Cpu starts primary
    if(mode==AcceleratorWriteDiscard){
      CpuDiscard(AccCache);
      AccCache.state  = AccDirty;   // Empty + AcceleratorWrite=> AccDirty
    } else if(mode==AcceleratorWrite){
      Clone(AccCache);
      AccCache.state  = AccDirty;   // Empty + AcceleratorWrite=> AccDirty
    } else {
      Clone(AccCache);
      AccCache.state  = Consistent; // Empty + AccRead => Consistent
    }
    AccCache.accLock= 1;
    dprintf("Copied Empty entry into device accLock= %d\n",AccCache.accLock);
  } else if(AccCache.state==CpuDirty ){
    if(mode==AcceleratorWriteDiscard) {
      CpuDiscard(AccCache);
      AccCache.state  = AccDirty;   // CpuDirty + AcceleratorWrite=> AccDirty
    } else if(mode==AcceleratorWrite) {
      Clone(AccCache);
      AccCache.state  = AccDirty;   // CpuDirty + AcceleratorWrite=> AccDirty
    } else {
      Clone(AccCache);
      AccCache.state  = Consistent; // CpuDirty + AccRead => Consistent
    }
    AccCache.accLock++;
    dprintf("CpuDirty entry into device ++accLock= %d\n",AccCache.accLock);
  } else if(AccCache.state==Consistent) {
    if((mode==AcceleratorWrite)||(mode==AcceleratorWriteDiscard))
      AccCache.state  = AccDirty;   // Consistent + AcceleratorWrite=> AccDirty
    else
      AccCache.state  = Consistent; // Consistent + AccRead => Consistent
    AccCache.accLock++;
    dprintf("Consistent entry into device ++accLock= %d\n",AccCache.accLock);
  } else if(AccCache.state==AccDirty) {
    if((mode==AcceleratorWrite)||(mode==AcceleratorWriteDiscard))
      AccCache.state  = AccDirty; // AccDirty + AcceleratorWrite=> AccDirty
    else
      AccCache.state  = AccDirty; // AccDirty + AccRead => AccDirty
    AccCache.accLock++;
    dprintf("AccDirty entry ++accLock= %d\n",AccCache.accLock);
  } else {
    assert(0);
  }

  assert(AccCache.accLock>0);
  // If view is opened on device must remove from LRU
  if(AccCache.LRU_valid==1){
    // must possibly remove from LRU as now locked on GPU
    dprintf("AccCache entry removed from LRU \n");
    LRUremove(AccCache);
  }

  int transient =hint;
  AccCache.transient= transient? EvictNext : 0;

  return AccCache.AccPtr;
}
////////////////////////////////////
// look up & decrement lock count
////////////////////////////////////
void MemoryManager::AcceleratorViewClose(uint64_t CpuPtr)
{
  auto AccCacheIterator = EntryLookup(CpuPtr);
  auto & AccCache = AccCacheIterator->second;

  assert(AccCache.cpuLock==0);
  assert(AccCache.accLock>0);

  AccCache.accLock--;
  // Move to LRU queue if not locked and close on device
  if(AccCache.accLock==0) {
    dprintf("AccleratorViewClose %lx AccLock decremented to %ld move to LRU queue\n",(uint64_t)CpuPtr,(uint64_t)AccCache.accLock);
    LRUinsert(AccCache);
  } else {
    dprintf("AccleratorViewClose %lx AccLock decremented to %ld\n",(uint64_t)CpuPtr,(uint64_t)AccCache.accLock);
  }
}
void MemoryManager::CpuViewClose(uint64_t CpuPtr)
{
  auto AccCacheIterator = EntryLookup(CpuPtr);
  auto & AccCache = AccCacheIterator->second;

  assert(AccCache.cpuLock>0);
  assert(AccCache.accLock==0);

  AccCache.cpuLock--;
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
uint64_t MemoryManager::CpuViewOpen(uint64_t CpuPtr,size_t bytes,ViewMode mode,ViewAdvise transient)
{
  ////////////////////////////////////////////////////////////////////////////
  // Find if present, otherwise get or force an empty
  ////////////////////////////////////////////////////////////////////////////
  if ( EntryPresent(CpuPtr)==0 ){
    EntryCreate(CpuPtr,bytes,mode,transient);
  }

  auto AccCacheIterator = EntryLookup(CpuPtr);
  auto & AccCache = AccCacheIterator->second;

  // CPU doesn't need to free space
  //  if (!AccCache.AccPtr) {
  //    EvictVictims(bytes);
  //  }

  assert((mode==CpuRead)||(mode==CpuWrite));
  assert(AccCache.accLock==0);  // Programming error

  if(AccCache.state!=Empty) {
    assert(AccCache.CpuPtr == CpuPtr);
    assert(AccCache.bytes==bytes);
  }

  if(AccCache.state==Empty) {
    AccCache.CpuPtr = CpuPtr;
    AccCache.AccPtr = (uint64_t)NULL;
    AccCache.bytes  = bytes;
    AccCache.state  = CpuDirty; // Empty + CpuRead/CpuWrite => CpuDirty
    AccCache.accLock= 0;
    AccCache.cpuLock= 1;
  } else if(AccCache.state==CpuDirty ){
    // AccPtr dont care, deferred allocate
    AccCache.state = CpuDirty; // CpuDirty +CpuRead/CpuWrite => CpuDirty
    AccCache.cpuLock++;
  } else if(AccCache.state==Consistent) {
    assert(AccCache.AccPtr != (uint64_t)NULL);
    if(mode==CpuWrite)
      AccCache.state = CpuDirty;   // Consistent +CpuWrite => CpuDirty
    else 
      AccCache.state = Consistent; // Consistent +CpuRead  => Consistent
    AccCache.cpuLock++;
  } else if(AccCache.state==AccDirty) {
    assert(AccCache.AccPtr != (uint64_t)NULL);
    Flush(AccCache);
    if(mode==CpuWrite) AccCache.state = CpuDirty;   // AccDirty +CpuWrite => CpuDirty, Flush
    else            AccCache.state = Consistent; // AccDirty +CpuRead  => Consistent, Flush
    AccCache.cpuLock++;
  } else {
    assert(0); // should be unreachable
  }

  AccCache.transient= transient? EvictNext : 0;

  return AccCache.CpuPtr;
}
void  MemoryManager::NotifyDeletion(void *_ptr)
{
  // Look up in ViewCache
  uint64_t ptr = (uint64_t)_ptr;
  if(EntryPresent(ptr)) {
    auto e = EntryLookup(ptr);
    AccDiscard(e->second);
  }
}
void  MemoryManager::Print(void)
{
  PrintBytes();
  std::cout << GridLogMessage << "--------------------------------------------" << std::endl;
  std::cout << GridLogMessage << "Memory Manager                             " << std::endl;
  std::cout << GridLogMessage << "--------------------------------------------" << std::endl;
  std::cout << GridLogMessage << DeviceBytes   << " bytes allocated on device " << std::endl;
  std::cout << GridLogMessage << DeviceLRUBytes<< " bytes evictable on device " << std::endl;
  std::cout << GridLogMessage << DeviceMaxBytes<< " bytes max on device       " << std::endl;
  std::cout << GridLogMessage << HostToDeviceXfer << " transfers        to   device " << std::endl;
  std::cout << GridLogMessage << DeviceToHostXfer << " transfers        from device " << std::endl;
  std::cout << GridLogMessage << HostToDeviceBytes<< " bytes transfered to   device " << std::endl;
  std::cout << GridLogMessage << DeviceToHostBytes<< " bytes transfered from device " << std::endl;
  std::cout << GridLogMessage << DeviceEvictions  << " Evictions from device " << std::endl;
  std::cout << GridLogMessage << DeviceDestroy    << " Destroyed vectors on device " << std::endl;
  std::cout << GridLogMessage << AccViewTable.size()<< " vectors " << LRU.size()<<" evictable"<< std::endl;
  acceleratorMem();
  std::cout << GridLogMessage << "--------------------------------------------" << std::endl;
}
void  MemoryManager::PrintAll(void)
{
  Print();
  std::cout << GridLogMessage << std::endl;
  std::cout << GridLogMessage << "--------------------------------------------" << std::endl;
  std::cout << GridLogMessage << "CpuAddr\t\tAccAddr\t\tState\t\tcpuLock\taccLock\tLRU_valid "<<std::endl;
  std::cout << GridLogMessage << "--------------------------------------------" << std::endl;
  for(auto it=AccViewTable.begin();it!=AccViewTable.end();it++){
    auto &AccCache = it->second;
    
    std::string str;
    if ( AccCache.state==Empty    ) str = std::string("Empty");
    if ( AccCache.state==CpuDirty ) str = std::string("CpuDirty");
    if ( AccCache.state==AccDirty ) str = std::string("AccDirty");
    if ( AccCache.state==Consistent)str = std::string("Consistent");

    std::cout << GridLogMessage << "0x"<<std::hex<<AccCache.CpuPtr<<std::dec
	      << "\t0x"<<std::hex<<AccCache.AccPtr<<std::dec<<"\t" <<str
	      << "\t" << AccCache.cpuLock
	      << "\t" << AccCache.accLock
	      << "\t" << AccCache.LRU_valid<<std::endl;
  }
  std::cout << GridLogMessage << "--------------------------------------------" << std::endl;

};
int   MemoryManager::isOpen   (void* _CpuPtr) 
{ 
  uint64_t CpuPtr = (uint64_t)_CpuPtr;
  if ( EntryPresent(CpuPtr) ){
    auto AccCacheIterator = EntryLookup(CpuPtr);
    auto & AccCache = AccCacheIterator->second;
    return AccCache.cpuLock+AccCache.accLock;
  } else { 
    return 0;
  }
}
void MemoryManager::Audit(std::string s)
{
  uint64_t CpuBytes=0;
  uint64_t AccBytes=0;
  uint64_t LruBytes1=0;
  uint64_t LruBytes2=0;
  uint64_t LruCnt=0;
  
  std::cout << " Memory Manager::Audit() from "<<s<<std::endl;
  for(auto it=LRU.begin();it!=LRU.end();it++){
    uint64_t cpuPtr = *it;
    assert(EntryPresent(cpuPtr));
    auto AccCacheIterator = EntryLookup(cpuPtr);
    auto & AccCache = AccCacheIterator->second;
    LruBytes2+=AccCache.bytes;
    assert(AccCache.LRU_valid==1);
    assert(AccCache.LRU_entry==it);
  }
  std::cout << " Memory Manager::Audit() LRU queue matches table entries "<<std::endl;

  for(auto it=AccViewTable.begin();it!=AccViewTable.end();it++){
    auto &AccCache = it->second;
    
    std::string str;
    if ( AccCache.state==Empty    ) str = std::string("Empty");
    if ( AccCache.state==CpuDirty ) str = std::string("CpuDirty");
    if ( AccCache.state==AccDirty ) str = std::string("AccDirty");
    if ( AccCache.state==Consistent)str = std::string("Consistent");

    CpuBytes+=AccCache.bytes;
    if( AccCache.AccPtr )    AccBytes+=AccCache.bytes;
    if( AccCache.LRU_valid ) LruBytes1+=AccCache.bytes;
    if( AccCache.LRU_valid ) LruCnt++;
    
    if ( AccCache.cpuLock || AccCache.accLock ) {
      assert(AccCache.LRU_valid==0);

      std::cout << GridLogError << s<< "\n\t 0x"<<std::hex<<AccCache.CpuPtr<<std::dec
		<< "\t0x"<<std::hex<<AccCache.AccPtr<<std::dec<<"\t" <<str
		<< "\t cpuLock  " << AccCache.cpuLock
		<< "\t accLock  " << AccCache.accLock
		<< "\t LRUvalid " << AccCache.LRU_valid<<std::endl;
    }

    assert( AccCache.cpuLock== 0 ) ;
    assert( AccCache.accLock== 0 ) ;
  }
  std::cout << " Memory Manager::Audit() no locked table entries "<<std::endl;
  assert(LruBytes1==LruBytes2);
  assert(LruBytes1==DeviceLRUBytes);
  std::cout << " Memory Manager::Audit() evictable bytes matches sum over table "<<std::endl;
  assert(AccBytes==DeviceBytes);
  std::cout << " Memory Manager::Audit() device bytes matches sum over table "<<std::endl;
  assert(LruCnt == LRU.size());
  std::cout << " Memory Manager::Audit() LRU entry count matches "<<std::endl;

}

void MemoryManager::PrintState(void* _CpuPtr)
{
  uint64_t CpuPtr = (uint64_t)_CpuPtr;

  if ( EntryPresent(CpuPtr) ){
    auto AccCacheIterator = EntryLookup(CpuPtr);
    auto & AccCache = AccCacheIterator->second;
    std::string str;
    if ( AccCache.state==Empty    ) str = std::string("Empty");
    if ( AccCache.state==CpuDirty ) str = std::string("CpuDirty");
    if ( AccCache.state==AccDirty ) str = std::string("AccDirty");
    if ( AccCache.state==Consistent)str = std::string("Consistent");
    if ( AccCache.state==EvictNext) str = std::string("EvictNext");

    std::cout << GridLogMessage << "CpuAddr\t\tAccAddr\t\tState\t\tcpuLock\taccLock\tLRU_valid "<<std::endl;
    std::cout << GridLogMessage << "\tx"<<std::hex<<AccCache.CpuPtr<<std::dec
    << "\tx"<<std::hex<<AccCache.AccPtr<<std::dec<<"\t" <<str
    << "\t" << AccCache.cpuLock
    << "\t" << AccCache.accLock
    << "\t" << AccCache.LRU_valid<<std::endl;

  } else {
    std::cout << GridLogMessage << "No Entry in AccCache table." << std::endl; 
  }
}

NAMESPACE_END(Grid);

#endif
