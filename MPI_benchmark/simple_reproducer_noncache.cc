// Minimal reproducer for the CXI uncached-registration defect.
//
// Hypothesis: with the libfabric MR cache disabled, a device pointer at a non-zero offset
// into a hipMalloc'd allocation is registered as the whole allocation and the offset is
// discarded, so a send transmits from the base and a receive lands at the base.
//
// Each test therefore prints, for both the data and its landing place, what a correct
// implementation must give, what the hypothesis predicts instead, and what was observed.
//
//   hipcc -O2 -o simple_reproducer_noncache simple_reproducer_noncache.cc \
//         -I$MPICH_DIR/include -L$MPICH_DIR/lib -lmpi -L$MPICH_DIR/gtl/lib -lmpi_gtl_hsa
//
//   MPICH_GPU_SUPPORT_ENABLED=1 FI_MR_CACHE_MAX_COUNT=0 \
//     srun -N2 -n2 --ntasks-per-node=1 ./simple_reproducer_noncache

#include <mpi.h>
#include <hip/hip_runtime.h>
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <vector>

// The message must be large enough to go rendezvous: an eager message is staged through a host
// bounce buffer with hipMemcpy, which honours the offset, and never registers the user pointer.
static const int NWORDS = 262144; // 2 MB buffers
static int       NMSG   = 16384;  // 128 KB messages          (argv[1], in words)
static int       OFF    = 16384;  // offending offset in words (argv[2]); clearest at >= NMSG,
                                  // which keeps the intended and predicted windows disjoint

static int rank, size, peer;
static uint64_t *A, *B;
static std::vector<uint64_t> h(NWORDS);

#define HIP(cmd) do { hipError_t e=(cmd); if ( e != hipSuccess ) {			\
      printf("hip error %s at line %d\n",hipGetErrorString(e),__LINE__);		\
      MPI_Abort(MPI_COMM_WORLD,1); } } while(0)

static uint64_t word(uint64_t tag,int r,int n)
{
  return (tag<<60) | ((uint64_t)r<<32) | (uint64_t)n;
}

static const char *verdict(int observed,int correct,int predicted)
{
  if ( observed == correct   ) return "CORRECT";
  if ( observed == predicted ) return "WRONG, as predicted";
  return "WRONG, and not as predicted";
}

static void refresh(void)
{
  for(int n=0;n<NWORDS;n++) h[n] = word(0xa,rank,n);
  HIP(hipMemcpy(A,h.data(),NWORDS*sizeof(uint64_t),hipMemcpyHostToDevice));
  for(int n=0;n<NWORDS;n++) h[n] = word(0xb,rank,n);
  HIP(hipMemcpy(B,h.data(),NWORDS*sizeof(uint64_t),hipMemcpyHostToDevice));
  MPI_Barrier(MPI_COMM_WORLD);
}

static void exchange(int sendoff,int recvoff)
{
  MPI_Request req[2];
  MPI_Irecv(B+recvoff,NMSG*sizeof(uint64_t),MPI_BYTE,peer,0,MPI_COMM_WORLD,&req[0]);
  MPI_Isend(A+sendoff,NMSG*sizeof(uint64_t),MPI_BYTE,peer,0,MPI_COMM_WORLD,&req[1]);
  MPI_Waitall(2,req,MPI_STATUSES_IGNORE);
}

// Report the state of one NMSG-word window of B: how much of it was overwritten, and by what.
static void window(const char *label,int at)
{
  int mod=0, m0=-1;
  for(int i=0;i<NMSG && at+i<NWORDS;i++) {
    if ( h[at+i] != word(0xb,rank,at+i) ) { if ( m0 < 0 ) m0 = at+i; mod++; }
  }
  if ( mod == 0 ) {
    printf("rank %d    %-9s B[%d..%d] : untouched\n",rank,label,at,at+NMSG-1);
    return;
  }
  uint64_t w = h[m0];
  printf("rank %d    %-9s B[%d..%d] : %d of %d words overwritten, B[%d] holds A[%d] from rank %d\n",
	 rank,label,at,at+NMSG-1,mod,NMSG,m0,(int)(w & 0xffffffff),(int)((w>>32) & 0xfffffff));
}

// Name a word: still the receive pattern, or a word of the neighbour's send buffer.
static const char *describe(uint64_t w,int idx)
{
  static char s[64];
  if      ( w == word(0xb,rank,idx) ) snprintf(s,sizeof(s),"untouched");
  else if ( (w>>60) == 0xa )          snprintf(s,sizeof(s),"A[%d] of rank %d",
					       (int)(w & 0xffffffff),(int)((w>>32) & 0xfffffff));
  else                                snprintf(s,sizeof(s),"unrecognised");
  return s;
}

static void report(const char *name,int sendoff,int recvoff)
{
  HIP(hipMemcpy(h.data(),B,NWORDS*sizeof(uint64_t),hipMemcpyDeviceToHost));

  int first=-1, last=-1;
  for(int n=0;n<NWORDS;n++) {
    if ( h[n] != word(0xb,rank,n) ) { if ( first < 0 ) first = n; last = n; }
  }

  printf("rank %d  %s : send &A[%d] -> recv &B[%d]\n",rank,name,sendoff,recvoff);
  if ( first < 0 ) { printf("rank %d    nothing was written\n",rank); return; }

  uint64_t w   = h[first];
  int      src = (int)(w & 0xffffffff);
  int      sr  = (int)((w>>32) & 0xfffffff);

  // The hypothesis predicts both offsets are discarded: the bytes at A[0], landing at B[0].
  printf("rank %d    data      : want A[%d]  predict A[0]  got A[%d] from rank %d   %s\n",
	 rank,sendoff,src,sr,verdict(src,sendoff,0));
  printf("rank %d    place     : want B[%d]  predict B[0]  got B[%d]                %s\n",
	 rank,recvoff,first,verdict(first,recvoff,0));

  if ( src == sendoff && first == recvoff ) return;

  // Both candidate landing sites in our receive buffer, and both candidate source words in the
  // neighbour's send buffer, which the fill pattern fixes by construction.
  printf("rank %d    Boff %d\n",rank,recvoff);
  printf("rank %d      B[0]      = 0x%016llx  %s\n",
	 rank,(unsigned long long)h[0],describe(h[0],0));
  if ( recvoff != 0 )
    printf("rank %d      B[%d]\t= 0x%016llx  %s\n",
	   rank,recvoff,(unsigned long long)h[recvoff],describe(h[recvoff],recvoff));
  printf("rank %d    Neighbour rank %d send buffer holds\n",rank,peer);
  printf("rank %d      A[0]\t= 0x%016llx\n",
	 rank,(unsigned long long)word(0xa,peer,0));
  if ( sendoff != 0 )
    printf("rank %d      A[%d]\t= 0x%016llx\n",
	   rank,sendoff,(unsigned long long)word(0xa,peer,sendoff));

  window("intended",recvoff);
  if ( recvoff != 0 ) window("predicted",0);
  if ( first < recvoff || last >= recvoff+NMSG ) {
    if ( first != 0 || last != NMSG-1 ) {
      printf("rank %d    stray     : modified words span B[%d..%d], outside both windows\n",
	     rank,first,last);
    }
  }

  int bad=0;
  for(int i=0;i<NMSG && first+i<NWORDS;i++) {
    if ( h[first+i] != word(0xa,peer,src+i) ) bad++;
  }
  if ( bad ) printf("rank %d    body      : %d of %d words are not a contiguous run from A[%d]\n",
		    rank,bad,NMSG,src);
}

// Let each rank print in turn: flush, then a barrier so rank r is on the page before r+1 starts.
static void ordered(void (*print)(void))
{
  for(int r=0;r<size;r++) {
    if ( r == rank ) { print(); fflush(stdout); }
    MPI_Barrier(MPI_COMM_WORLD);
  }
}

static const char *test_name;
static int         test_sendoff, test_recvoff;
static void        print_report(void) { report(test_name,test_sendoff,test_recvoff); }
static void        print_setup (void)
{
  printf("rank %d  A %p  B %p  message %d words (%d bytes)  offset %d words (%d bytes)\n",
	 rank,(void *)A,(void *)B,NMSG,(int)(NMSG*sizeof(uint64_t)),OFF,(int)(OFF*sizeof(uint64_t)));
}

int main(int argc,char **argv)
{
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  if ( size != 2 ) {
    if ( rank == 0 ) printf("run with two ranks, one per node\n");
    MPI_Finalize();
    return 1;
  }
  peer = 1-rank;

  if ( argc > 1 ) NMSG = atoi(argv[1]);
  if ( argc > 2 ) OFF  = atoi(argv[2]);
  if ( NMSG < 1 || OFF < 0 || OFF+NMSG > NWORDS ) {
    if ( rank == 0 ) printf("message %d and offset %d words do not fit in %d\n",NMSG,OFF,NWORDS);
    MPI_Finalize();
    return 1;
  }

  int ndev;
  HIP(hipGetDeviceCount(&ndev));
  if ( ndev < 1 ) {
    printf("rank %d sees no GPU\n",rank);
    MPI_Abort(MPI_COMM_WORLD,1);
  }
  HIP(hipSetDevice(rank%ndev));
  HIP(hipMalloc(&A,NWORDS*sizeof(uint64_t)));
  HIP(hipMalloc(&B,NWORDS*sizeof(uint64_t)));

  ordered(print_setup);

  const int  offsets[4][2] = { {0,0}, {OFF,0}, {0,OFF}, {OFF,OFF} };
  const char *names[4]     = { "test 1","test 2","test 3","test 4" };

  for(int t=0;t<4;t++) {
    refresh();
    exchange(offsets[t][0],offsets[t][1]);
    test_name = names[t]; test_sendoff = offsets[t][0]; test_recvoff = offsets[t][1];
    ordered(print_report);
  }

  HIP(hipFree(A));
  HIP(hipFree(B));
  MPI_Finalize();
  return 0;
}
