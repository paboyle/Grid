#ifndef _GRID_STAT_H
#define _GRID_STAT_H

#ifdef AVX512
#define _KNIGHTS_LANDING_ROOTONLY
#endif

namespace Grid { 

///////////////////////////////////////////////////////////////////////////////
// Extra KNL counters from MCDRAM
///////////////////////////////////////////////////////////////////////////////
#ifdef _KNIGHTS_LANDING_
#define NMC 6
#define NEDC 8
struct ctrs
{
    uint64_t mcrd[NMC];
    uint64_t mcwr[NMC];
    uint64_t edcrd[NEDC]; 
    uint64_t edcwr[NEDC];
    uint64_t edchite[NEDC];
    uint64_t edchitm[NEDC];
    uint64_t edcmisse[NEDC];
    uint64_t edcmissm[NEDC];
};
// Peter/Azusa:
// Our modification of a code provided by Larry Meadows from Intel
// Verified by email exchange non-NDA, ok for github. Should be as uses /sys/devices/ FS
// so is already public and in the linux kernel for KNL.
struct knl_gbl_
{
  int mc_rd[NMC];
  int mc_wr[NMC];
  int edc_rd[NEDC];
  int edc_wr[NEDC];
  int edc_hite[NEDC];
  int edc_hitm[NEDC];
  int edc_misse[NEDC];
  int edc_missm[NEDC];
};
#endif
///////////////////////////////////////////////////////////////////////////////

class PmuStat
{
    uint64_t counters[8][256];
#ifdef _KNIGHTS_LANDING_
    static struct knl_gbl_ gbl;
#endif
    const char *name;

    uint64_t reads;     // memory reads
    uint64_t writes;    // memory writes
    uint64_t mrstart;   // memory read counter at start of parallel region
    uint64_t mrend;     // memory read counter at end of parallel region
    uint64_t mwstart;   // memory write counter at start of parallel region
    uint64_t mwend;     // memory write counter at end of parallel region

    // cumulative counters
    uint64_t count;     // number of invocations
    uint64_t tregion;   // total time in parallel region (from thread 0)
    uint64_t tcycles;   // total cycles inside parallel region
    uint64_t inst, ref, cyc;   // fixed counters
    uint64_t pmc0, pmc1;// pmu
    // add memory counters here
    // temp variables
    uint64_t tstart;    // tsc at start of parallel region
    uint64_t tend;      // tsc at end of parallel region
    // map for ctrs values
    // 0 pmc0 start
    // 1 pmc0 end
    // 2 pmc1 start
    // 3 pmc1 end
    // 4 tsc start
    // 5 tsc end
    static bool pmu_initialized;
public:
    static bool is_init(void){ return pmu_initialized;}
    static void pmu_init(void);
    static void pmu_fini(void);
    static void pmu_start(void);
    static void pmu_stop(void);
    void accum(int nthreads);
    static void xmemctrs(uint64_t *mr, uint64_t *mw);
    void start(void);
    void enter(int t);
    void exit(int t);
    void print(void);
    void init(const char *regname);
    void clear(void);
#ifdef _KNIGHTS_LANDING_
    static void     KNLsetup(void);
    static uint64_t KNLreadctr(int fd);
    static void     KNLreadctrs(ctrs &c);
    static void     KNLevsetup(const char *ename, int &fd, int event, int umask);
#endif
    
  };

}
#endif


