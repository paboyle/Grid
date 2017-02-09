#include <Grid/Grid.h>
#include <Grid/PerfCount.h>
#include <Grid/Stat.h>


namespace Grid { 


bool PmuStat::pmu_initialized=false;


void PmuStat::init(const char *regname)
{
#ifdef __x86_64__
  name = regname;
  if (!pmu_initialized)
    {
      std::cout<<"initialising pmu"<<std::endl;
      pmu_initialized = true;
      pmu_init();
    }
  clear();
#endif
}
void PmuStat::clear(void)
{
#ifdef __x86_64__
  count = 0;
  tregion = 0;
  pmc0 = 0;
  pmc1 = 0;
  inst = 0;
  cyc = 0;
  ref = 0;
  tcycles = 0;
  reads = 0;
  writes = 0;
#endif
}
void PmuStat::print(void)
{
#ifdef __x86_64__
  std::cout <<"Reg "<<std::string(name)<<":\n";
  std::cout <<"  region "<<tregion<<std::endl;
  std::cout <<"  cycles "<<tcycles<<std::endl;
  std::cout <<"  inst   "<<inst   <<std::endl;
  std::cout <<"  cyc    "<<cyc    <<std::endl;
  std::cout <<"  ref    "<<ref    <<std::endl;
  std::cout <<"  pmc0   "<<pmc0   <<std::endl;
  std::cout <<"  pmc1   "<<pmc1   <<std::endl;
  std::cout <<"  count  "<<count  <<std::endl;
  std::cout <<"  reads  "<<reads  <<std::endl;
  std::cout <<"  writes "<<writes <<std::endl;
#endif
}
void PmuStat::start(void)
{
#ifdef __x86_64__
  pmu_start();
  ++count;
  xmemctrs(&mrstart, &mwstart);
  tstart = __rdtsc();
#endif
}
void PmuStat::enter(int t)
{
#ifdef __x86_64__
  counters[0][t] = __rdpmc(0);
  counters[1][t] = __rdpmc(1);
  counters[2][t] = __rdpmc((1<<30)|0);
  counters[3][t] = __rdpmc((1<<30)|1);
  counters[4][t] = __rdpmc((1<<30)|2);
  counters[5][t] = __rdtsc();
#endif
}
void PmuStat::exit(int t)
{
#ifdef __x86_64__
  counters[0][t] = __rdpmc(0) - counters[0][t];
  counters[1][t] = __rdpmc(1) - counters[1][t];
  counters[2][t] = __rdpmc((1<<30)|0) - counters[2][t];
  counters[3][t] = __rdpmc((1<<30)|1) - counters[3][t];
  counters[4][t] = __rdpmc((1<<30)|2) - counters[4][t];
  counters[5][t] = __rdtsc() - counters[5][t];
#endif
}
void PmuStat::accum(int nthreads)
{
#ifdef __x86_64__
  tend = __rdtsc();
  xmemctrs(&mrend, &mwend);
  pmu_stop();
  for (int t = 0; t < nthreads; ++t) {
    pmc0 += counters[0][t];
    pmc1 += counters[1][t];
    inst += counters[2][t];
    cyc += counters[3][t];
    ref += counters[4][t];
    tcycles += counters[5][t];
  }
  uint64_t region = tend - tstart;
  tregion += region;
  uint64_t mreads = mrend - mrstart;
  reads += mreads;
  uint64_t mwrites = mwend - mwstart;
  writes += mwrites;
#endif
}


void PmuStat::pmu_fini(void) {}
void PmuStat::pmu_start(void) {};
void PmuStat::pmu_stop(void) {};
void PmuStat::pmu_init(void)
{
#ifdef _KNIGHTS_LANDING_
  KNLsetup();
#endif
}
void PmuStat::xmemctrs(uint64_t *mr, uint64_t *mw)
{
#ifdef _KNIGHTS_LANDING_
  ctrs c;
  KNLreadctrs(c);
  uint64_t emr = 0, emw = 0;
  for (int i = 0; i < NEDC; ++i)
    {
      emr += c.edcrd[i];
      emw += c.edcwr[i];
    }
  *mr = emr;
  *mw = emw;
#else
  *mr = *mw = 0;
#endif
}

#ifdef _KNIGHTS_LANDING_

struct knl_gbl_ PmuStat::gbl;

#define PMU_MEM

void PmuStat::KNLevsetup(const char *ename, int &fd, int event, int umask)
{
  char fname[1024];
  snprintf(fname, sizeof(fname), "%s/type", ename);
  FILE *fp = fopen(fname, "r");
  if (fp == 0) {
    ::printf("open %s", fname);
    ::exit(0);
  }
  int type;
  int ret = fscanf(fp, "%d", &type);
  assert(ret == 1);
  fclose(fp);
  //  std::cout << "Using PMU type "<<type<<" from " << std::string(ename) <<std::endl;

  struct perf_event_attr hw = {};
  hw.size = sizeof(hw);
  hw.type = type;
  // see /sys/devices/uncore_*/format/*
  // All of the events we are interested in are configured the same way, but
  // that isn't always true. Proper code would parse the format files
  hw.config = event | (umask << 8);
  //hw.read_format = PERF_FORMAT_GROUP;
  // unfortunately the above only works within a single PMU; might
  // as well just read them one at a time
  int cpu = 0;
  fd = perf_event_open(&hw, -1, cpu, -1, 0);
  if (fd == -1) {
    ::printf("CPU %d, box %s, event 0x%lx", cpu, ename, hw.config);
    ::exit(0);
  } else { 
    //    std::cout << "event "<<std::string(ename)<<" set up for fd "<<fd<<" hw.config "<<hw.config <<std::endl;
  }
}


 void PmuStat::KNLsetup(void){

   int ret;
   char fname[1024];

   // MC RPQ inserts and WPQ inserts (reads & writes)
   for (int mc = 0; mc < NMC; ++mc)
     {
       ::snprintf(fname, sizeof(fname), "/sys/devices/uncore_imc_%d",mc);
       // RPQ Inserts
       KNLevsetup(fname, gbl.mc_rd[mc], 0x1, 0x1);
       // WPQ Inserts
       KNLevsetup(fname, gbl.mc_wr[mc], 0x2, 0x1);
     }
   // EDC RPQ inserts and WPQ inserts
   for (int edc=0; edc < NEDC; ++edc)
     {
       ::snprintf(fname, sizeof(fname), "/sys/devices/uncore_edc_eclk_%d",edc);
       // RPQ inserts
       KNLevsetup(fname, gbl.edc_rd[edc], 0x1, 0x1);
       // WPQ inserts
       KNLevsetup(fname, gbl.edc_wr[edc], 0x2, 0x1);
     }
   // EDC HitE, HitM, MissE, MissM
   for (int edc=0; edc < NEDC; ++edc)
     {
       ::snprintf(fname, sizeof(fname), "/sys/devices/uncore_edc_uclk_%d", edc);
       KNLevsetup(fname, gbl.edc_hite[edc], 0x2, 0x1);
       KNLevsetup(fname, gbl.edc_hitm[edc], 0x2, 0x2);
       KNLevsetup(fname, gbl.edc_misse[edc], 0x2, 0x4);
       KNLevsetup(fname, gbl.edc_missm[edc], 0x2, 0x8);
     }
 }

uint64_t PmuStat::KNLreadctr(int fd)
{
  uint64_t data;
  size_t s = ::read(fd, &data, sizeof(data));
  if (s != sizeof(uint64_t)){
    ::printf("read counter %lu", s);
    ::exit(0);
  }
  return data;
}

void PmuStat::KNLreadctrs(ctrs &c)
{
  for (int i = 0; i < NMC; ++i)
    {
      c.mcrd[i] = KNLreadctr(gbl.mc_rd[i]);
      c.mcwr[i] = KNLreadctr(gbl.mc_wr[i]);
    }
  for (int i = 0; i < NEDC; ++i)
    {
      c.edcrd[i] = KNLreadctr(gbl.edc_rd[i]);
      c.edcwr[i] = KNLreadctr(gbl.edc_wr[i]);
    }
  for (int i = 0; i < NEDC; ++i)
    {
      c.edchite[i] = KNLreadctr(gbl.edc_hite[i]);
      c.edchitm[i] = KNLreadctr(gbl.edc_hitm[i]);
      c.edcmisse[i] = KNLreadctr(gbl.edc_misse[i]);
      c.edcmissm[i] = KNLreadctr(gbl.edc_missm[i]);
    }
}

#endif
}
