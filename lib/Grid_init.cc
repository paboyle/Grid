/****************************************************************************/
/* PAB: Signal magic. Processor state dump is x86-64 specific               */
/****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h> 
#include <sys/time.h>
#include <signal.h>
#include <iostream>
#include <Grid.h>

#undef __X86_64
#define MAC

#ifdef MAC
#include <execinfo.h>
#endif

namespace Grid {

  std::streambuf *Grid_saved_stream_buf;
#if 0
  void Grid_quiesce_nodes(void)
  {
#ifdef GRID_COMMS_MPI
    int me;
    MPI_Comm_rank(MPI_COMM_WORLD,&me);
    std::streambuf* Grid_saved_stream_buf = std::cout.rdbuf();
    if ( me ) { 
      std::ofstream file("log.node");
      std::cout.rdbuf(file.rdbuf());
    }
#endif
  }
#endif
  void Grid_quiesce_nodes(void)
  {
#ifdef GRID_COMMS_MPI
    int me;
    MPI_Comm_rank(MPI_COMM_WORLD,&me);
    if ( me ) { 
      std::cout.setstate(std::ios::badbit);
    }
#endif
  }
  void Grid_unquiesce_nodes(void)
  {
#ifdef GRID_COMMS_MPI
    std::cout.clear();
#endif
  }

void Grid_init(int *argc,char ***argv)
{
#ifdef GRID_COMMS_MPI
  MPI_Init(argc,argv);
#endif
  Grid_debug_handler_init();
  Grid_quiesce_nodes();
}
void Grid_finalize(void)
{
#ifdef GRID_COMMS_MPI
  MPI_Finalize();
#endif
}
double usecond(void) {
  struct timeval tv;
  gettimeofday(&tv,NULL);
  return 1.0*tv.tv_usec + 1.0e6*tv.tv_sec;
}


#define _NBACKTRACE (256)
void * Grid_backtrace_buffer[_NBACKTRACE];

void Grid_sa_signal_handler(int sig,siginfo_t *si,void * ptr)
{
  printf("Caught signal %d\n",si->si_signo);
  printf("  mem address %lx\n",(uint64_t)si->si_addr);
  printf("         code %d\n",si->si_code);

#ifdef __X86_64
    ucontext_t * uc= (ucontext_t *)ptr;
  struct sigcontext *sc = (struct sigcontext *)&uc->uc_mcontext;
  printf("  instruction %llx\n",(uint64_t)sc->rip);
#define REG(A)  printf("  %s %lx\n",#A, sc-> A);
  REG(rdi);
  REG(rsi);
  REG(rbp);
  REG(rbx);
  REG(rdx);
  REG(rax);
  REG(rcx);
  REG(rsp);
  REG(rip);


  REG(r8);
  REG(r9);
  REG(r10);
  REG(r11);
  REG(r12);
  REG(r13);
  REG(r14);
  REG(r15);
#endif
#ifdef MAC
  int symbols    = backtrace        (Grid_backtrace_buffer,_NBACKTRACE);
  char **strings = backtrace_symbols(Grid_backtrace_buffer,symbols);
  for (int i = 0; i < symbols; i++){
    printf ("%s\n", strings[i]);
  }
#endif
  exit(0);
  return;
};

void Grid_debug_handler_init(void)
{
  struct sigaction sa,osa;
  sigemptyset (&sa.sa_mask);
  sa.sa_sigaction= Grid_sa_signal_handler;
  sa.sa_flags    = SA_SIGINFO;
  sigaction(SIGSEGV,&sa,NULL);
  sigaction(SIGTRAP,&sa,NULL);
}
}
