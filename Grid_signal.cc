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

#include "Grid.h"

#undef __X86_64
namespace dpo {

  void Grid_sa_signal_handler(int sig,siginfo_t *si,void * ptr);
  void Grid_debug_handler_init(void);

  void Grid_init(void)
  {
    Grid_debug_handler_init();
  }

void Grid_sa_signal_handler(int sig,siginfo_t *si,void * ptr)
{
         ucontext_t * uc= (ucontext_t *)ptr;

  printf("Caught signal %d\n",si->si_signo);
  printf("  mem address %llx\n",(uint64_t)si->si_addr);
  printf("         code %d\n",si->si_code);

#ifdef __X86_64
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

  fflush(stdout);

  if ( si->si_signo == SIGSEGV ) {
    printf("Grid_sa_signal_handler: Oops... this was a sigsegv you naughty naughty programmer. Goodbye\n");
    fflush(stdout);
    exit(-1);
  }
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
