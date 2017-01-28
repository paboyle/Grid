    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/Init.cc

    Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: Peter Boyle <peterboyle@MacBook-Pro.local>
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
/****************************************************************************/
/* pab: Signal magic. Processor state dump is x86-64 specific               */
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
#include <iterator>
#include <Grid/Grid.h>
#include <algorithm>
#include <iterator>
#include <cstdlib>
#include <memory>


#include <fenv.h>
#ifdef __APPLE__
static int
feenableexcept (unsigned int excepts)
{
  static fenv_t fenv;
  unsigned int new_excepts = excepts & FE_ALL_EXCEPT,
    old_excepts;  // previous masks

  if ( fegetenv (&fenv) ) return -1;
  old_excepts = fenv.__control & FE_ALL_EXCEPT;

  // unmask
  fenv.__control &= ~new_excepts;
  fenv.__mxcsr   &= ~(new_excepts << 7);

  return ( fesetenv (&fenv) ? -1 : old_excepts );
}
#endif

namespace Grid {


//////////////////////////////////////////////////////
// Convenience functions to access stadard command line arg
// driven parallelism controls
//////////////////////////////////////////////////////
static std::vector<int> Grid_default_latt;
static std::vector<int> Grid_default_mpi;

int GridThread::_threads =1;
int GridThread::_hyperthreads=1;
int GridThread::_cores=1;

const std::vector<int> &GridDefaultLatt(void)     {return Grid_default_latt;};
const std::vector<int> &GridDefaultMpi(void)      {return Grid_default_mpi;};
const std::vector<int> GridDefaultSimd(int dims,int nsimd)
{
    std::vector<int> layout(dims);
    int nn=nsimd;
    for(int d=dims-1;d>=0;d--){
      if ( nn>=2) {
	layout[d]=2;
	nn/=2;
      } else { 
	layout[d]=1;
      }
    }
    assert(nn==1);
    return layout;
}
  
////////////////////////////////////////////////////////////
// Command line parsing assist for stock controls
////////////////////////////////////////////////////////////
std::string GridCmdOptionPayload(char ** begin, char ** end, const std::string & option)
{
  char ** itr = std::find(begin, end, option);
  if (itr != end && ++itr != end) {
    std::string payload(*itr);
    return payload;
  }
  return std::string("");
}
bool GridCmdOptionExists(char** begin, char** end, const std::string& option)
{
  return std::find(begin, end, option) != end;
}
  // Comma separated list
void GridCmdOptionCSL(std::string str,std::vector<std::string> & vec)
{
  size_t pos = 0;
  std::string token;
  std::string delimiter(",");

  vec.resize(0);
  while ((pos = str.find(delimiter)) != std::string::npos) {
    token = str.substr(0, pos);
    vec.push_back(token);
    str.erase(0, pos + delimiter.length());
  }
  token = str;
  vec.push_back(token);
  return;
}

void GridCmdOptionIntVector(std::string &str,std::vector<int> & vec)
{
  vec.resize(0);
  std::stringstream ss(str);
  int i;
  while (ss >> i){
    vec.push_back(i);
    if(std::ispunct(ss.peek()))
      ss.ignore();
  }    
  return;
}

void GridCmdOptionInt(std::string &str,int & val)
{
  std::stringstream ss(str);
  ss>>val;
  return;
}


void GridParseLayout(char **argv,int argc,
		     std::vector<int> &latt,
		     std::vector<int> &mpi)
{
  mpi =std::vector<int>({1,1,1,1});
  latt=std::vector<int>({8,8,8,8});

  GridThread::SetMaxThreads();

  std::string arg;
  if( GridCmdOptionExists(argv,argv+argc,"--mpi") ){
    arg = GridCmdOptionPayload(argv,argv+argc,"--mpi");
    GridCmdOptionIntVector(arg,mpi);
  }
  if( GridCmdOptionExists(argv,argv+argc,"--grid") ){
    arg= GridCmdOptionPayload(argv,argv+argc,"--grid");
    GridCmdOptionIntVector(arg,latt);
  }
  if( GridCmdOptionExists(argv,argv+argc,"--threads") ){
    std::vector<int> ompthreads(0);
#ifndef GRID_OMP
    std::cout << GridLogWarning << "'--threads' option used but Grid was"
              << " not compiled with thread support" << std::endl;
#endif
    arg= GridCmdOptionPayload(argv,argv+argc,"--threads");
    GridCmdOptionIntVector(arg,ompthreads);
    assert(ompthreads.size()==1);
    GridThread::SetThreads(ompthreads[0]);
  }
  if( GridCmdOptionExists(argv,argv+argc,"--cores") ){
    int cores;
    arg= GridCmdOptionPayload(argv,argv+argc,"--cores");
    GridCmdOptionInt(arg,cores);
    GridThread::SetCores(cores);
  }
}

std::string GridCmdVectorIntToString(const std::vector<int> & vec){
  std::ostringstream oss;
  std::copy(vec.begin(), vec.end(),std::ostream_iterator<int>(oss, " "));
  return oss.str();
}
/////////////////////////////////////////////////////////
// Reinit guard
/////////////////////////////////////////////////////////
static int Grid_is_initialised = 0;


void Grid_init(int *argc,char ***argv)
{
  GridLogger::StopWatch.Start();

  std::string arg;

  ////////////////////////////////////
  // Shared memory block size
  ////////////////////////////////////
  if( GridCmdOptionExists(*argv,*argv+*argc,"--shm") ){
    int MB;
    arg= GridCmdOptionPayload(*argv,*argv+*argc,"--shm");
    GridCmdOptionInt(arg,MB);
    CartesianCommunicator::MAX_MPI_SHM_BYTES = MB*1024*1024;
  }

  CartesianCommunicator::Init(argc,argv);

  ////////////////////////////////////
  // Logging
  ////////////////////////////////////

  std::vector<std::string> logstreams;
  std::string defaultLog("Error,Warning,Message,Performance");
  GridCmdOptionCSL(defaultLog,logstreams);
  GridLogConfigure(logstreams);

  if( !GridCmdOptionExists(*argv,*argv+*argc,"--debug-stdout") ){
    Grid_quiesce_nodes();
  }

  if( GridCmdOptionExists(*argv,*argv+*argc,"--log") ){
    arg = GridCmdOptionPayload(*argv,*argv+*argc,"--log");
    GridCmdOptionCSL(arg,logstreams);
    GridLogConfigure(logstreams);
  }

  ////////////////////////////////////
  // Help message
  ////////////////////////////////////

  if( GridCmdOptionExists(*argv,*argv+*argc,"--help") ){
    std::cout<<GridLogMessage<<"  --help : this message"<<std::endl;
    std::cout<<GridLogMessage<<std::endl;
    std::cout<<GridLogMessage<<"Geometry:"<<std::endl;
    std::cout<<GridLogMessage<<"  --mpi n.n.n.n   : default MPI decomposition"<<std::endl;    
    std::cout<<GridLogMessage<<"  --threads n     : default number of OMP threads"<<std::endl;
    std::cout<<GridLogMessage<<"  --grid n.n.n.n  : default Grid size"<<std::endl;    
    std::cout<<GridLogMessage<<"  --shm  M        : allocate M megabytes of shared memory for comms"<<std::endl;    
    std::cout<<GridLogMessage<<std::endl;
    std::cout<<GridLogMessage<<"Verbose and debug:"<<std::endl;
    std::cout<<GridLogMessage<<"  --log list      : comma separted list of streams from Error,Warning,Message,Performance,Iterative,Integrator,Debug,Colours"<<std::endl;
    std::cout<<GridLogMessage<<"  --decomposition : report on default omp,mpi and simd decomposition"<<std::endl;    
    std::cout<<GridLogMessage<<"  --debug-signals : catch sigsegv and print a blame report"<<std::endl;
    std::cout<<GridLogMessage<<"  --debug-stdout  : print stdout from EVERY node"<<std::endl;    
    std::cout<<GridLogMessage<<"  --notimestamp   : suppress millisecond resolution stamps"<<std::endl;    
    std::cout<<GridLogMessage<<std::endl;
    std::cout<<GridLogMessage<<"Performance:"<<std::endl;
    std::cout<<GridLogMessage<<"  --dslash-generic: Wilson kernel for generic Nc"<<std::endl;    
    std::cout<<GridLogMessage<<"  --dslash-unroll : Wilson kernel for Nc=3"<<std::endl;    
    std::cout<<GridLogMessage<<"  --dslash-asm    : Wilson kernel for AVX512"<<std::endl;    
    std::cout<<GridLogMessage<<"  --lebesgue      : Cache oblivious Lebesgue curve/Morton order/Z-graph stencil looping"<<std::endl;    
    std::cout<<GridLogMessage<<"  --cacheblocking n.m.o.p : Hypercuboidal cache blocking"<<std::endl;    
    std::cout<<GridLogMessage<<std::endl;
    exit(EXIT_SUCCESS);
  }

  ////////////////////////////////////
  // Banner
  ////////////////////////////////////

  std::string COL_RED    = GridLogColours.colour["RED"];
  std::string COL_PURPLE = GridLogColours.colour["PURPLE"];
  std::string COL_BLACK  = GridLogColours.colour["BLACK"];
  std::string COL_GREEN  = GridLogColours.colour["GREEN"];
  std::string COL_BLUE   = GridLogColours.colour["BLUE"];
  std::string COL_YELLOW = GridLogColours.colour["YELLOW"];
  std::string COL_BACKGROUND = GridLogColours.colour["NORMAL"];
  
  std::cout <<std::endl;
  std::cout <<COL_RED  << "__|__|__|__|__"<<             "|__|__|_"<<COL_PURPLE<<"_|__|__|"<<                "__|__|__|__|__"<<std::endl; 
  std::cout <<COL_RED  << "__|__|__|__|__"<<             "|__|__|_"<<COL_PURPLE<<"_|__|__|"<<                "__|__|__|__|__"<<std::endl; 
  std::cout <<COL_RED  << "__|_ |  |  |  "<<             "|  |  | "<<COL_PURPLE<<" |  |  |"<<                "  |  |  | _|__"<<std::endl; 
  std::cout <<COL_RED  << "__|_          "<<             "        "<<COL_PURPLE<<"        "<<                "          _|__"<<std::endl; 
  std::cout <<COL_RED  << "__|_  "<<COL_GREEN<<" GGGG   "<<COL_RED<<" RRRR   "<<COL_BLUE  <<" III    "<<COL_PURPLE<<"DDDD  "<<COL_PURPLE<<"    _|__"<<std::endl;
  std::cout <<COL_RED  << "__|_  "<<COL_GREEN<<"G       "<<COL_RED<<" R   R  "<<COL_BLUE  <<"  I     "<<COL_PURPLE<<"D   D "<<COL_PURPLE<<"    _|__"<<std::endl;
  std::cout <<COL_RED  << "__|_  "<<COL_GREEN<<"G       "<<COL_RED<<" R   R  "<<COL_BLUE  <<"  I     "<<COL_PURPLE<<"D    D"<<COL_PURPLE<<"    _|__"<<std::endl;
  std::cout <<COL_BLUE << "__|_  "<<COL_GREEN<<"G  GG   "<<COL_RED<<" RRRR   "<<COL_BLUE  <<"  I     "<<COL_PURPLE<<"D    D"<<COL_GREEN <<"    _|__"<<std::endl;
  std::cout <<COL_BLUE << "__|_  "<<COL_GREEN<<"G   G   "<<COL_RED<<" R  R   "<<COL_BLUE  <<"  I     "<<COL_PURPLE<<"D   D "<<COL_GREEN <<"    _|__"<<std::endl;
  std::cout <<COL_BLUE << "__|_  "<<COL_GREEN<<" GGGG   "<<COL_RED<<" R   R  "<<COL_BLUE  <<" III    "<<COL_PURPLE<<"DDDD  "<<COL_GREEN <<"    _|__"<<std::endl;
  std::cout <<COL_BLUE << "__|_          "<<             "        "<<COL_GREEN <<"        "<<                "          _|__"<<std::endl; 
  std::cout <<COL_BLUE << "__|__|__|__|__"<<             "|__|__|_"<<COL_GREEN <<"_|__|__|"<<                "__|__|__|__|__"<<std::endl; 
  std::cout <<COL_BLUE << "__|__|__|__|__"<<             "|__|__|_"<<COL_GREEN <<"_|__|__|"<<                "__|__|__|__|__"<<std::endl; 
  std::cout <<COL_BLUE << "  |  |  |  |  "<<             "|  |  | "<<COL_GREEN <<" |  |  |"<<                "  |  |  |  |  "<<std::endl; 
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout <<COL_YELLOW<< std::endl;
  std::cout << "Copyright (C) 2015 Peter Boyle, Azusa Yamaguchi, Guido Cossu, Antonin Portelli and other authors"<<std::endl;
  std::cout << std::endl;
  std::cout << "This program is free software; you can redistribute it and/or modify"<<std::endl;
  std::cout << "it under the terms of the GNU General Public License as published by"<<std::endl;
  std::cout << "the Free Software Foundation; either version 2 of the License, or"<<std::endl;
  std::cout << "(at your option) any later version."<<std::endl;
  std::cout << std::endl;
  std::cout << "This program is distributed in the hope that it will be useful,"<<std::endl;
  std::cout << "but WITHOUT ANY WARRANTY; without even the implied warranty of"<<std::endl;
  std::cout << "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the"<<std::endl;
  std::cout << "GNU General Public License for more details."<<std::endl;
  std::cout << COL_BACKGROUND <<std::endl;
  std::cout << std::endl;

  ////////////////////////////////////
  // Debug and performance options
  ////////////////////////////////////


  if( GridCmdOptionExists(*argv,*argv+*argc,"--debug-signals") ){
    Grid_debug_handler_init();
  }
  if( GridCmdOptionExists(*argv,*argv+*argc,"--dslash-unroll") ){
    QCD::WilsonKernelsStatic::Opt=QCD::WilsonKernelsStatic::OptHandUnroll;
  }
  if( GridCmdOptionExists(*argv,*argv+*argc,"--dslash-asm") ){
    QCD::WilsonKernelsStatic::Opt=QCD::WilsonKernelsStatic::OptInlineAsm;
  }
  if( GridCmdOptionExists(*argv,*argv+*argc,"--dslash-generic") ){
    QCD::WilsonKernelsStatic::Opt=QCD::WilsonKernelsStatic::OptGeneric;
  }
  if( GridCmdOptionExists(*argv,*argv+*argc,"--lebesgue") ){
    LebesgueOrder::UseLebesgueOrder=1;
  }


  if( GridCmdOptionExists(*argv,*argv+*argc,"--cacheblocking") ){
    arg= GridCmdOptionPayload(*argv,*argv+*argc,"--cacheblocking");
    GridCmdOptionIntVector(arg,LebesgueOrder::Block);
  }
  if( GridCmdOptionExists(*argv,*argv+*argc,"--notimestamp") ){
    GridLogTimestamp(0);
  } else { 
    GridLogTimestamp(1);
  }

  GridParseLayout(*argv,*argc,
		  Grid_default_latt,
		  Grid_default_mpi);

  std::cout << GridLogMessage << "Requesting "<< CartesianCommunicator::MAX_MPI_SHM_BYTES <<" byte stencil comms buffers "<<std::endl;

  if( GridCmdOptionExists(*argv,*argv+*argc,"--decomposition") ){
    std::cout<<GridLogMessage<<"Grid Decomposition\n";
    std::cout<<GridLogMessage<<"\tOpenMP threads : "<<GridThread::GetThreads()<<std::endl;
    std::cout<<GridLogMessage<<"\tMPI tasks      : "<<GridCmdVectorIntToString(GridDefaultMpi())<<std::endl;
    std::cout<<GridLogMessage<<"\tvRealF         : "<<sizeof(vRealF)*8    <<"bits ; " <<GridCmdVectorIntToString(GridDefaultSimd(4,vRealF::Nsimd()))<<std::endl;
    std::cout<<GridLogMessage<<"\tvRealD         : "<<sizeof(vRealD)*8    <<"bits ; " <<GridCmdVectorIntToString(GridDefaultSimd(4,vRealD::Nsimd()))<<std::endl;
    std::cout<<GridLogMessage<<"\tvComplexF      : "<<sizeof(vComplexF)*8 <<"bits ; " <<GridCmdVectorIntToString(GridDefaultSimd(4,vComplexF::Nsimd()))<<std::endl;
    std::cout<<GridLogMessage<<"\tvComplexD      : "<<sizeof(vComplexD)*8 <<"bits ; " <<GridCmdVectorIntToString(GridDefaultSimd(4,vComplexD::Nsimd()))<<std::endl;
  }


  Grid_is_initialised = 1;
}

  
void Grid_finalize(void)
{
#if defined (GRID_COMMS_MPI) || defined (GRID_COMMS_MPI3) 
  MPI_Finalize();
  Grid_unquiesce_nodes();
#endif
}

void * Grid_backtrace_buffer[_NBACKTRACE];

void Grid_sa_signal_handler(int sig,siginfo_t *si,void * ptr)
{
  printf("Caught signal %d\n",si->si_signo);
  printf("  mem address %llx\n",(unsigned long long)si->si_addr);
  printf("         code %d\n",si->si_code);

  // Linux/Posix
#ifdef __linux__
  // And x86 64bit
#ifdef __x86_64__
  ucontext_t * uc= (ucontext_t *)ptr;
  struct sigcontext *sc = (struct sigcontext *)&uc->uc_mcontext;
  printf("  instruction %llx\n",(unsigned long long)sc->rip);
#define REG(A)  printf("  %s %lx\n",#A,sc-> A);
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
#endif
  BACKTRACE();
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

  feenableexcept( FE_INVALID|FE_OVERFLOW|FE_DIVBYZERO);

  sigaction(SIGFPE,&sa,NULL);
}
}
