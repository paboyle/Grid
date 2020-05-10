    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./benchmarks/Benchmark_dwf.cc

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
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
#include <Grid/Grid.h>

using namespace std;
using namespace Grid;

  Gamma::Algebra Gmu [] = {
    Gamma::Algebra::GammaX,
    Gamma::Algebra::GammaY,
    Gamma::Algebra::GammaZ,
    Gamma::Algebra::GammaT
  };

void benchDw(std::vector<int> & L, int Ls);

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);


  const int Ls=12;
  std::vector< std::vector<int> > latts;
#if 1
  latts.push_back(std::vector<int> ({24,24,24,24}) );
  latts.push_back(std::vector<int> ({48,24,24,24}) );
  latts.push_back(std::vector<int> ({96,24,24,24}) );
  latts.push_back(std::vector<int> ({96,48,24,24}) );
  //  latts.push_back(std::vector<int> ({96,48,48,24}) );
  //  latts.push_back(std::vector<int> ({96,48,48,48}) );
#else
  //  latts.push_back(std::vector<int> ({96,48,48,48}) );
  latts.push_back(std::vector<int> ({96,96,96,192}) );
#endif

  std::cout << GridLogMessage<< "*****************************************************************" <<std::endl;
  std::cout << GridLogMessage<< "* Kernel options --dslash-generic, --dslash-unroll, --dslash-asm" <<std::endl;
  std::cout << GridLogMessage<< "*****************************************************************" <<std::endl;

  if ( WilsonKernelsStatic::Opt == WilsonKernelsStatic::OptGeneric   ) std::cout << GridLogMessage<< "* Using GENERIC Nc WilsonKernels" <<std::endl;
  if ( WilsonKernelsStatic::Opt == WilsonKernelsStatic::OptHandUnroll) std::cout << GridLogMessage<< "* Using Nc=3       WilsonKernels" <<std::endl;
  if ( WilsonKernelsStatic::Opt == WilsonKernelsStatic::OptInlineAsm ) std::cout << GridLogMessage<< "* Using Asm Nc=3   WilsonKernels" <<std::endl;
  std::cout << GridLogMessage<< "*****************************************************************" <<std::endl;

  int threads = GridThread::GetThreads();
  std::cout<<GridLogMessage << "Grid is setup to use "<<threads<<" threads"<<std::endl;
  std::cout<<GridLogMessage << "=========================================================================="<<std::endl;
  std::cout<<GridLogMessage << "= Benchmarking DWF"<<std::endl;
  std::cout<<GridLogMessage << "=========================================================================="<<std::endl;
  std::cout<<GridLogMessage << "Volume \t\t\tProcs \t SchurDiagOne "<<std::endl;
  std::cout<<GridLogMessage << "=========================================================================="<<std::endl;

  for (int l=0;l<latts.size();l++){
    std::vector<int> latt4 = latts[l];
    std::cout << GridLogMessage <<"\t";
    for(int d=0;d<Nd;d++){
      std::cout<<latt4[d]<<"x";
    }
    std::cout <<Ls<<"\t" ;
    benchDw (latt4,Ls);
  }
  std::cout<<GridLogMessage << "=========================================================================="<<std::endl;
  Grid_finalize();
}


void benchDw(std::vector<int> & latt4, int Ls)
{
  /////////////////////////////////////////////////////////////////////////////////////
  // for Nc=3
  /////////////////////////////////////////////////////////////////////////////////////
  // Dw :  Ls*24*(7+48)= Ls*1320 
  //
  // M5D:  Ls*(4*2*Nc mul + 4*2*Nc madd ) = 3*4*2*Nc*Ls = Ls*72
  // Meo:  Ls*24*(7+48) + Ls*72 = Ls*1392 
  //
  // Mee:  3*Ns*2*Nc*Ls  // Chroma 6*N5*Nc*Ns 
  //
  // LeemInv : 2*2*Nc*madd*Ls
  // LeeInv  : 2*2*Nc*madd*Ls
  // DeeInv  : 4*2*Nc*mul *Ls
  // UeeInv  : 2*2*Nc*madd*Ls
  // UeemInv : 2*2*Nc*madd*Ls = Nc*Ls*(8+8+8+8+8) = 40*Nc*Ls// Chroma (10*N5 - 8)*Nc*Ns ~ (40 N5 - 32)Nc flops
  // QUDA counts as dense LsxLs real matrix x Ls x NcNsNreim => Nc*4*2 x Ls^2 FMA = 16Nc Ls^2 flops
  // Mpc => 1452*cbvol*2*Ls flops // 
  //     => (1344+Ls*48)*Ls*cbvol*2 flops QUDA = 1920 @Ls=12 and 2112 @Ls=16
  /////////////////////////////////////////////////////////////////////////////////////
  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(latt4, GridDefaultSimd(Nd,vComplexF::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);
  //  long unsigned int single_site_flops     = 8*Nc*(7+16*Nc)*Ls;
  long unsigned int single_site_mpc_flops = 8*Nc*(7+16*Nc)*2*Ls + 40*Nc*2*Ls + 4*Nc*2*Ls;
  long unsigned int single_site_quda_flops = 8*Nc*(7+16*Nc)*2*Ls + 16*Nc*Ls*Ls + 4*Nc*2*Ls;
  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});


  ColourMatrixF cm = ComplexF(1.0,0.0);

  int ncall=300;
  RealD mass=0.1;
  RealD M5  =1.8;
  RealD NP = UGrid->_Nprocessors;
  double volume=1;  for(int mu=0;mu<Nd;mu++) volume=volume*latt4[mu];

  LatticeGaugeFieldF Umu(UGrid); Umu=Zero();
  MobiusFermionF Dw(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,1.5,0.5);
  
  LatticeFermionF src_o (FrbGrid); src_o=1.0;
  LatticeFermionF r_o   (FrbGrid); r_o=Zero();

  int order =151;
  SchurDiagOneOperator<MobiusFermionF,LatticeFermionF>  Mpc(Dw);
  Chebyshev<LatticeFermionF>      Cheby(0.0,60.0,order);

  {
    Mpc.Mpc(src_o,r_o);
    Mpc.Mpc(src_o,r_o);
    Mpc.Mpc(src_o,r_o);

    double t0=usecond();
    for(int i=0;i<ncall;i++){
      Mpc.Mpc(src_o,r_o);
    }
    double t1=usecond();

    double flops=(single_site_mpc_flops*volume*ncall); // Mpc has 1 - Moo^-1 Moe Mee^-1 Meo  so CB cancels.
    std::cout <<"\t"<<NP<< "\t"<<flops/(t1-t0);
    flops=(single_site_quda_flops*volume*ncall);
    std::cout <<"\t"<<flops/(t1-t0)<<"\t"<<(t1-t0)/1000./1000.<<" s\t";

    // Cheby uses MpcDagMpc so 2x flops
    for(int i=0;i<1;i++){
    Cheby(Mpc,src_o,r_o);
    t0=usecond();
    Cheby(Mpc,src_o,r_o);
    t1=usecond();
    flops=(single_site_mpc_flops*volume*2*order);
    std::cout <<"\t"<<flops/(t1-t0);
    flops=(single_site_quda_flops*volume*2*order);
    std::cout <<"\t"<<flops/(t1-t0) << "\t" << (t1-t0)/1000./1000. <<" s";
    std::cout <<std::endl;
    }
  }
  //  Dw.Report();
}



