    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_memory_manager.cc

    Copyright (C) 2022

Author: Peter Boyle <pboyle@bnl.gov>

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

void  MemoryTest(GridCartesian         * FGrid,int N);

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());

  int N=100;
  for(int i=0;i<N;i++){
    std::cout << "============================"<<std::endl;
    std::cout << "Epoch "<<i<<"/"<<N<<std::endl;
    std::cout << "============================"<<std::endl;
    MemoryTest(UGrid,256);
    MemoryManager::Print();
    AUDIT();
  }
  Grid_finalize();
}

void  MemoryTest(GridCartesian         * FGrid, int N)
{
  LatticeComplexD zero(FGrid); zero=Zero();
  std::vector<LatticeComplexD> A(N,zero);//FGrid);

  std::vector<ComplexD> B(N,ComplexD(0.0)); // Update sequentially on host

  for(int v=0;v<N;v++) A[v] = Zero();

  uint64_t counter = 0;
  for(int epoch = 0;epoch<10000;epoch++){

    int v  = random() %N; // Which vec
    int w  = random() %2; // Write or read
    int e  = random() %3; // expression or for loop
    int dev= random() %2; // On device?
    //    int e=1;
    ComplexD zc = counter++;
    
    if ( w ) {
      B[v] = B[v] + zc;
      if ( e == 0 ) {
	A[v] = A[v] + zc - A[v] + A[v];
      } else {
	if ( dev ) { 
	  autoView(A_v,A[v],AcceleratorWrite);
	  accelerator_for(ss,FGrid->oSites(),1,{
	    A_v[ss] = A_v[ss] + zc;
	    });
	} else {
	  autoView(A_v,A[v],CpuWrite);
	  thread_for(ss,FGrid->oSites(),{
	      A_v[ss] = A_v[ss] + zc;
	    });
	}
      }
    } else {
      if ( e == 0 ) {
	A[v] = A[v] + A[v] - A[v];
      } else { 
	if ( dev ) { 
	  autoView(A_v,A[v],AcceleratorRead);
	  accelerator_for(ss,FGrid->oSites(),1,{
	      //	      assert(B[v]==A_v[ss]()()().getlane(0));
	    });
	  //	std::cout << "["<<v<<"] checked on GPU"<<B[v]<<std::endl;
	} else {
	  autoView(A_v,A[v],CpuRead);
	  thread_for(ss,FGrid->oSites(),{
	      assert(B[v]==A_v[ss]()()().getlane(0));
	    });
	  //	std::cout << "["<<v<<"] checked on CPU"<<B[v]<<std::endl;
	}
      }    
    }
  }

}
