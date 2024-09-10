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

const int64_t Pages=32;
const int64_t PageWords=4096/sizeof(ComplexD);
const int64_t VecWords=PageWords*Pages;
const int64_t N=10000;

class Tester {
public:
  Vector<ComplexD>      zero_uvm;
  std::vector<ComplexD> zero_host;
  std::vector<Vector<ComplexD> > A;
  std::vector<std::vector<ComplexD> > B;
  uint64_t counter;

  Tester() :
    zero_uvm(VecWords,ComplexD(0.0)),
    zero_host(VecWords,ComplexD(0.0)),
    A(N,zero_uvm),
    B(N,zero_host)
  { counter = 0; }

  void  MemoryTest(int N)
  {
    for(int epoch = 0;epoch<100000;epoch++){
      
      int p  = random() %Pages; // Which address/page to hit
      int v  = random() %N; // Which vec
      int w  = random() %2; // Write or read
      int dev= random() %2; // On device?
      //    int e=1;
      ComplexD zc = counter++;
      
      if ( w ) {
	B[v][p*PageWords] = B[v][p*PageWords] + zc;
	if ( dev ) { 
	  ComplexD *A_v=&A[v][0];
	  accelerator_for(ss,1,1,{
	      A_v[p*PageWords] = A_v[p*PageWords] + zc;
	    });
	} else {
	  A[v][p*PageWords] = A[v][p*PageWords] + zc;
	}
      } else {
	if ( dev ) { 
	  ComplexD *A_v=&A[v][0];
	  ComplexD ref = B[v][p*PageWords];
	  std::cout << "Device compare "<<B[v][p*PageWords]<<std::endl;
	  accelerator_for(ss,1,1,{
	      assert(ref==A_v[p*PageWords]);
	    });
	} else {
	  std::cout << "Host compare "<<B[v][p*PageWords]<<std::endl;
	  assert(B[v][p*PageWords]==A[v][p*PageWords]);
	}
      }
    }
    
  }

};
int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  Tester test;

  for(int i=0;i<N;i++){
    std::cout << "============================"<<std::endl;
    std::cout << "Epoch "<<i<<"/"<<N<<std::endl;
    std::cout << "============================"<<std::endl;
    test.MemoryTest(32);
  }
  Grid_finalize();
}

