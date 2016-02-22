    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/algorithms/iterative/EigenSort.h

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>

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
#ifndef GRID_EIGENSORT_H
#define GRID_EIGENSORT_H


namespace Grid {
    /////////////////////////////////////////////////////////////
    // Eigen sorter to begin with
    /////////////////////////////////////////////////////////////

template<class Field>
class SortEigen {
 private:
  
  static bool less_lmd(RealD left,RealD right){
    return fabs(left) < fabs(right);
  }  
  static bool less_pair(std::pair<RealD,Field>& left,
		 std::pair<RealD,Field>& right){
    return fabs(left.first) < fabs(right.first);
  }  
  
 public:

  void push(DenseVector<RealD>& lmd,
	    DenseVector<Field>& evec,int N) {

    DenseVector<std::pair<RealD, Field> > emod;
    typename DenseVector<std::pair<RealD, Field> >::iterator it;
    
    for(int i=0;i<lmd.size();++i){
      emod.push_back(std::pair<RealD,Field>(lmd[i],evec[i]));
    }

    partial_sort(emod.begin(),emod.begin()+N,emod.end(),less_pair);

    it=emod.begin();
    for(int i=0;i<N;++i){
      lmd[i]=it->first;
      evec[i]=it->second;
      ++it;
    }
  }
  void push(DenseVector<RealD>& lmd,int N) {
    std::partial_sort(lmd.begin(),lmd.begin()+N,lmd.end(),less_lmd);
  }
  bool saturated(RealD lmd, RealD thrs) {
    return fabs(lmd) > fabs(thrs);
  }
};

}
#endif
