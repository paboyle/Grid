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
