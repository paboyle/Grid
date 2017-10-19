#ifndef GRID_LEXICOGRAPHIC_H
#define GRID_LEXICOGRAPHIC_H


namespace Grid{

  class Lexicographic {
  public:

    static inline void CoorFromIndex (std::vector<int>& coor,int index,const std::vector<int> &dims){
      int nd= dims.size();
      coor.resize(nd);
      for(int d=0;d<nd;d++){
	coor[d] = index % dims[d];
	index   = index / dims[d];
      }
    }

    static inline void IndexFromCoor (const std::vector<int>& coor,int &index,const std::vector<int> &dims){
      int nd=dims.size();
      if(nd > coor.size())  {
	std::cout<< "coor.size "<<coor.size()<<" >dims.size "<<dims.size()<<std::endl; 
	assert(0);
	}
      int stride=1;
      index=0;
      for(int d=0;d<nd;d++){
	index = index+stride*coor[d];
	stride=stride*dims[d];
      }
    }

  };

}
#endif
