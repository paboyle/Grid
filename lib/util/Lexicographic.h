#ifndef GRID_LEXICOGRAPHIC_H
#define GRID_LEXICOGRAPHIC_H


namespace Grid{

  class Lexicographic {
  public:

    static inline void CoorFromIndex (std::vector<int>& coor,int index,std::vector<int> &dims){
      int nd= dims.size();
      coor.resize(nd);
      for(int d=0;d<nd;d++){
	coor[d] = index % dims[d];
	index   = index / dims[d];
      }
    }

    static inline void IndexFromCoor (std::vector<int>& coor,int &index,std::vector<int> &dims){
      int nd=dims.size();
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
