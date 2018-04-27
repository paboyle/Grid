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
      int stride=1;
      index=0;
      for(int d=0;d<nd;d++){
	index = index+stride*coor[d];
	stride=stride*dims[d];
      }
    }

    static inline void IndexFromCoorReversed (const std::vector<int>& coor,int &index,const std::vector<int> &dims){
      int nd=dims.size();
      int stride=1;
      index=0;
      for(int d=nd-1;d>=0;d--){
	index = index+stride*coor[d];
	stride=stride*dims[d];
      }
    }
    static inline void CoorFromIndexReversed (std::vector<int>& coor,int index,const std::vector<int> &dims){
      int nd= dims.size();
      coor.resize(nd);
      for(int d=nd-1;d>=0;d--){
	coor[d] = index % dims[d];
	index   = index / dims[d];
      }
    }


  };

}
#endif
