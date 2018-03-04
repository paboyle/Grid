#ifndef GRID_LEXICOGRAPHIC_H
#define GRID_LEXICOGRAPHIC_H


namespace Grid{

  class Lexicographic {
  public:

    template<class coor_t>
    static accelerator_inline void CoorFromIndex (coor_t& coor,int index,const coor_t &dims){
      int nd= dims.size();
      coor.resize(nd);
      for(int d=0;d<nd;d++){
	coor[d] = index % dims[d];
	index   = index / dims[d];
      }
    }

    template<class coor_t>
    static accelerator_inline void IndexFromCoor (const coor_t& coor,int &index,const coor_t &dims){
      int nd=dims.size();
      int stride=1;
      index=0;
      for(int d=0;d<nd;d++){
	index = index+stride*coor[d];
	stride=stride*dims[d];
      }
    }

    template<class coor_t>
    static accelerator_inline void IndexFromCoorReversed (const coor_t& coor,int &index,const coor_t &dims){
      int nd=dims.size();
      int stride=1;
      index=0;
      for(int d=nd-1;d>=0;d--){
	index = index+stride*coor[d];
	stride=stride*dims[d];
      }
    }
    template<class coor_t>
    static accelerator_inline void CoorFromIndexReversed (coor_t& coor,int index,const coor_t &dims){
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
