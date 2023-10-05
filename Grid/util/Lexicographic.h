#ifndef GRID_LEXICOGRAPHIC_H
#define GRID_LEXICOGRAPHIC_H


namespace Grid{

  class Lexicographic {
  public:

    template<class coor_t>
    static accelerator_inline void CoorFromIndex (coor_t& coor,int64_t index,const coor_t &dims){
      int nd= dims.size();
      coor.resize(nd);
      for(int d=0;d<nd;d++){
	coor[d] = index % dims[d];
	index   = index / dims[d];
      }
    }

    template<class coor_t>
    static accelerator_inline void IndexFromCoor (const coor_t& coor,int64_t &index,const coor_t &dims){
      int nd=dims.size();
      int stride=1;
      index=0;
      for(int d=0;d<nd;d++){
	index = index+(int64_t)stride*coor[d];
	stride=stride*dims[d];
      }
    }
    template<class coor_t>
    static accelerator_inline void IndexFromCoor (const coor_t& coor,int &index,const coor_t &dims){
      int64_t index64;
      IndexFromCoor(coor,index64,dims);
      assert(index64<2*1024*1024*1024LL);
      index = (int) index64;
    }

    template<class coor_t>
    static inline void IndexFromCoorReversed (const coor_t& coor,int64_t &index,const coor_t &dims){
      int nd=dims.size();
      int stride=1;
      index=0;
      for(int d=nd-1;d>=0;d--){
	index = index+(int64_t)stride*coor[d];
	stride=stride*dims[d];
      }
    }
    template<class coor_t>
    static inline void IndexFromCoorReversed (const coor_t& coor,int &index,const coor_t &dims){
      int64_t index64;
      IndexFromCoorReversed(coor,index64,dims);
      if ( index64>=2*1024*1024*1024LL ){
	std::cout << " IndexFromCoorReversed " << coor<<" index " << index64<< " dims "<<dims<<std::endl;
      }
      assert(index64<2*1024*1024*1024LL);
      index = (int) index64;
    }
    template<class coor_t>
    static inline void CoorFromIndexReversed (coor_t& coor,int64_t index,const coor_t &dims){
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
