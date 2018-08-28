#ifndef _STENCIL_SIMPLE_COMPRESSOR_H_
#define _STENCIL_SIMPLE_COMPRESSOR_H_

namespace Grid {

template<class vobj>
class SimpleCompressor {
public:
  void Point(int) {};
  inline int  CommDatumSize(void) { return sizeof(vobj); }
  inline bool DecompressionStep(void) { return false; }
  inline void Compress(vobj *buf,int o,const vobj &in) { buf[o]=in; }
  inline void Exchange(vobj *mp,vobj *vp0,vobj *vp1,Integer type,Integer o){
    exchange(mp[2*o],mp[2*o+1],vp0[o],vp1[o],type);
  }
  inline void Decompress(vobj *out,vobj *in, int o){ assert(0); }
  inline void CompressExchange(vobj *out0,vobj *out1,const vobj *in,
			       int j,int k, int m,int type){
    exchange(out0[j],out1[j],in[k],in[m],type);
  }
  // For cshift. Cshift should drop compressor coupling altogether 
  // because I had to decouple the code from the Stencil anyway
  inline vobj operator() (const vobj &arg) {
    return arg;
  }
};

}
#endif
