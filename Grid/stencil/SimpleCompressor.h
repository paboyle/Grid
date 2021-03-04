#ifndef _STENCIL_SIMPLE_COMPRESSOR_H_
#define _STENCIL_SIMPLE_COMPRESSOR_H_

NAMESPACE_BEGIN(Grid);

template<class vobj>
class SimpleCompressor {
public:
  void Point(int) {};
  accelerator_inline int  CommDatumSize(void) const { return sizeof(vobj); }
  accelerator_inline bool DecompressionStep(void) const { return false; }
  template<class cobj> accelerator_inline void Compress(cobj *buf,int o,const cobj &in) const { buf[o]=in; }
  accelerator_inline void Exchange(vobj *mp,vobj *vp0,vobj *vp1,Integer type,Integer o) const {
    exchange(mp[2*o],mp[2*o+1],vp0[o],vp1[o],type);
  }
  accelerator_inline void Decompress(vobj *out,vobj *in, int o) const { assert(0); }
  accelerator_inline void CompressExchange(vobj *out0,vobj *out1,const vobj *in,
			       int j,int k, int m,int type) const {
    exchange(out0[j],out1[j],in[k],in[m],type);
  }
  // For cshift. Cshift should drop compressor coupling altogether 
  // because I had to decouple the code from the Stencil anyway
  accelerator_inline vobj operator() (const vobj &arg) const {
    return arg;
  }
};

NAMESPACE_END(Grid);

#endif
