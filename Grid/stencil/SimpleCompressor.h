#ifndef _STENCIL_SIMPLE_COMPRESSOR_H_
#define _STENCIL_SIMPLE_COMPRESSOR_H_

NAMESPACE_BEGIN(Grid);

template<class vobj>
accelerator_inline void exchangeSIMT(vobj &mp0,vobj &mp1,const vobj &vp0,const vobj &vp1,Integer type)
{
    typedef decltype(coalescedRead(mp0)) sobj;
    unsigned int Nsimd = vobj::Nsimd();
    unsigned int mask = Nsimd >> (type + 1);
    int lane = acceleratorSIMTlane(Nsimd);
    int j0 = lane &(~mask); // inner coor zero
    int j1 = lane |(mask) ; // inner coor one
    const vobj *vpa = &vp0;
    const vobj *vpb = &vp1;
    const vobj *vp = (lane&mask) ? (vpb) : (vpa);
    auto sa = coalescedRead(vp[0],j0);
    auto sb = coalescedRead(vp[0],j1);
    coalescedWrite(mp0,sa);
    coalescedWrite(mp1,sb);
}

template<class vobj>
class SimpleCompressor {
public:
  void Point(int) {};
  accelerator_inline int  CommDatumSize(void) const { return sizeof(vobj); }
  accelerator_inline bool DecompressionStep(void) const { return false; }
  accelerator_inline void Compress(vobj &buf,const vobj &in) const {
    coalescedWrite(buf,coalescedRead(in));
  }
  accelerator_inline void Exchange(vobj *mp,vobj *vp0,vobj *vp1,Integer type,Integer o) const {
#ifdef GRID_SIMT
    exchangeSIMT(mp[2*o],mp[2*o+1],vp0[o],vp1[o],type);
#else
    exchange(mp[2*o],mp[2*o+1],vp0[o],vp1[o],type);
#endif
  }
  accelerator_inline void Decompress(vobj *out,vobj *in, int o) const { assert(0); }
  accelerator_inline void CompressExchange(vobj *out0,vobj *out1,const vobj *in,
					   int j,int k, int m,int type) const {
#ifdef GRID_SIMT
    exchangeSIMT(out0[j],out1[j],in[k],in[m],type);
#else
    exchange(out0[j],out1[j],in[k],in[m],type);
#endif
  }
  // For cshift. Cshift should drop compressor coupling altogether 
  // because I had to decouple the code from the Stencil anyway
  accelerator_inline vobj operator() (const vobj &arg) const {
    return arg;
  }
};

NAMESPACE_END(Grid);

#endif
