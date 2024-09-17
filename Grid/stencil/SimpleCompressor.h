#ifndef _STENCIL_SIMPLE_COMPRESSOR_H_
#define _STENCIL_SIMPLE_COMPRESSOR_H_

NAMESPACE_BEGIN(Grid);

class SimpleStencilParams{
public:
  Coordinate dirichlet;
  int partialDirichlet;
  SimpleStencilParams() { partialDirichlet = 0; };
};


// Compressors will inherit buffer management policies
// Standard comms buffer management
class FaceGatherSimple
{
public:
  static int PartialCompressionFactor(GridBase *grid) {return 1;};
  // Decompress is after merge so ok
  template<class vobj,class cobj,class compressor> 
  static void Gather_plane_simple (deviceVector<std::pair<int,int> >& table,
				   const Lattice<vobj> &rhs,
				   cobj *buffer,
				   compressor &compress,
				   int off,int so,int partial)
  {
    int num=table.size();
    std::pair<int,int> *table_v = & table[0];
    
    auto rhs_v = rhs.View(AcceleratorRead);
    accelerator_forNB( i,num, vobj::Nsimd(), {
	compress.Compress(buffer[off+table_v[i].first],rhs_v[so+table_v[i].second]);
    });
    rhs_v.ViewClose();
  }
  template<class vobj,class cobj,class compressor>
  static void Gather_plane_exchange(deviceVector<std::pair<int,int> >& table,const Lattice<vobj> &rhs,
				    std::vector<cobj *> pointers,int dimension,int plane,int cbmask,
				    compressor &compress,int type,int partial)
  {
    assert( (table.size()&0x1)==0);
    int num=table.size()/2;
    int so  = plane*rhs.Grid()->_ostride[dimension]; // base offset for start of plane
    
    auto rhs_v = rhs.View(AcceleratorRead);
    auto p0=&pointers[0][0];
    auto p1=&pointers[1][0];
    auto tp=&table[0];
    auto rhs_p = &rhs_v[0];
    accelerator_forNB(j, num, vobj::Nsimd(), {
	compress.CompressExchange(p0[j],p1[j],
				  rhs_p[so+tp[2*j  ].second],
				  rhs_p[so+tp[2*j+1].second],
				  type);
    });
    rhs_v.ViewClose();
  }

  template<class decompressor,class Decompression>
  static void DecompressFace(decompressor decompress,Decompression &dd)
  {
    auto kp = dd.kernel_p;
    auto mp = dd.mpi_p;
    accelerator_forNB(o,dd.buffer_size,1,{
      decompress.Decompress(kp[o],mp[o]);
    });
  }
  template<class decompressor,class Merger>
  static void MergeFace(decompressor decompress,Merger &mm)
  {
    auto mp = &mm.mpointer[0];
    auto vp0= &mm.vpointers[0][0];
    auto vp1= &mm.vpointers[1][0];
    auto type= mm.type;
    accelerator_forNB(o,mm.buffer_size/2,Merger::Nsimd,{
	decompress.Exchange(mp[2*o],mp[2*o+1],vp0[o],vp1[o],type);
    });
  }
};

////////////////////////////////////
// Wilson compressor will add alternate policies for Dirichlet
// and possibly partial Dirichlet for DWF
////////////////////////////////////

template<class vobj,class FaceGather>
class SimpleCompressorGather : public FaceGather {
public:
  void Point(int) {};
  accelerator_inline int  CommDatumSize(void) const { return sizeof(vobj); }
  accelerator_inline bool DecompressionStep(void) const { return false; }
  accelerator_inline void Compress(vobj &buf,const vobj &in) const {
    coalescedWrite(buf,coalescedRead(in));
  }
  accelerator_inline void Exchange(vobj &mp0,vobj &mp1,vobj &vp0,vobj &vp1,Integer type) const {
#ifdef GRID_SIMT
    exchangeSIMT(mp0,mp1,vp0,vp1,type);
#else
    exchange(mp0,mp1,vp0,vp1,type);
#endif
  }
  accelerator_inline void Decompress(vobj &out,vobj &in) const {  };
  accelerator_inline void CompressExchange(vobj &out0,vobj &out1,const vobj &in0,const vobj &in1,int type) const {
#ifdef GRID_SIMT
    exchangeSIMT(out0,out1,in0,in1,type);
#else
    exchange(out0,out1,in0,in1,type);
#endif
  }
  // For cshift. Cshift should drop compressor coupling altogether 
  // because I had to decouple the code from the Stencil anyway
  accelerator_inline vobj operator() (const vobj &arg) const {
    return arg;
  }
};

// Standard compressor never needs dirichlet.
//
// Get away with a local period wrap and rely on dirac operator to use a zero gauge link as it is faster
//
// Compressors that inherit Dirichlet and Non-dirichlet behaviour.
//
// Currently run-time behaviour through StencilParameters paramaters, p.dirichlet
// combined with the FaceGatherSimple behaviour

template <class vobj> using SimpleCompressor = SimpleCompressorGather<vobj,FaceGatherSimple>;
//template <class vobj> using SimpleCompressorDirichlet = SimpleCompressorGather<vobj,FaceGatherDirichlet>;

NAMESPACE_END(Grid);

#endif
