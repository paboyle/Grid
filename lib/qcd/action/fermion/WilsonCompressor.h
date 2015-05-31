#ifndef  GRID_QCD_WILSON_COMPRESSOR_H
#define  GRID_QCD_WILSON_COMPRESSOR_H

namespace Grid {
namespace QCD {

  class WilsonCompressor {
  public:
    int mu;
    int dag;

    WilsonCompressor(int _dag){
      mu=0;
      dag=_dag;
      assert((dag==0)||(dag==1));
    }
    void Point(int p) { 
      mu=p;
    };

    vHalfSpinColourVector operator () (const vSpinColourVector &in)
    {
      vHalfSpinColourVector ret;
      int mudag=mu;
      if (dag) {
	mudag=(mu+Nd)%(2*Nd);
      }
      switch(mudag) {
      case Xp:
	spProjXp(ret,in);
	break;
      case Yp:
	spProjYp(ret,in);
	break;
      case Zp:
	spProjZp(ret,in);
	break;
      case Tp:
	spProjTp(ret,in);
	break;
      case Xm:
	spProjXm(ret,in);
	break;
      case Ym:
	spProjYm(ret,in);
	break;
      case Zm:
	spProjZm(ret,in);
	break;
      case Tm:
	spProjTm(ret,in);
	break;
      default: 
	assert(0);
	break;
      }
      return ret;
    }
  };
}} // namespace close
#endif
