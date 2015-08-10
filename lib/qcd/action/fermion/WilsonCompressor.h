#ifndef  GRID_QCD_WILSON_COMPRESSOR_H
#define  GRID_QCD_WILSON_COMPRESSOR_H

namespace Grid {
namespace QCD {

  template<class SiteHalfSpinor,class SiteSpinor>
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

    virtual SiteHalfSpinor operator () (const SiteSpinor &in) {
      return spinproject(in);
    }

    SiteHalfSpinor spinproject(const SiteSpinor &in)
    {
      SiteHalfSpinor ret;
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


  template<class SiteHalfSpinor,class SiteSpinor>
    class GparityWilsonCompressor : public WilsonCompressor<SiteHalfSpinor,SiteSpinor>{
  public:

    GparityWilsonCompressor(int _dag) : WilsonCompressor<SiteHalfSpinor,SiteSpinor> (_dag){};

    SiteHalfSpinor operator () (const SiteSpinor &in)
    {
      SiteHalfSpinor tmp = spinproject(in);
      if( 0 ) tmp = flavourflip(tmp);
      return tmp;
    }

    SiteHalfSpinor flavourflip(const SiteHalfSpinor &in) {

      SiteHalfSpinor ret;
      for(int f=0;f<Ngp;f++){
	ret(0) = in(1);
	ret(1) = in(0);
      }
      return ret;
    }

  };


}} // namespace close
#endif
