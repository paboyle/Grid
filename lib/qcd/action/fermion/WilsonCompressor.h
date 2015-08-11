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

    virtual SiteHalfSpinor operator () (const SiteSpinor &in,int dim,int plane,int osite,GridBase *grid) {
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

    SiteHalfSpinor operator () (const SiteSpinor &in,int dim,int plane,int osite,GridBase *grid)
    {
      std::vector<int> Gbcs({1,0,0,0});

      typedef typename SiteHalfSpinor::scalar_object scalar_object;

      const int Nsimd = grid->Nsimd();

      int ocb=grid->CheckerBoardFromOindex(osite);

      SiteHalfSpinor tmp = spinproject(in); // spin projected

      //////////////////////////////////////////////////////////////
      // Check whether we must flavour flip
      // do this if source is plane 0 on processor 0 in dimension dim
      //////////////////////////////////////////////////////////////
      if ( (this->mu==Xp)||(this->mu==Yp)||(this->mu==Zp)||(this->mu==Tp) ) {      

	if ( Gbcs[this->mu] && (grid->_processor_coor[dim] == 0) && (plane==0) && (ocb==0) ) {

	  std::vector<scalar_object>  flat(Nsimd);

	  extract(tmp,flat);

	  for(int i=0;i<Nsimd;i++) {
	    std::vector<int> coor;

	    grid->iCoorFromIindex(coor,i);

	    if ( coor[dim]==0 ) {
	      scalar_object stmp;
	      stmp(0) = flat[i](1);
	      stmp(1) = flat[i](0);
	      flat[i] = stmp;
	    }
	  }

	  merge(tmp,flat);

	} // coor & bc guard
      } // shift guard

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
