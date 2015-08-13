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

  /*

  template<class SiteHalfSpinor,class SiteSpinor>
    class GparityWilsonCompressor : public WilsonCompressor<SiteHalfSpinor,SiteSpinor>{
  public:
    GparityWilsonCompressor(int _dag) : WilsonCompressor<SiteHalfSpinor,SiteSpinor> (_dag){};

    SiteHalfSpinor operator () (const SiteSpinor &in,int dim,int plane,int osite,GridBase *grid)
    {
      std::vector<int> Gbcs({1,0,0,0});

      typedef typename SiteHalfSpinor::scalar_object scalar_object;

      const int Nsimd = grid->Nsimd();

      int checkered=grid->CheckerBoarded(dim);
      int       ocb=grid->CheckerBoardFromOindex(osite);

      SiteHalfSpinor tmp = this->spinproject(in); // spin projected

      //////////////////////////////////////////////////////////////
      // Check whether we must flavour flip
      // do this if source is plane 0 on processor 0 in dimension dim
      //////////////////////////////////////////////////////////////
      int do_flip  = 0;
      int flipicoor= 0;
      if(Gbcs[this->mu]){

	std::cout << "Applying Gparity BC's in direction "<<this->mu<<std::endl;

	if ( (this->mu==Xp)||(this->mu==Yp)||(this->mu==Zp)||(this->mu==Tp) ) {      
	  if ( (grid->_processor_coor[dim] == 0) 
	       && (plane==0) 
	       && ( (!checkered)||(ocb==0) ) ) {
	    do_flip=1;
	    flipicoor=0;
	  }
	}
	if ( (this->mu==Xm)||(this->mu==Ym)||(this->mu==Zm)||(this->mu==Tm) ) {      
	  if ( (grid->_processor_coor[dim] == (grid->_processors[dim]-1) ) 
	    && (plane==grid->_rdimensions[dim]-1) 
            && ( (!checkered)||(ocb==1) ) ) {
	    do_flip=1;
	    flipicoor=grid->_simd_layout[dim]-1;
	  }
	}
      }
      if ( do_flip ) {

	std::cout << "Applying Gparity BC's in direction "<<this->mu<< " osite " << osite << " plane "<<plane <<std::endl;

	std::vector<scalar_object>  flat(Nsimd);
	std::vector<int> coor;

	extract(tmp,flat);
	for(int i=0;i<Nsimd;i++) {
	  grid->iCoorFromIindex(coor,i);
	  if ( coor[dim]==flipicoor ) {
	    scalar_object stmp;
	    stmp(0) = flat[i](1);
	    stmp(1) = flat[i](0);
	    flat[i] = stmp;
	  }
	}
	merge(tmp,flat);

      }

      return tmp;
    }

  };

  */

}} // namespace close
#endif
