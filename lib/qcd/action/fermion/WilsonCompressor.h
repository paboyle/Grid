    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/action/fermion/WilsonCompressor.h

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: Peter Boyle <peterboyle@Peters-MacBook-Pro-2.local>
Author: paboyle <paboyle@ph.ed.ac.uk>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

    See the full license in the file "LICENSE" in the top level distribution directory
    *************************************************************************************/
    /*  END LEGAL */
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

    inline SiteHalfSpinor operator () (const SiteSpinor &in) {
      SiteHalfSpinor ret;
      int mudag=mu;
      if (!dag) {
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

  /////////////////////////
  // optimised versions
  /////////////////////////

  template<class SiteHalfSpinor,class SiteSpinor>
  class WilsonXpCompressor {
  public:
    inline SiteHalfSpinor operator () (const SiteSpinor &in) {
      SiteHalfSpinor ret;
      spProjXp(ret,in);
      return ret;
    }
  };
  template<class SiteHalfSpinor,class SiteSpinor>
  class WilsonYpCompressor {
  public:
    inline SiteHalfSpinor operator () (const SiteSpinor &in) {
      SiteHalfSpinor ret;
      spProjYp(ret,in);
      return ret;
    }
  };
  template<class SiteHalfSpinor,class SiteSpinor>
  class WilsonZpCompressor {
  public:
    inline SiteHalfSpinor operator () (const SiteSpinor &in) {
      SiteHalfSpinor ret;
      spProjZp(ret,in);
      return ret;
    }
  };
  template<class SiteHalfSpinor,class SiteSpinor>
  class WilsonTpCompressor {
  public:
    inline SiteHalfSpinor operator () (const SiteSpinor &in) {
      SiteHalfSpinor ret;
      spProjTp(ret,in);
      return ret;
    }
  };

  template<class SiteHalfSpinor,class SiteSpinor>
  class WilsonXmCompressor {
  public:
    inline SiteHalfSpinor operator () (const SiteSpinor &in) {
      SiteHalfSpinor ret;
      spProjXm(ret,in);
      return ret;
    }
  };
  template<class SiteHalfSpinor,class SiteSpinor>
  class WilsonYmCompressor {
  public:
    inline SiteHalfSpinor operator () (const SiteSpinor &in) {
      SiteHalfSpinor ret;
      spProjYm(ret,in);
      return ret;
    }
  };
  template<class SiteHalfSpinor,class SiteSpinor>
  class WilsonZmCompressor {
  public:
    inline SiteHalfSpinor operator () (const SiteSpinor &in) {
      SiteHalfSpinor ret;
      spProjZm(ret,in);
      return ret;
    }
  };
  template<class SiteHalfSpinor,class SiteSpinor>
  class WilsonTmCompressor {
  public:
    inline SiteHalfSpinor operator () (const SiteSpinor &in) {
      SiteHalfSpinor ret;
      spProjTm(ret,in);
      return ret;
    }
  };

    // Fast comms buffer manipulation which should inline right through (avoid direction
    // dependent logic that prevents inlining
  template<class vobj,class cobj>
  class WilsonStencil : public CartesianStencil<vobj,cobj> {
  public:

    WilsonStencil(GridBase *grid,
		int npoints,
		int checkerboard,
		const std::vector<int> &directions,
		const std::vector<int> &distances)  : CartesianStencil<vobj,cobj> (grid,npoints,checkerboard,directions,distances) 
      {    };

    template < class compressor>
    std::thread HaloExchangeOptBegin(const Lattice<vobj> &source,compressor &compress) {
      this->Mergers.resize(0); 
      this->Packets.resize(0);
      this->HaloGatherOpt(source,compress);
      return std::thread([&] { this->Communicate(); });
    }

    template < class compressor>
    void HaloExchangeOpt(const Lattice<vobj> &source,compressor &compress) 
    {
      auto thr = this->HaloExchangeOptBegin(source,compress);
      this->HaloExchangeOptComplete(thr);
    }

    void HaloExchangeOptComplete(std::thread &thr) 
    {
	this->CommsMerge(); // spins
	this->jointime-=usecond();
	thr.join();
	this->jointime+=usecond();
    }

    template < class compressor>
    void HaloGatherOpt(const Lattice<vobj> &source,compressor &compress)
    {
	// conformable(source._grid,_grid);
	assert(source._grid==this->_grid);
	this->halogtime-=usecond();

	assert (this->comm_buf.size() == this->_unified_buffer_size );
	this->u_comm_offset=0;

	int dag = compress.dag;
	static std::vector<int> dirs(Nd*2);
	for(int mu=0;mu<Nd;mu++){
	  if ( dag ) {
	    dirs[mu]  =mu;
	    dirs[mu+4]=mu+Nd;
	  } else { 
	    dirs[mu]  =mu+Nd;
	    dirs[mu+Nd]=mu;
	  }
	}


	WilsonXpCompressor<cobj,vobj> XpCompress;
	this->HaloGatherDir(source,XpCompress,dirs[0]);

	WilsonYpCompressor<cobj,vobj> YpCompress;
	this->HaloGatherDir(source,YpCompress,dirs[1]);

	WilsonZpCompressor<cobj,vobj> ZpCompress;
	this->HaloGatherDir(source,ZpCompress,dirs[2]);

	WilsonTpCompressor<cobj,vobj> TpCompress;
	this->HaloGatherDir(source,TpCompress,dirs[3]);

	WilsonXmCompressor<cobj,vobj> XmCompress;
	this->HaloGatherDir(source,XmCompress,dirs[4]);

	WilsonYmCompressor<cobj,vobj> YmCompress;
	this->HaloGatherDir(source,YmCompress,dirs[5]);

	WilsonZmCompressor<cobj,vobj> ZmCompress;
	this->HaloGatherDir(source,ZmCompress,dirs[6]);

	WilsonTmCompressor<cobj,vobj> TmCompress;
	this->HaloGatherDir(source,TmCompress,dirs[7]);

	assert(this->u_comm_offset==this->_unified_buffer_size);
	this->halogtime+=usecond();
      }

  };


}} // namespace close
#endif
