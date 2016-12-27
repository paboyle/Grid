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

    typedef CartesianCommunicator::CommsRequest_t CommsRequest_t;

    WilsonStencil(GridBase *grid,
		int npoints,
		int checkerboard,
		const std::vector<int> &directions,
		const std::vector<int> &distances)  : CartesianStencil<vobj,cobj> (grid,npoints,checkerboard,directions,distances) 
      {    };


    template < class compressor>
    void HaloExchangeOpt(const Lattice<vobj> &source,compressor &compress) 
    {
      std::vector<std::vector<CommsRequest_t> > reqs;
      this->Mergers.resize(0); 
      this->Packets.resize(0);
      this->HaloGatherOpt(source,compress);
      this->CommunicateBegin(reqs);
      this->CommunicateComplete(reqs);
      this->CommsMerge(); // spins
      this->calls++;
    }


    template < class compressor>
    void HaloGatherOpt(const Lattice<vobj> &source,compressor &compress)
    {
      int face_idx=0;

      // conformable(source._grid,_grid);
      assert(source._grid==this->_grid);
      this->halogtime-=usecond();
      
      this->u_comm_offset=0;
      
      int dag = compress.dag;
      
      WilsonXpCompressor<cobj,vobj> XpCompress; 
      WilsonYpCompressor<cobj,vobj> YpCompress; 
      WilsonZpCompressor<cobj,vobj> ZpCompress; 
      WilsonTpCompressor<cobj,vobj> TpCompress;
      WilsonXmCompressor<cobj,vobj> XmCompress;
      WilsonYmCompressor<cobj,vobj> YmCompress;
      WilsonZmCompressor<cobj,vobj> ZmCompress;
      WilsonTmCompressor<cobj,vobj> TmCompress;

      // Gather all comms buffers
      //    for(int point = 0 ; point < _npoints; point++) {
      //      compress.Point(point);
      //      HaloGatherDir(source,compress,point,face_idx);
      //    }
      if ( dag ) { 
	this->HaloGatherDir(source,XpCompress,Xp,face_idx);
	this->HaloGatherDir(source,YpCompress,Yp,face_idx);
	this->HaloGatherDir(source,ZpCompress,Zp,face_idx);
	this->HaloGatherDir(source,TpCompress,Tp,face_idx);
	this->HaloGatherDir(source,XmCompress,Xm,face_idx);
	this->HaloGatherDir(source,YmCompress,Ym,face_idx);
	this->HaloGatherDir(source,ZmCompress,Zm,face_idx);
	this->HaloGatherDir(source,TmCompress,Tm,face_idx);
      } else {
	this->HaloGatherDir(source,XmCompress,Xp,face_idx);
	this->HaloGatherDir(source,YmCompress,Yp,face_idx);
	this->HaloGatherDir(source,ZmCompress,Zp,face_idx);
	this->HaloGatherDir(source,TmCompress,Tp,face_idx);
	this->HaloGatherDir(source,XpCompress,Xm,face_idx);
	this->HaloGatherDir(source,YpCompress,Ym,face_idx);
	this->HaloGatherDir(source,ZpCompress,Zm,face_idx);
	this->HaloGatherDir(source,TpCompress,Tm,face_idx);
      }
      this->face_table_computed=1;
      assert(this->u_comm_offset==this->_unified_buffer_size);
      this->halogtime+=usecond();
    }

  };


}} // namespace close
#endif
