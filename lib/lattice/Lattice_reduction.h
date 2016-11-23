    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/lattice/Lattice_reduction.h

    Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
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
#ifndef GRID_LATTICE_REDUCTION_H
#define GRID_LATTICE_REDUCTION_H

namespace Grid {

template <class T>
struct ReproducibilityState {
  typedef typename T::vector_type vector_type;
  unsigned int n_call;
  bool do_check;
  bool enable_reprocheck;
  std::vector<std::vector<vector_type, alignedAllocator<vector_type> > >
      th_states;

  void reset() {
    th_states.clear();
    do_check = false;
    enable_reprocheck = false;
    n_call = 0;
  }

  void reset_counter() { n_call = 0; }

  ReproducibilityState() { reset(); }
};

#ifdef GRID_WARN_SUBOPTIMAL
#warning "Optimisation alert all these reduction loops are NOT threaded "
#endif     

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Deterministic Reduction operations
    ////////////////////////////////////////////////////////////////////////////////////////////////////
  template<class vobj> inline RealD norm2(const Lattice<vobj> &arg){
        ReproducibilityState<vobj> repr;
        ComplexD nrm = innerProduct(arg, arg, repr);
        return std::real(nrm); 
    }

    template <class vobj>
    inline RealD norm2(const Lattice<vobj> &arg, ReproducibilityState<vobj>& rep) {
      ComplexD nrm = innerProduct(arg, arg, rep);
      return std::real(nrm);
    }

    template<class vobj>
    inline ComplexD innerProduct(const Lattice<vobj> &left,const Lattice<vobj> &right){
      ReproducibilityState<vobj> repr;
      return innerProduct(left, right, repr);
    } 
    


    template<class vobj>
    inline ComplexD innerProduct(const Lattice<vobj> &left,const Lattice<vobj> &right, ReproducibilityState<vobj>& repr) 
    {
      typedef typename vobj::scalar_type scalar_type;
      typedef typename vobj::vector_type vector_type;
      scalar_type  nrm;

      GridBase *grid = left._grid;

      std::vector<vector_type,alignedAllocator<vector_type> > sumarray(grid->SumArraySize());
      for(int i=0;i<grid->SumArraySize();i++){
        sumarray[i]=zero;
      }

      // accumulation done in the same precision ad vobj...
      // may need to froce higher precision
      PARALLEL_FOR_LOOP_STATIC //request statically scheduled threads for reproducibility
      for(int thr=0;thr<grid->SumArraySize();thr++){
        int nwork, mywork, myoff;
        GridThread::GetWork(left._grid->oSites(),thr,mywork,myoff);
        
        decltype(innerProduct(left._odata[0],right._odata[0])) vnrm = zero; // private to thread; sub summation
        for(int ss=myoff;ss<mywork+myoff; ss++){
          vnrm = vnrm + innerProduct(left._odata[ss],right._odata[ss]);// accumulate here in higher precision
        }
        sumarray[thr]=TensorRemove(vnrm) ;
      }
      

      ///////////////////////  Reproducibility section
      if (repr.enable_reprocheck) {
        if (repr.do_check) {
          //std::cout << GridLogMessage << "Checking thread state for inner product. Call n. " << repr.n_call << std::endl;
          for (int thread = 0; thread < sumarray.size(); thread++) {
            if (sumarray[thread] != repr.th_states[repr.n_call][thread]) {
              std::cout << GridLogMessage << "Reproducibility failure on node " << grid->ThisRank() << std::endl;
              std::cout << GridLogMessage << "Call: "<< repr.n_call << " Thread: " << thread << std::endl;
              std::cout << GridLogMessage << "Size of states: " << repr.th_states.size() << std::endl;
              std::cout << GridLogMessage << sumarray[thread] << std::endl;
              std::cout << GridLogMessage << repr.th_states[repr.n_call][thread] << std::endl;
              //exit(1);
            }
          }
          repr.n_call++;
        } else 
        {
          //std::cout << GridLogMessage << "Saving thread state for inner product. Call n. " << repr.n_call << std::endl;
          repr.th_states.resize(repr.n_call+1);
          repr.th_states[repr.n_call].resize(grid->SumArraySize());
          repr.th_states[repr.n_call] = sumarray;  // save threads state
          repr.n_call++;
        }
      }
      ////////////////////////////////////////////////////////


      vector_type vvnrm; vvnrm=zero;  // sum across threads
      for(int i=0;i<grid->SumArraySize();i++){
        vvnrm = vvnrm+sumarray[i];
      } 
      nrm = Reduce(vvnrm);// sum across simd
      right._grid->GlobalSum(nrm);
      return nrm;
    }

    template<class Op,class T1>
    inline auto sum(const LatticeUnaryExpression<Op,T1> & expr)
    ->typename decltype(expr.first.func(eval(0,std::get<0>(expr.second))))::scalar_object
    {
      return sum(closure(expr));
    }

    template<class Op,class T1,class T2>
    inline auto sum(const LatticeBinaryExpression<Op,T1,T2> & expr)
    ->typename decltype(expr.first.func(eval(0,std::get<0>(expr.second)),eval(0,std::get<1>(expr.second))))::scalar_object
    {
      return sum(closure(expr));
    }


    template<class Op,class T1,class T2,class T3>
    inline auto sum(const LatticeTrinaryExpression<Op,T1,T2,T3> & expr)
    ->typename decltype(expr.first.func(eval(0,std::get<0>(expr.second)),
     eval(0,std::get<1>(expr.second)),
     eval(0,std::get<2>(expr.second))
     ))::scalar_object
    {
      return sum(closure(expr));
    }

    template<class vobj>
    inline typename vobj::scalar_object sum(const Lattice<vobj> &arg){

      GridBase *grid=arg._grid;
      int Nsimd = grid->Nsimd();

      std::vector<vobj,alignedAllocator<vobj> > sumarray(grid->SumArraySize());
      for(int i=0;i<grid->SumArraySize();i++){
        sumarray[i]=zero;
      }

      PARALLEL_FOR_LOOP
      for(int thr=0;thr<grid->SumArraySize();thr++){
        int nwork, mywork, myoff;
        GridThread::GetWork(grid->oSites(),thr,mywork,myoff);

        vobj vvsum=zero;
        for(int ss=myoff;ss<mywork+myoff; ss++){
          vvsum = vvsum + arg._odata[ss];
        }
        sumarray[thr]=vvsum;
      }

      vobj vsum=zero;  // sum across threads
      for(int i=0;i<grid->SumArraySize();i++){
        vsum = vsum+sumarray[i];
      } 

      typedef typename vobj::scalar_object sobj;
      sobj ssum=zero;

      std::vector<sobj>               buf(Nsimd);
      extract(vsum,buf);

      for(int i=0;i<Nsimd;i++) ssum = ssum + buf[i];
        arg._grid->GlobalSum(ssum);

      return ssum;
    }



template<class vobj> inline void sliceSum(const Lattice<vobj> &Data,std::vector<typename vobj::scalar_object> &result,int orthogdim)
    {
      typedef typename vobj::scalar_object sobj;
      GridBase  *grid = Data._grid;
      assert(grid!=NULL);

  // FIXME
  // std::cout<<GridLogMessage<<"WARNING ! SliceSum is unthreaded "<<grid->SumArraySize()<<" threads "<<std::endl;

      const int    Nd = grid->_ndimension;
      const int Nsimd = grid->Nsimd();

      assert(orthogdim >= 0);
      assert(orthogdim < Nd);

      int fd=grid->_fdimensions[orthogdim];
      int ld=grid->_ldimensions[orthogdim];
      int rd=grid->_rdimensions[orthogdim];

  std::vector<vobj,alignedAllocator<vobj> > lvSum(rd); // will locally sum vectors first
  std::vector<sobj> lsSum(ld,zero); // sum across these down to scalars
  std::vector<sobj> extracted(Nsimd);     // splitting the SIMD

  result.resize(fd); // And then global sum to return the same vector to every node for IO to file
  for(int r=0;r<rd;r++){
    lvSum[r]=zero;
  }

  std::vector<int>  coor(Nd);  

  // sum over reduced dimension planes, breaking out orthog dir

  for(int ss=0;ss<grid->oSites();ss++){
    Lexicographic::CoorFromIndex(coor,ss,grid->_rdimensions);
    int r = coor[orthogdim];
    lvSum[r]=lvSum[r]+Data._odata[ss];
  }  

  // Sum across simd lanes in the plane, breaking out orthog dir.
  std::vector<int> icoor(Nd);

  for(int rt=0;rt<rd;rt++){

    extract(lvSum[rt],extracted);

    for(int idx=0;idx<Nsimd;idx++){

      grid->iCoorFromIindex(icoor,idx);

      int ldx =rt+icoor[orthogdim]*rd;

      lsSum[ldx]=lsSum[ldx]+extracted[idx];

    }
  }
  
  // sum over nodes.
  sobj gsum;
  for(int t=0;t<fd;t++){
    int pt = t/ld; // processor plane
    int lt = t%ld;
    if ( pt == grid->_processor_coor[orthogdim] ) {
      gsum=lsSum[lt];
    } else {
      gsum=zero;
    }

    grid->GlobalSum(gsum);

    result[t]=gsum;
  }

}


}
#endif

