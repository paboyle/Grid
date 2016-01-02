    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/spin/TwoSpinor.h

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: Peter Boyle <peterboyle@Peters-MacBook-Pro-2.local>

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
#ifndef GRID_QCD_TWOSPIN_H
#define GRID_QCD_TWOSPIN_H
namespace Grid{
namespace QCD {

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  // Normalisation alert; the g5   project is 1/2(1+-G5) 
  //                      the xyzt projects are (1+-Gxyzt)
  //
  // * xyzt project
  //
  // This is because this is how the Wilson operator is normally written as
  // (m+4r) - \frac{1}{2} D_{hop}
  // and / or
  // 1 - \frac{1}{2 m+8r} D_{hop} = 1 - kappa D_{hop}
  //
  // Note that the free, critical hopping parameter kappa is then 1/8 th for r=1.
  //
  // However, the xyzt 2spin "projectors" are not really projectors because they do not
  // square to 1, however the ChiralProjector is a true projector.
  //
  // For this reason there is NO provision in Grid of a four spinor result from the
  // xyzt projectors. They are intended to be used only in combination with "reconstruct" in the
  // wilson dirac operator and closely related actions.
  //
  // I also do NOT provide lattice wide operators of these, since the dirac operator is best implemented
  // via Stencils and single site variants will be needed only for the cache friendly high perf dirac operator.
  //
  // * chiral project
  //
  // Both four spinor and two spinor result variants are provided.
  //
  // The four spinor project will be recursively provided to Lattice wide routines, and likely used in 
  // the domain wall and mobius implementations.
  //
  //////////////////////////////////////////////////////////////////////////////////////////////////////

  /* Gx
   *  0 0  0  i    [0]+-i[3]
   *  0 0  i  0    [1]+-i[2]
   *  0 -i 0  0
   * -i 0  0  0
   */

  // To fail is not to err (Cryptic clue: suggest to Google SFINAE ;) )
  template<class vtype,IfSpinor<iVector<vtype,Ns> > = 0> strong_inline void spProjXp (iVector<vtype,Nhs> &hspin,const iVector<vtype,Ns> &fspin)
    {
      hspin(0)=fspin(0)+timesI(fspin(3));
      hspin(1)=fspin(1)+timesI(fspin(2));
    }
  template<class vtype,IfSpinor<iVector<vtype,Ns> > = 0> strong_inline void spProjXm (iVector<vtype,Nhs> &hspin,const iVector<vtype,Ns> &fspin)
    {
      hspin(0)=fspin(0)-timesI(fspin(3));
      hspin(1)=fspin(1)-timesI(fspin(2));
    }

      //  0 0  0  -1  [0] -+ [3]
      //  0 0  1  0   [1] +- [2]
      //  0 1  0  0
      // -1 0  0  0
  template<class vtype,IfSpinor<iVector<vtype,Ns> > = 0> strong_inline void spProjYp (iVector<vtype,Nhs> &hspin,const iVector<vtype,Ns> &fspin)
    {
      hspin(0)=fspin(0)-fspin(3);
      hspin(1)=fspin(1)+fspin(2);
    }
  template<class vtype,IfSpinor<iVector<vtype,Ns> > = 0> strong_inline void spProjYm (iVector<vtype,Nhs> &hspin,const iVector<vtype,Ns> &fspin)
    {
      hspin(0)=fspin(0)+fspin(3);
      hspin(1)=fspin(1)-fspin(2);
    }
	    /*Gz
	     *  0 0  i  0   [0]+-i[2]
	     *  0 0  0 -i   [1]-+i[3]
	     * -i 0  0  0
	     *  0 i  0  0
	     */
  template<class vtype,IfSpinor<iVector<vtype,Ns> > = 0> strong_inline void spProjZp (iVector<vtype,Nhs> &hspin,const iVector<vtype,Ns> &fspin)
    {
      hspin(0)=fspin(0)+timesI(fspin(2));
      hspin(1)=fspin(1)-timesI(fspin(3));
    }
  template<class vtype,IfSpinor<iVector<vtype,Ns> > = 0> strong_inline void spProjZm (iVector<vtype,Nhs> &hspin,const iVector<vtype,Ns> &fspin)
    {
      //typename std::enable_if<matchGridTensorIndex<iVector<vtype,Ns>,SpinorIndex>::value,iVector<vtype,Ns> >::type *SFINAE;
      hspin(0)=fspin(0)-timesI(fspin(2));
      hspin(1)=fspin(1)+timesI(fspin(3));
    }
	    /*Gt
	     *  0 0  1  0 [0]+-[2]
	     *  0 0  0  1 [1]+-[3]
	     *  1 0  0  0
	     *  0 1  0  0
	     */
  template<class vtype,IfSpinor<iVector<vtype,Ns> > = 0> strong_inline void spProjTp (iVector<vtype,Nhs> &hspin,const iVector<vtype,Ns> &fspin)
    {
      //typename std::enable_if<matchGridTensorIndex<iVector<vtype,Ns>,SpinorIndex>::value,iVector<vtype,Ns> >::type *SFINAE;
      hspin(0)=fspin(0)+fspin(2);
      hspin(1)=fspin(1)+fspin(3);
    }
  template<class vtype,IfSpinor<iVector<vtype,Ns> > = 0> strong_inline void spProjTm (iVector<vtype,Nhs> &hspin,const iVector<vtype,Ns> &fspin)
    {
      //typename std::enable_if<matchGridTensorIndex<iVector<vtype,Ns>,SpinorIndex>::value,iVector<vtype,Ns> >::type *SFINAE;
      hspin(0)=fspin(0)-fspin(2);
      hspin(1)=fspin(1)-fspin(3);
    }
	    /*G5
	     *  1 0  0  0 
	     *  0 1  0  0 
	     *  0 0 -1  0
	     *  0 0  0 -1
	     */

  template<class vtype,IfSpinor<iVector<vtype,Ns> > = 0> strong_inline void spProj5p (iVector<vtype,Nhs> &hspin,const iVector<vtype,Ns> &fspin)
    {
      //typename std::enable_if<matchGridTensorIndex<iVector<vtype,Ns>,SpinorIndex>::value,iVector<vtype,Ns> >::type *SFINAE;
      hspin(0)=fspin(0);
      hspin(1)=fspin(1);
    }

  template<class vtype,IfSpinor<iVector<vtype,Ns> > = 0> strong_inline void spProj5m (iVector<vtype,Nhs> &hspin,const iVector<vtype,Ns> &fspin)
    {
      //typename std::enable_if<matchGridTensorIndex<iVector<vtype,Ns>,SpinorIndex>::value,iVector<vtype,Ns> >::type *SFINAE;
      hspin(0)=fspin(2);
      hspin(1)=fspin(3);
    }
  
  //  template<class vtype> strong_inline void fspProj5p (iVector<vtype,Ns> &rfspin,const iVector<vtype,Ns> &fspin)
  template<class vtype,IfSpinor<iVector<vtype,Ns> > = 0> strong_inline void spProj5p (iVector<vtype,Ns> &rfspin,const iVector<vtype,Ns> &fspin)
    {
      //typename std::enable_if<matchGridTensorIndex<iVector<vtype,Ns>,SpinorIndex>::value,iVector<vtype,Ns> >::type *SFINAE;
      rfspin(0)=fspin(0);
      rfspin(1)=fspin(1);
      rfspin(2)=zero;
      rfspin(3)=zero;
    }
  //  template<class vtype> strong_inline void fspProj5m (iVector<vtype,Ns> &rfspin,const iVector<vtype,Ns> &fspin)
  template<class vtype,IfSpinor<iVector<vtype,Ns> > = 0> strong_inline void spProj5m (iVector<vtype,Ns> &rfspin,const iVector<vtype,Ns> &fspin)
    {
      //typename std::enable_if<matchGridTensorIndex<iVector<vtype,Ns>,SpinorIndex>::value,iVector<vtype,Ns> >::type *SFINAE;
      rfspin(0)=zero;
      rfspin(1)=zero;
      rfspin(2)=fspin(2);
      rfspin(3)=fspin(3);
    }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Reconstruction routines to move back again to four spin
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /* Gx
   *  0 0  0  i    [0]+-i[3]
   *  0 0  i  0    [1]+-i[2]
   *  0 -i 0  0  -i[1]+-[2]   == -i ([0]+-i[3]) = -i (1)
   * -i 0  0  0  -i[0]+-[3]   == -i ([1]+-i[2]) = -i (0)
   */
  template<class vtype,IfSpinor<iVector<vtype,Ns> > = 0> strong_inline void spReconXp (iVector<vtype,Ns> &fspin,const iVector<vtype,Nhs> &hspin)
    {
      //typename std::enable_if<matchGridTensorIndex<iVector<vtype,Ns>,SpinorIndex>::value,iVector<vtype,Ns> >::type *SFINAE;
      fspin(0)=hspin(0);
      fspin(1)=hspin(1);
      fspin(2)=timesMinusI(hspin(1));
      fspin(3)=timesMinusI(hspin(0));
    }
  template<class vtype,IfSpinor<iVector<vtype,Ns> > = 0> strong_inline void spReconXm (iVector<vtype,Ns> &fspin,const iVector<vtype,Nhs> &hspin)
    {
      //typename std::enable_if<matchGridTensorIndex<iVector<vtype,Ns>,SpinorIndex>::value,iVector<vtype,Ns> >::type *SFINAE;
      fspin(0)=hspin(0);
      fspin(1)=hspin(1);
      fspin(2)=timesI(hspin(1));
      fspin(3)=timesI(hspin(0));
    }
  template<class vtype,IfSpinor<iVector<vtype,Ns> > = 0> strong_inline void accumReconXp (iVector<vtype,Ns> &fspin,const iVector<vtype,Nhs> &hspin)
    {
      //typename std::enable_if<matchGridTensorIndex<iVector<vtype,Ns>,SpinorIndex>::value,iVector<vtype,Ns> >::type *SFINAE;
      fspin(0)+=hspin(0);
      fspin(1)+=hspin(1);
      fspin(2)-=timesI(hspin(1));
      fspin(3)-=timesI(hspin(0));
    }
  template<class vtype,IfSpinor<iVector<vtype,Ns> > = 0> strong_inline void accumReconXm (iVector<vtype,Ns> &fspin,const iVector<vtype,Nhs> &hspin)
    {
      //typename std::enable_if<matchGridTensorIndex<iVector<vtype,Ns>,SpinorIndex>::value,iVector<vtype,Ns> >::type *SFINAE;
      fspin(0)+=hspin(0);
      fspin(1)+=hspin(1);
      fspin(2)+=timesI(hspin(1));
      fspin(3)+=timesI(hspin(0));
    }

      //  0 0  0  -1  [0] -+ [3]
      //  0 0  1  0   [1] +- [2]
      //  0 1  0  0              == 1(1)
      // -1 0  0  0              ==-1(0)

  template<class vtype,IfSpinor<iVector<vtype,Ns> > = 0> strong_inline void spReconYp (iVector<vtype,Ns> &fspin,const iVector<vtype,Nhs> &hspin)
    {
      //typename std::enable_if<matchGridTensorIndex<iVector<vtype,Ns>,SpinorIndex>::value,iVector<vtype,Ns> >::type *SFINAE;
      fspin(0)=hspin(0);
      fspin(1)=hspin(1);
      fspin(2)= hspin(1);
      fspin(3)=-hspin(0);//Unary minus?
    }
  template<class vtype,IfSpinor<iVector<vtype,Ns> > = 0> strong_inline void spReconYm (iVector<vtype,Ns> &fspin,const iVector<vtype,Nhs> &hspin)
    {
      //typename std::enable_if<matchGridTensorIndex<iVector<vtype,Ns>,SpinorIndex>::value,iVector<vtype,Ns> >::type *SFINAE;
      fspin(0)=hspin(0);
      fspin(1)=hspin(1);
      fspin(2)=-hspin(1);
      fspin(3)= hspin(0);
    }
  template<class vtype,IfSpinor<iVector<vtype,Ns> > = 0> strong_inline void accumReconYp (iVector<vtype,Ns> &fspin,const iVector<vtype,Nhs> &hspin)
    {
      //typename std::enable_if<matchGridTensorIndex<iVector<vtype,Ns>,SpinorIndex>::value,iVector<vtype,Ns> >::type *SFINAE;
      fspin(0)+=hspin(0);
      fspin(1)+=hspin(1);
      fspin(2)+=hspin(1);
      fspin(3)-=hspin(0);
    }
  template<class vtype,IfSpinor<iVector<vtype,Ns> > = 0> strong_inline void accumReconYm (iVector<vtype,Ns> &fspin,const iVector<vtype,Nhs> &hspin)
    {
      //typename std::enable_if<matchGridTensorIndex<iVector<vtype,Ns>,SpinorIndex>::value,iVector<vtype,Ns> >::type *SFINAE;
      fspin(0)+=hspin(0);
      fspin(1)+=hspin(1);
      fspin(2)-=hspin(1);
      fspin(3)+=hspin(0);
    }

	    /*Gz
	     *  0 0  i  0   [0]+-i[2]
	     *  0 0  0 -i   [1]-+i[3]
	     * -i 0  0  0     => -i (0)
	     *  0 i  0  0     =>  i (1)
	     */
  template<class vtype,IfSpinor<iVector<vtype,Ns> > = 0> strong_inline void spReconZp (iVector<vtype,Ns> &fspin,const iVector<vtype,Nhs> &hspin)
    {
      //typename std::enable_if<matchGridTensorIndex<iVector<vtype,Ns>,SpinorIndex>::value,iVector<vtype,Ns> >::type *SFINAE;
      fspin(0)=hspin(0);
      fspin(1)=hspin(1);
      fspin(2)=timesMinusI(hspin(0));
      fspin(3)=timesI(hspin(1));
    }
  template<class vtype,IfSpinor<iVector<vtype,Ns> > = 0> strong_inline void spReconZm (iVector<vtype,Ns> &fspin,const iVector<vtype,Nhs> &hspin)
    {
      //typename std::enable_if<matchGridTensorIndex<iVector<vtype,Ns>,SpinorIndex>::value,iVector<vtype,Ns> >::type *SFINAE;
      fspin(0)=hspin(0);
      fspin(1)=hspin(1);
      fspin(2)=     timesI(hspin(0));
      fspin(3)=timesMinusI(hspin(1));
    }
  template<class vtype,IfSpinor<iVector<vtype,Ns> > = 0> strong_inline void accumReconZp (iVector<vtype,Ns> &fspin,const iVector<vtype,Nhs> &hspin)
    {
      //typename std::enable_if<matchGridTensorIndex<iVector<vtype,Ns>,SpinorIndex>::value,iVector<vtype,Ns> >::type *SFINAE;
      fspin(0)+=hspin(0);
      fspin(1)+=hspin(1);
      fspin(2)-=timesI(hspin(0));
      fspin(3)+=timesI(hspin(1));
    }
  template<class vtype,IfSpinor<iVector<vtype,Ns> > = 0> strong_inline void accumReconZm (iVector<vtype,Ns> &fspin,const iVector<vtype,Nhs> &hspin)
    {
      //typename std::enable_if<matchGridTensorIndex<iVector<vtype,Ns>,SpinorIndex>::value,iVector<vtype,Ns> >::type *SFINAE;
      fspin(0)+=hspin(0);
      fspin(1)+=hspin(1);
      fspin(2)+=timesI(hspin(0));
      fspin(3)-=timesI(hspin(1));
    }
	    /*Gt
	     *  0 0  1  0 [0]+-[2]
	     *  0 0  0  1 [1]+-[3]
	     *  1 0  0  0    => (0)
	     *  0 1  0  0    => (1)
	     */
  template<class vtype,IfSpinor<iVector<vtype,Ns> > = 0> strong_inline void spReconTp (iVector<vtype,Ns> &fspin,const iVector<vtype,Nhs> &hspin)
    {
      //typename std::enable_if<matchGridTensorIndex<iVector<vtype,Ns>,SpinorIndex>::value,iVector<vtype,Ns> >::type *SFINAE;
      fspin(0)=hspin(0);
      fspin(1)=hspin(1);
      fspin(2)=hspin(0);
      fspin(3)=hspin(1);
    }
  template<class vtype,IfSpinor<iVector<vtype,Ns> > = 0> strong_inline void spReconTm (iVector<vtype,Ns> &fspin,const iVector<vtype,Nhs> &hspin)
    {
      //typename std::enable_if<matchGridTensorIndex<iVector<vtype,Ns>,SpinorIndex>::value,iVector<vtype,Ns> >::type *SFINAE;
      fspin(0)=hspin(0);
      fspin(1)=hspin(1);
      fspin(2)=-hspin(0);
      fspin(3)=-hspin(1);
    }
  template<class vtype,IfSpinor<iVector<vtype,Ns> > = 0> strong_inline void accumReconTp (iVector<vtype,Ns> &fspin,const iVector<vtype,Nhs> &hspin)
    {
      //typename std::enable_if<matchGridTensorIndex<iVector<vtype,Ns>,SpinorIndex>::value,iVector<vtype,Ns> >::type *SFINAE;
      fspin(0)+=hspin(0);
      fspin(1)+=hspin(1);
      fspin(2)+=hspin(0);
      fspin(3)+=hspin(1);
    }
  template<class vtype,IfSpinor<iVector<vtype,Ns> > = 0> strong_inline void accumReconTm (iVector<vtype,Ns> &fspin,const iVector<vtype,Nhs> &hspin)
    {
      //typename std::enable_if<matchGridTensorIndex<iVector<vtype,Ns>,SpinorIndex>::value,iVector<vtype,Ns> >::type *SFINAE;
      fspin(0)+=hspin(0);
      fspin(1)+=hspin(1);
      fspin(2)-=hspin(0);
      fspin(3)-=hspin(1);
    }
	    /*G5
	     *  1 0  0  0 
	     *  0 1  0  0 
	     *  0 0 -1  0
	     *  0 0  0 -1
	     */
  template<class vtype,IfSpinor<iVector<vtype,Ns> > = 0> strong_inline void spRecon5p (iVector<vtype,Ns> &fspin,const iVector<vtype,Nhs> &hspin)
    {
      //typename std::enable_if<matchGridTensorIndex<iVector<vtype,Ns>,SpinorIndex>::value,iVector<vtype,Ns> >::type *SFINAE;
      fspin(0)=hspin(0)+hspin(0); // add is lower latency than mul
      fspin(1)=hspin(1)+hspin(1); // probably no measurable diffence though
      fspin(2)=zero;
      fspin(3)=zero;
    }
  template<class vtype,IfSpinor<iVector<vtype,Ns> > = 0> strong_inline void spRecon5m (iVector<vtype,Ns> &fspin,const iVector<vtype,Nhs> &hspin)
    {
      //typename std::enable_if<matchGridTensorIndex<iVector<vtype,Ns>,SpinorIndex>::value,iVector<vtype,Ns> >::type *SFINAE;
      fspin(0)=zero;
      fspin(1)=zero;
      fspin(2)=hspin(0)+hspin(0);
      fspin(3)=hspin(1)+hspin(1);
    }
  template<class vtype,IfSpinor<iVector<vtype,Ns> > = 0> strong_inline void accumRecon5p (iVector<vtype,Ns> &fspin,const iVector<vtype,Nhs> &hspin)
    {
      //typename std::enable_if<matchGridTensorIndex<iVector<vtype,Ns>,SpinorIndex>::value,iVector<vtype,Ns> >::type *SFINAE;
      fspin(0)+=hspin(0)+hspin(0);
      fspin(1)+=hspin(1)+hspin(1);
    }
  template<class vtype,IfSpinor<iVector<vtype,Ns> > = 0> strong_inline void accumRecon5m (iVector<vtype,Ns> &fspin,const iVector<vtype,Nhs> &hspin)
    {
      //typename std::enable_if<matchGridTensorIndex<iVector<vtype,Ns>,SpinorIndex>::value,iVector<vtype,Ns> >::type *SFINAE;
      fspin(2)+=hspin(0)+hspin(0);
      fspin(3)+=hspin(1)+hspin(1);
    }

  //////////////////////////////////////////////////////////////////////////////////////////////
  // Recursively apply these until we hit the spin index
  //////////////////////////////////////////////////////////////////////////////////////////////

  //////////
  // Xp
  //////////
  template<class rtype,class vtype,int N,IfNotSpinor<iVector<vtype,N> > = 0> strong_inline void spProjXp (iVector<rtype,N> &hspin,const iVector<vtype,N> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iVector<vtype,N>,SpinorIndex>::notvalue,iVector<vtype,N> >::type *temp;
    for(int i=0;i<N;i++) {
      spProjXp(hspin._internal[i],fspin._internal[i]);
    }
  }
  template<class rtype,class vtype> strong_inline void spProjXp (iScalar<rtype> &hspin,const iScalar<vtype> &fspin)
  {
    spProjXp(hspin._internal,fspin._internal);
  }
  template<class rtype,class vtype,int N> strong_inline void spProjXp (iMatrix<rtype,N> &hspin,const iMatrix<vtype,N> &fspin)
  {
    for(int i=0;i<N;i++){ 
    for(int j=0;j<N;j++){
      spProjXp(hspin._internal[i][j],fspin._internal[i][j]);
    }}
  }


  template<class rtype,class vtype> strong_inline void spReconXp (iScalar<rtype> &hspin,const iScalar<vtype> &fspin)
  {
    spReconXp(hspin._internal,fspin._internal);
  }
  template<class rtype,class vtype,int N,IfNotSpinor<iVector<vtype,N> > = 0> strong_inline void spReconXp (iVector<rtype,N> &hspin,const iVector<vtype,N> &fspin)
  {
    for(int i=0;i<N;i++) {
      spReconXp(hspin._internal[i],fspin._internal[i]);
    }
  }
  template<class rtype,class vtype,int N> strong_inline void spReconXp (iMatrix<rtype,N> &hspin,const iMatrix<vtype,N> &fspin)
  {
    for(int i=0;i<N;i++){ 
    for(int j=0;j<N;j++){
      spReconXp(hspin._internal[i][j],fspin._internal[i][j]);
    }}
  }

  template<class rtype,class vtype> strong_inline void accumReconXp (iScalar<rtype> &hspin,const iScalar<vtype> &fspin)
  {
    accumReconXp(hspin._internal,fspin._internal);
  }
  template<class rtype,class vtype,int N,IfNotSpinor<iVector<vtype,N> > = 0> strong_inline void accumReconXp (iVector<rtype,N> &hspin,const iVector<vtype,N> &fspin)
  {
    for(int i=0;i<N;i++) {
      accumReconXp(hspin._internal[i],fspin._internal[i]);
    }
  }
  template<class rtype,class vtype,int N> strong_inline void accumReconXp (iMatrix<rtype,N> &hspin,const iMatrix<vtype,N> &fspin)
  {
    for(int i=0;i<N;i++){ 
    for(int j=0;j<N;j++){
      accumReconXp(hspin._internal[i][j],fspin._internal[i][j]);
    }}
  }



  ////////
  // Xm
  ////////
  template<class rtype,class vtype> strong_inline void spProjXm (iScalar<rtype> &hspin,const iScalar<vtype> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iScalar<vtype>,SpinorIndex>::notvalue,iScalar<vtype> >::type *temp;
    spProjXm(hspin._internal,fspin._internal);
  }
  template<class rtype,class vtype,int N,IfNotSpinor<iVector<vtype,N> > = 0> strong_inline void spProjXm (iVector<rtype,N> &hspin,const iVector<vtype,N> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iVector<vtype,N>,SpinorIndex>::notvalue,iVector<vtype,N> >::type *temp;
    for(int i=0;i<N;i++) {
      spProjXm(hspin._internal[i],fspin._internal[i]);
    }
  }
  template<class rtype,class vtype,int N> strong_inline void spProjXm (iMatrix<rtype,N> &hspin,const iMatrix<vtype,N> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iMatrix<vtype,N>,SpinorIndex>::notvalue,iMatrix<vtype,N> >::type *temp;
    for(int i=0;i<N;i++){ 
    for(int j=0;j<N;j++){
      spProjXm(hspin._internal[i][j],fspin._internal[i][j]);
    }}
  }


  template<class rtype,class vtype> strong_inline void spReconXm (iScalar<rtype> &hspin,const iScalar<vtype> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iScalar<vtype>,SpinorIndex>::notvalue,iScalar<vtype> >::type *temp;
    spReconXm(hspin._internal,fspin._internal);
  }
  template<class rtype,class vtype,int N,IfNotSpinor<iVector<vtype,N> > = 0> strong_inline void spReconXm (iVector<rtype,N> &hspin,const iVector<vtype,N> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iVector<vtype,N>,SpinorIndex>::notvalue,iVector<vtype,N> >::type *temp;
    for(int i=0;i<N;i++) {
      spReconXm(hspin._internal[i],fspin._internal[i]);
    }
  }
  template<class rtype,class vtype,int N> strong_inline void spReconXm (iMatrix<rtype,N> &hspin,const iMatrix<vtype,N> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iMatrix<vtype,N>,SpinorIndex>::notvalue,iMatrix<vtype,N> >::type *temp;
    for(int i=0;i<N;i++){ 
    for(int j=0;j<N;j++){
      spReconXm(hspin._internal[i][j],fspin._internal[i][j]);
    }}
  }

  template<class rtype,class vtype> strong_inline void accumReconXm (iScalar<rtype> &hspin,const iScalar<vtype> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iScalar<vtype>,SpinorIndex>::notvalue,iScalar<vtype> >::type *temp;
    accumReconXm(hspin._internal,fspin._internal);
  }
  template<class rtype,class vtype,int N,IfNotSpinor<iVector<vtype,N> > = 0> strong_inline void accumReconXm (iVector<rtype,N> &hspin,const iVector<vtype,N> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iVector<vtype,N>,SpinorIndex>::notvalue,iVector<vtype,N> >::type *temp;
    for(int i=0;i<N;i++) {
      accumReconXm(hspin._internal[i],fspin._internal[i]);
    }
  }
  template<class rtype,class vtype,int N> strong_inline void accumReconXm (iMatrix<rtype,N> &hspin,const iMatrix<vtype,N> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iMatrix<vtype,N>,SpinorIndex>::notvalue,iMatrix<vtype,N> >::type *temp;
    for(int i=0;i<N;i++){ 
    for(int j=0;j<N;j++){
      accumReconXm(hspin._internal[i][j],fspin._internal[i][j]);
    }}
  }



  ////////
  // Yp
  ////////
  template<class rtype,class vtype> strong_inline void spProjYp (iScalar<rtype> &hspin,const iScalar<vtype> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iScalar<vtype>,SpinorIndex>::notvalue,iScalar<vtype> >::type *temp;
    spProjYp(hspin._internal,fspin._internal);
  }
  template<class rtype,class vtype,int N,IfNotSpinor<iVector<vtype,N> > = 0> strong_inline void spProjYp (iVector<rtype,N> &hspin,const iVector<vtype,N> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iVector<vtype,N>,SpinorIndex>::notvalue,iVector<vtype,N> >::type *temp;
    for(int i=0;i<N;i++) {
      spProjYp(hspin._internal[i],fspin._internal[i]);
    }
  }
  template<class rtype,class vtype,int N> strong_inline void spProjYp (iMatrix<rtype,N> &hspin,const iMatrix<vtype,N> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iMatrix<vtype,N>,SpinorIndex>::notvalue,iMatrix<vtype,N> >::type *temp;
    for(int i=0;i<N;i++){ 
    for(int j=0;j<N;j++){
      spProjYp(hspin._internal[i][j],fspin._internal[i][j]);
    }}
  }


  template<class rtype,class vtype> strong_inline void spReconYp (iScalar<rtype> &hspin,const iScalar<vtype> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iScalar<vtype>,SpinorIndex>::notvalue,iScalar<vtype> >::type *temp;
    spReconYp(hspin._internal,fspin._internal);
  }
  template<class rtype,class vtype,int N,IfNotSpinor<iVector<vtype,N> > = 0> strong_inline void spReconYp (iVector<rtype,N> &hspin,const iVector<vtype,N> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iVector<vtype,N>,SpinorIndex>::notvalue,iVector<vtype,N> >::type *temp;
    for(int i=0;i<N;i++) {
      spReconYp(hspin._internal[i],fspin._internal[i]);
    }
  }
  template<class rtype,class vtype,int N> strong_inline void spReconYp (iMatrix<rtype,N> &hspin,const iMatrix<vtype,N> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iMatrix<vtype,N>,SpinorIndex>::notvalue,iMatrix<vtype,N> >::type *temp;
    for(int i=0;i<N;i++){ 
    for(int j=0;j<N;j++){
      spReconYp(hspin._internal[i][j],fspin._internal[i][j]);
    }}
  }

  template<class rtype,class vtype> strong_inline void accumReconYp (iScalar<rtype> &hspin,const iScalar<vtype> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iScalar<vtype>,SpinorIndex>::notvalue,iScalar<vtype> >::type *temp;
    accumReconYp(hspin._internal,fspin._internal);
  }
  template<class rtype,class vtype,int N,IfNotSpinor<iVector<vtype,N> > = 0> strong_inline void accumReconYp (iVector<rtype,N> &hspin,const iVector<vtype,N> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iVector<vtype,N>,SpinorIndex>::notvalue,iVector<vtype,N> >::type *temp;
    for(int i=0;i<N;i++) {
      accumReconYp(hspin._internal[i],fspin._internal[i]);
    }
  }
  template<class rtype,class vtype,int N> strong_inline void accumReconYp (iMatrix<rtype,N> &hspin,const iMatrix<vtype,N> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iMatrix<vtype,N>,SpinorIndex>::notvalue,iMatrix<vtype,N> >::type *temp;
    for(int i=0;i<N;i++){ 
    for(int j=0;j<N;j++){
      accumReconYp(hspin._internal[i][j],fspin._internal[i][j]);
    }}
  }


  ////////
  // Ym
  ////////
  template<class rtype,class vtype> strong_inline void spProjYm (iScalar<rtype> &hspin,const iScalar<vtype> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iScalar<vtype>,SpinorIndex>::notvalue,iScalar<vtype> >::type *temp;
    spProjYm(hspin._internal,fspin._internal);
  }
  template<class rtype,class vtype,int N,IfNotSpinor<iVector<vtype,N> > = 0> strong_inline void spProjYm (iVector<rtype,N> &hspin,const iVector<vtype,N> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iVector<vtype,N>,SpinorIndex>::notvalue,iVector<vtype,N> >::type *temp;
    for(int i=0;i<N;i++) {
      spProjYm(hspin._internal[i],fspin._internal[i]);
    }
  }
  template<class rtype,class vtype,int N> strong_inline void spProjYm (iMatrix<rtype,N> &hspin,const iMatrix<vtype,N> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iMatrix<vtype,N>,SpinorIndex>::notvalue,iMatrix<vtype,N> >::type *temp;
    for(int i=0;i<N;i++){ 
    for(int j=0;j<N;j++){
      spProjYm(hspin._internal[i][j],fspin._internal[i][j]);
    }}
  }


  template<class rtype,class vtype> strong_inline void spReconYm (iScalar<rtype> &hspin,const iScalar<vtype> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iScalar<vtype>,SpinorIndex>::notvalue,iScalar<vtype> >::type *temp;
    spReconYm(hspin._internal,fspin._internal);
  }
  template<class rtype,class vtype,int N,IfNotSpinor<iVector<vtype,N> > = 0> strong_inline void spReconYm (iVector<rtype,N> &hspin,const iVector<vtype,N> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iVector<vtype,N>,SpinorIndex>::notvalue,const iVector<vtype,N> >::type *temp;
    for(int i=0;i<N;i++) {
      spReconYm(hspin._internal[i],fspin._internal[i]);
    }
  }
  template<class rtype,class vtype,int N> strong_inline void spReconYm (iMatrix<rtype,N> &hspin,const iMatrix<vtype,N> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iMatrix<vtype,N>,SpinorIndex>::notvalue,iMatrix<vtype,N> >::type *temp;
    for(int i=0;i<N;i++){ 
    for(int j=0;j<N;j++){
      spReconYm(hspin._internal[i][j],fspin._internal[i][j]);
    }}
  }

  template<class rtype,class vtype> strong_inline void accumReconYm (iScalar<rtype> &hspin,const iScalar<vtype> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iScalar<vtype>,SpinorIndex>::notvalue,iScalar<vtype> >::type *temp;
    accumReconYm(hspin._internal,fspin._internal);
  }
  template<class rtype,class vtype,int N,IfNotSpinor<iVector<vtype,N> > = 0> strong_inline void accumReconYm (iVector<rtype,N> &hspin,const iVector<vtype,N> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iVector<vtype,N>,SpinorIndex>::notvalue,iVector<vtype,N> >::type *temp;
    for(int i=0;i<N;i++) {
      accumReconYm(hspin._internal[i],fspin._internal[i]);
    }
  }
  template<class rtype,class vtype,int N> strong_inline void accumReconYm (iMatrix<rtype,N> &hspin,const iMatrix<vtype,N> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iMatrix<vtype,N>,SpinorIndex>::notvalue,iMatrix<vtype,N> >::type *temp;
    for(int i=0;i<N;i++){ 
    for(int j=0;j<N;j++){
      accumReconYm(hspin._internal[i][j],fspin._internal[i][j]);
    }}
  }


  ////////
  // Zp
  ////////
  template<class rtype,class vtype> strong_inline void spProjZp (iScalar<rtype> &hspin,const iScalar<vtype> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iScalar<vtype>,SpinorIndex>::notvalue,iScalar<vtype> >::type *temp;
    spProjZp(hspin._internal,fspin._internal);
  }
  template<class rtype,class vtype,int N,IfNotSpinor<iVector<vtype,N> > = 0> strong_inline void spProjZp (iVector<rtype,N> &hspin,const iVector<vtype,N> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iVector<vtype,N>,SpinorIndex>::notvalue,iVector<vtype,N> >::type *temp;
    for(int i=0;i<N;i++) {
      spProjZp(hspin._internal[i],fspin._internal[i]);
    }
  }
  template<class rtype,class vtype,int N> strong_inline void spProjZp (iMatrix<rtype,N> &hspin,const iMatrix<vtype,N> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iMatrix<vtype,N>,SpinorIndex>::notvalue,iMatrix<vtype,N> >::type *temp;
    for(int i=0;i<N;i++){ 
    for(int j=0;j<N;j++){
      spProjZp(hspin._internal[i][j],fspin._internal[i][j]);
    }}
  }


  template<class rtype,class vtype> strong_inline void spReconZp (iScalar<rtype> &hspin,const iScalar<vtype> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iScalar<vtype>,SpinorIndex>::notvalue,iScalar<vtype> >::type *temp;
    spReconZp(hspin._internal,fspin._internal);
  }
  template<class rtype,class vtype,int N,IfNotSpinor<iVector<vtype,N> > = 0> strong_inline void spReconZp (iVector<rtype,N> &hspin,const iVector<vtype,N> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iVector<vtype,N>,SpinorIndex>::notvalue,iVector<vtype,N> >::type *temp;
    for(int i=0;i<N;i++) {
      spReconZp(hspin._internal[i],fspin._internal[i]);
    }
  }
  template<class rtype,class vtype,int N> strong_inline void spReconZp (iMatrix<rtype,N> &hspin,const iMatrix<vtype,N> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iMatrix<vtype,N>,SpinorIndex>::notvalue,iMatrix<vtype,N> >::type *temp;
    for(int i=0;i<N;i++){ 
    for(int j=0;j<N;j++){
      spReconZp(hspin._internal[i][j],fspin._internal[i][j]);
    }}
  }

  template<class rtype,class vtype> strong_inline void accumReconZp (iScalar<rtype> &hspin,const iScalar<vtype> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iScalar<vtype>,SpinorIndex>::notvalue,iScalar<vtype> >::type *temp;
    accumReconZp(hspin._internal,fspin._internal);
  }
  template<class rtype,class vtype,int N,IfNotSpinor<iVector<vtype,N> > = 0> strong_inline void accumReconZp (iVector<rtype,N> &hspin,const iVector<vtype,N> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iVector<vtype,N>,SpinorIndex>::notvalue,iVector<vtype,N> >::type *temp;
    for(int i=0;i<N;i++) {
      accumReconZp(hspin._internal[i],fspin._internal[i]);
    }
  }
  template<class rtype,class vtype,int N> strong_inline void accumReconZp (iMatrix<rtype,N> &hspin,const iMatrix<vtype,N> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iMatrix<vtype,N>,SpinorIndex>::notvalue,iMatrix<vtype,N> >::type *temp;
    for(int i=0;i<N;i++){ 
    for(int j=0;j<N;j++){
      accumReconZp(hspin._internal[i][j],fspin._internal[i][j]);
    }}
  }


  ////////
  // Zm
  ////////
  template<class rtype,class vtype> strong_inline void spProjZm (iScalar<rtype> &hspin,const iScalar<vtype> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iScalar<vtype>,SpinorIndex>::notvalue,iScalar<vtype> >::type *temp;
    spProjZm(hspin._internal,fspin._internal);
  }
  template<class rtype,class vtype,int N,IfNotSpinor<iVector<vtype,N> > = 0> strong_inline void spProjZm (iVector<rtype,N> &hspin,const iVector<vtype,N> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iVector<vtype,N>,SpinorIndex>::notvalue,iVector<vtype,N> >::type *temp;
    for(int i=0;i<N;i++) {
      spProjZm(hspin._internal[i],fspin._internal[i]);
    }
  }
  template<class rtype,class vtype,int N> strong_inline void spProjZm (iMatrix<rtype,N> &hspin,const iMatrix<vtype,N> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iMatrix<vtype,N>,SpinorIndex>::notvalue,iMatrix<vtype,N> >::type *temp;
    for(int i=0;i<N;i++){ 
    for(int j=0;j<N;j++){
      spProjZm(hspin._internal[i][j],fspin._internal[i][j]);
    }}
  }


  template<class rtype,class vtype> strong_inline void spReconZm (iScalar<rtype> &hspin,const iScalar<vtype> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iScalar<vtype>,SpinorIndex>::notvalue,iScalar<vtype> >::type *temp;
    spReconZm(hspin._internal,fspin._internal);
  }
  template<class rtype,class vtype,int N,IfNotSpinor<iVector<vtype,N> > = 0> strong_inline void spReconZm (iVector<rtype,N> &hspin,const iVector<vtype,N> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iVector<vtype,N>,SpinorIndex>::notvalue,iVector<vtype,N> >::type *temp;
    for(int i=0;i<N;i++) {
      spReconZm(hspin._internal[i],fspin._internal[i]);
    }
  }
  template<class rtype,class vtype,int N> strong_inline void spReconZm (iMatrix<rtype,N> &hspin,const iMatrix<vtype,N> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iMatrix<vtype,N>,SpinorIndex>::notvalue,iMatrix<vtype,N> >::type *temp;
    for(int i=0;i<N;i++){ 
    for(int j=0;j<N;j++){
      spReconZm(hspin._internal[i][j],fspin._internal[i][j]);
    }}
  }

  template<class rtype,class vtype> strong_inline void accumReconZm (iScalar<rtype> &hspin,const iScalar<vtype> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iScalar<vtype>,SpinorIndex>::notvalue,iScalar<vtype> >::type *temp;
    accumReconZm(hspin._internal,fspin._internal);
  }
  template<class rtype,class vtype,int N,IfNotSpinor<iVector<vtype,N> > = 0> strong_inline void accumReconZm (iVector<rtype,N> &hspin,const iVector<vtype,N> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iVector<vtype,N>,SpinorIndex>::notvalue,iVector<vtype,N> >::type *temp;
    for(int i=0;i<N;i++) {
      accumReconZm(hspin._internal[i],fspin._internal[i]);
    }
  }
  template<class rtype,class vtype,int N> strong_inline void accumReconZm (iMatrix<rtype,N> &hspin,const iMatrix<vtype,N> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iMatrix<vtype,N>,SpinorIndex>::notvalue,iMatrix<vtype,N> >::type *temp;
    for(int i=0;i<N;i++){ 
    for(int j=0;j<N;j++){
      accumReconZm(hspin._internal[i][j],fspin._internal[i][j]);
    }}
  }


  ////////
  // Tp
  ////////
  template<class rtype,class vtype> strong_inline void spProjTp (iScalar<rtype> &hspin,const iScalar<vtype> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iScalar<vtype>,SpinorIndex>::notvalue,iScalar<vtype> >::type *temp;
    spProjTp(hspin._internal,fspin._internal);
  }
  template<class rtype,class vtype,int N,IfNotSpinor<iVector<vtype,N> > = 0> strong_inline void spProjTp (iVector<rtype,N> &hspin,const iVector<vtype,N> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iVector<vtype,N>,SpinorIndex>::notvalue,iVector<vtype,N> >::type *temp;
    for(int i=0;i<N;i++) {
      spProjTp(hspin._internal[i],fspin._internal[i]);
    }
  }
  template<class rtype,class vtype,int N> strong_inline void spProjTp (iMatrix<rtype,N> &hspin,const iMatrix<vtype,N> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iMatrix<vtype,N>,SpinorIndex>::notvalue,iMatrix<vtype,N> >::type *temp;
    for(int i=0;i<N;i++){ 
    for(int j=0;j<N;j++){
      spProjTp(hspin._internal[i][j],fspin._internal[i][j]);
    }}
  }


  template<class rtype,class vtype> strong_inline void spReconTp (iScalar<rtype> &hspin,const iScalar<vtype> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iScalar<vtype>,SpinorIndex>::notvalue,iScalar<vtype> >::type *temp;
    spReconTp(hspin._internal,fspin._internal);
  }
  template<class rtype,class vtype,int N,IfNotSpinor<iVector<vtype,N> > = 0> strong_inline void spReconTp (iVector<rtype,N> &hspin,const iVector<vtype,N> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iVector<vtype,N>,SpinorIndex>::notvalue,iVector<vtype,N> >::type *temp;
    for(int i=0;i<N;i++) {
      spReconTp(hspin._internal[i],fspin._internal[i]);
    }
  }
  template<class rtype,class vtype,int N> strong_inline void spReconTp (iMatrix<rtype,N> &hspin,const iMatrix<vtype,N> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iMatrix<vtype,N>,SpinorIndex>::notvalue,iMatrix<vtype,N> >::type *temp;
    for(int i=0;i<N;i++){ 
    for(int j=0;j<N;j++){
      spReconTp(hspin._internal[i][j],fspin._internal[i][j]);
    }}
  }

  template<class rtype,class vtype> strong_inline void accumReconTp (iScalar<rtype> &hspin, iScalar<vtype> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iScalar<vtype>,SpinorIndex>::notvalue,iScalar<vtype> >::type *temp;
    accumReconTp(hspin._internal,fspin._internal);
  }
  template<class rtype,class vtype,int N,IfNotSpinor<iVector<vtype,N> > = 0> strong_inline void accumReconTp (iVector<rtype,N> &hspin, const iVector<vtype,N> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iVector<vtype,N>,SpinorIndex>::notvalue,iVector<vtype,N> >::type *temp;
    for(int i=0;i<N;i++) {
      accumReconTp(hspin._internal[i],fspin._internal[i]);
    }
  }
  template<class rtype,class vtype,int N> strong_inline void accumReconTp (iMatrix<rtype,N> &hspin, const iMatrix<vtype,N> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iMatrix<vtype,N>,SpinorIndex>::notvalue,iMatrix<vtype,N> >::type *temp;
    for(int i=0;i<N;i++){ 
    for(int j=0;j<N;j++){
      accumReconTp(hspin._internal[i][j],fspin._internal[i][j]);
    }}
  }


  ////////
  // Tm
  ////////
  template<class rtype,class vtype> strong_inline void spProjTm (iScalar<rtype> &hspin,const iScalar<vtype> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iScalar<vtype>,SpinorIndex>::notvalue,iScalar<vtype> >::type *temp;
    spProjTm(hspin._internal,fspin._internal);
  }
  template<class rtype,class vtype,int N,IfNotSpinor<iVector<vtype,N> > = 0> strong_inline void spProjTm (iVector<rtype,N> &hspin,const iVector<vtype,N> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iVector<vtype,N>,SpinorIndex>::notvalue,iVector<vtype,N> >::type *temp;
    for(int i=0;i<N;i++) {
      spProjTm(hspin._internal[i],fspin._internal[i]);
    }
  }
  template<class rtype,class vtype,int N> strong_inline void spProjTm (iMatrix<rtype,N> &hspin,const iMatrix<vtype,N> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iMatrix<vtype,N>,SpinorIndex>::notvalue,iMatrix<vtype,N> >::type *temp;
    for(int i=0;i<N;i++){ 
    for(int j=0;j<N;j++){
      spProjTm(hspin._internal[i][j],fspin._internal[i][j]);
    }}
  }


  template<class rtype,class vtype> strong_inline void spReconTm (iScalar<rtype> &hspin, const iScalar<vtype> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iScalar<vtype>,SpinorIndex>::notvalue,iScalar<vtype> >::type *temp;
    spReconTm(hspin._internal,fspin._internal);
  }
  template<class rtype,class vtype,int N,IfNotSpinor<iVector<vtype,N> > = 0> strong_inline void spReconTm (iVector<rtype,N> &hspin, const iVector<vtype,N> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iVector<vtype,N>,SpinorIndex>::notvalue,iVector<vtype,N> >::type *temp;
    for(int i=0;i<N;i++) {
      spReconTm(hspin._internal[i],fspin._internal[i]);
    }
  }
  template<class rtype,class vtype,int N> strong_inline void spReconTm (iMatrix<rtype,N> &hspin, const iMatrix<vtype,N> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iMatrix<vtype,N>,SpinorIndex>::notvalue,iMatrix<vtype,N> >::type *temp;
    for(int i=0;i<N;i++){ 
    for(int j=0;j<N;j++){
      spReconTm(hspin._internal[i][j],fspin._internal[i][j]);
    }}
  }

  template<class rtype,class vtype> strong_inline void accumReconTm (iScalar<rtype> &hspin, const iScalar<vtype> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iScalar<vtype>,SpinorIndex>::notvalue,iScalar<vtype> >::type *temp;
    accumReconTm(hspin._internal,fspin._internal);
  }
  template<class rtype,class vtype,int N,IfNotSpinor<iVector<vtype,N> > = 0> strong_inline void accumReconTm (iVector<rtype,N> &hspin, const iVector<vtype,N> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iVector<vtype,N>,SpinorIndex>::notvalue,iVector<vtype,N> >::type *temp;
    for(int i=0;i<N;i++) {
      accumReconTm(hspin._internal[i],fspin._internal[i]);
    }
  }
  template<class rtype,class vtype,int N> strong_inline void accumReconTm (iMatrix<rtype,N> &hspin, const iMatrix<vtype,N> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iMatrix<vtype,N>,SpinorIndex>::notvalue,iMatrix<vtype,N> >::type *temp;
    for(int i=0;i<N;i++){ 
    for(int j=0;j<N;j++){
      accumReconTm(hspin._internal[i][j],fspin._internal[i][j]);
    }}
  }


  ////////
  // 5p
  ////////
  template<class rtype,class vtype> strong_inline void spProj5p (iScalar<rtype> &hspin,const iScalar<vtype> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iScalar<vtype>,SpinorIndex>::notvalue,iScalar<vtype> >::type *temp;
    spProj5p(hspin._internal,fspin._internal);
  }
  template<class rtype,class vtype,int N,IfNotSpinor<iVector<vtype,N> > = 0> strong_inline void spProj5p (iVector<rtype,N> &hspin,const iVector<vtype,N> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iVector<vtype,N>,SpinorIndex>::notvalue,iVector<vtype,N> >::type *temp;
    for(int i=0;i<N;i++) {
      spProj5p(hspin._internal[i],fspin._internal[i]);
    }
  }
  template<class rtype,class vtype,int N> strong_inline void spProj5p (iMatrix<rtype,N> &hspin,const iMatrix<vtype,N> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iMatrix<vtype,N>,SpinorIndex>::notvalue,iMatrix<vtype,N> >::type *temp;
    for(int i=0;i<N;i++){ 
    for(int j=0;j<N;j++){
      spProj5p(hspin._internal[i][j],fspin._internal[i][j]);
    }}
  }

  template<class rtype,class vtype> strong_inline void spRecon5p (iScalar<rtype> &hspin,const iScalar<vtype> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iScalar<vtype>,SpinorIndex>::notvalue,iScalar<vtype> >::type *temp;
    spRecon5p(hspin._internal,fspin._internal);
  }
  template<class rtype,class vtype,int N,IfNotSpinor<iVector<vtype,N> > = 0> strong_inline void spRecon5p (iVector<rtype,N> &hspin,const iVector<vtype,N> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iVector<vtype,N>,SpinorIndex>::notvalue,iVector<vtype,N> >::type *temp;
    for(int i=0;i<N;i++) {
      spRecon5p(hspin._internal[i],fspin._internal[i]);
    }
  }
  template<class rtype,class vtype,int N> strong_inline void spRecon5p (iMatrix<rtype,N> &hspin,const iMatrix<vtype,N> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iMatrix<vtype,N>,SpinorIndex>::notvalue,iMatrix<vtype,N> >::type *temp;
    for(int i=0;i<N;i++){ 
    for(int j=0;j<N;j++){
      spRecon5p(hspin._internal[i][j],fspin._internal[i][j]);
    }}
  }

  template<class rtype,class vtype> strong_inline void accumRecon5p (iScalar<rtype> &hspin,const iScalar<vtype> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iScalar<vtype>,SpinorIndex>::notvalue,iScalar<vtype> >::type *temp;
    accumRecon5p(hspin._internal,fspin._internal);
  }
  template<class rtype,class vtype,int N,IfNotSpinor<iVector<vtype,N> > = 0> strong_inline void accumRecon5p (iVector<rtype,N> &hspin,const iVector<vtype,N> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iVector<vtype,N>,SpinorIndex>::notvalue,iVector<vtype,N> >::type *temp;
    for(int i=0;i<N;i++) {
      accumRecon5p(hspin._internal[i],fspin._internal[i]);
    }
  }
  template<class rtype,class vtype,int N> strong_inline void accumRecon5p (iMatrix<rtype,N> &hspin,const iMatrix<vtype,N> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iMatrix<vtype,N>,SpinorIndex>::notvalue,iMatrix<vtype,N> >::type *temp;
    for(int i=0;i<N;i++){ 
    for(int j=0;j<N;j++){
      accumRecon5p(hspin._internal[i][j],fspin._internal[i][j]);
    }}
  }

  // four spinor projectors for chiral proj
  //  template<class vtype> strong_inline void fspProj5p (iScalar<vtype> &hspin,const iScalar<vtype> &fspin)
  template<class vtype> strong_inline void spProj5p (iScalar<vtype> &hspin,const iScalar<vtype> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iScalar<vtype>,SpinorIndex>::notvalue,iScalar<vtype> >::type *temp;
    spProj5p(hspin._internal,fspin._internal);
  }
  //  template<class vtype,int N> strong_inline void fspProj5p (iVector<vtype,N> &hspin,iVector<vtype,N> &fspin)
  template<class vtype,int N,IfNotSpinor<iVector<vtype,N> > = 0> strong_inline void spProj5p (iVector<vtype,N> &hspin,const iVector<vtype,N> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iVector<vtype,N>,SpinorIndex>::notvalue,iVector<vtype,N> >::type *temp;
    for(int i=0;i<N;i++) {
      spProj5p(hspin._internal[i],fspin._internal[i]);
    }
  }
  //  template<class vtype,int N> strong_inline void fspProj5p (iMatrix<vtype,N> &hspin,iMatrix<vtype,N> &fspin)
  template<class vtype,int N> strong_inline void spProj5p (iMatrix<vtype,N> &hspin,const iMatrix<vtype,N> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iMatrix<vtype,N>,SpinorIndex>::notvalue,iMatrix<vtype,N> >::type *temp;
    for(int i=0;i<N;i++){ 
    for(int j=0;j<N;j++){
      spProj5p(hspin._internal[i][j],fspin._internal[i][j]);
    }}
  }


  ////////
  // 5m
  ////////

  template<class rtype,class vtype> strong_inline void spProj5m (iScalar<rtype> &hspin,const iScalar<vtype> &fspin)
  {
    spProj5m(hspin._internal,fspin._internal);
  }
  template<class rtype,class vtype,int N,IfNotSpinor<iVector<rtype,N> > = 0> strong_inline void spProj5m (iVector<rtype,N> &hspin,const iVector<vtype,N> &fspin)
  {
    for(int i=0;i<N;i++) {
      spProj5m(hspin._internal[i],fspin._internal[i]);
    }
  }
  template<class rtype,class vtype,int N> strong_inline void spProj5m (iMatrix<rtype,N> &hspin,const iMatrix<vtype,N> &fspin)
  {
    for(int i=0;i<N;i++){ 
    for(int j=0;j<N;j++){
      spProj5m(hspin._internal[i][j],fspin._internal[i][j]);
    }}
  }

  template<class rtype,class vtype> strong_inline void spRecon5m (iScalar<rtype> &hspin,const iScalar<vtype> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iScalar<vtype>,SpinorIndex>::notvalue,iScalar<vtype> >::type *temp;
    spRecon5m(hspin._internal,fspin._internal);
  }
  template<class rtype,class vtype,int N,IfNotSpinor<iVector<vtype,N> > = 0> strong_inline void spRecon5m (iVector<rtype,N> &hspin,const iVector<vtype,N> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iVector<vtype,N>,SpinorIndex>::notvalue,iVector<vtype,N> >::type *temp;
    for(int i=0;i<N;i++) {
      spRecon5m(hspin._internal[i],fspin._internal[i]);
    }
  }
  template<class rtype,class vtype,int N> strong_inline void spRecon5m (iMatrix<rtype,N> &hspin,const iMatrix<vtype,N> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iMatrix<vtype,N>,SpinorIndex>::notvalue,iMatrix<vtype,N> >::type *temp;
    for(int i=0;i<N;i++){ 
    for(int j=0;j<N;j++){
      spRecon5m(hspin._internal[i][j],fspin._internal[i][j]);
    }}
  }

  template<class rtype,class vtype> strong_inline void accumRecon5m (iScalar<rtype> &hspin,const iScalar<vtype> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iScalar<vtype>,SpinorIndex>::notvalue,iScalar<vtype> >::type *temp;
    accumRecon5m(hspin._internal,fspin._internal);
  }
  template<class rtype,class vtype,int N,IfNotSpinor<iVector<vtype,N> > = 0> strong_inline void accumRecon5m (iVector<rtype,N> &hspin,const iVector<vtype,N> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iVector<vtype,N>,SpinorIndex>::notvalue,iVector<vtype,N> >::type *temp;
    for(int i=0;i<N;i++) {
      accumRecon5m(hspin._internal[i],fspin._internal[i]);
    }
  }
  template<class rtype,class vtype,int N> strong_inline void accumRecon5m (iMatrix<rtype,N> &hspin,const iMatrix<vtype,N> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iMatrix<vtype,N>,SpinorIndex>::notvalue,iMatrix<vtype,N> >::type *temp;
    for(int i=0;i<N;i++){ 
    for(int j=0;j<N;j++){
      accumRecon5m(hspin._internal[i][j],fspin._internal[i][j]);
    }}
  }


  // four spinor projectors for chiral proj
  //  template<class vtype> strong_inline void fspProj5m (iScalar<vtype> &hspin,const iScalar<vtype> &fspin)
  template<class vtype> strong_inline void spProj5m (iScalar<vtype> &hspin,const iScalar<vtype> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iScalar<vtype>,SpinorIndex>::notvalue,iScalar<vtype> >::type *temp;
    spProj5m(hspin._internal,fspin._internal);
  }
  //  template<class vtype,int N> strong_inline void fspProj5m (iVector<vtype,N> &hspin,iVector<vtype,N> &fspin)
  template<class vtype,int N,IfNotSpinor<iVector<vtype,N> > = 0> strong_inline void spProj5m (iVector<vtype,N> &hspin,const iVector<vtype,N> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iVector<vtype,N>,SpinorIndex>::notvalue,iVector<vtype,N> >::type *temp;
    for(int i=0;i<N;i++) {
      spProj5m(hspin._internal[i],fspin._internal[i]);
    }
  }
  //  template<class vtype,int N> strong_inline void fspProj5m (iMatrix<vtype,N> &hspin,iMatrix<vtype,N> &fspin)
  template<class vtype,int N> strong_inline void spProj5m (iMatrix<vtype,N> &hspin,const iMatrix<vtype,N> &fspin)
  {
    //typename std::enable_if<matchGridTensorIndex<iMatrix<vtype,N>,SpinorIndex>::notvalue,iMatrix<vtype,N> >::type *temp;
    for(int i=0;i<N;i++){ 
    for(int j=0;j<N;j++){
      spProj5m(hspin._internal[i][j],fspin._internal[i][j]);
    }}
  }
}   //namespace QCD
} // Grid
#endif
