/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/lattice/Lattice_base.h

Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: paboyle <paboyle@ph.ed.ac.uk>
Author: Christoph Lehner <christoph@lhnr.de>

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

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
			   /*  END LEGAL */

#pragma once 

#define STREAMING_STORES

NAMESPACE_BEGIN(Grid);

extern int GridCshiftPermuteMap[4][16];

/////////////////////////////////////////////////////////////////////////////////////////
// The real lattice class, with normal copy and assignment semantics.
// This contains extra (host resident) grid pointer data that may be accessed by host code
/////////////////////////////////////////////////////////////////////////////////////////
template<class vobj>
class Lattice : public LatticeAccelerator<vobj>
{
public:
  GridBase *Grid(void) const { return this->_grid; }
  ///////////////////////////////////////////////////
  // Member types
  ///////////////////////////////////////////////////
  typedef typename vobj::scalar_type scalar_type;
  typedef typename vobj::vector_type vector_type;
  typedef typename vobj::scalar_object scalar_object;
  typedef vobj vector_object;

private:
  void dealloc(void)
  {
    if( this->_odata_size ) {
      alignedAllocator<vobj> alloc;
      alloc.deallocate(this->_odata,this->_odata_size);
      this->_odata=nullptr;
      this->_odata_size=0;
    }
  }
  void resize(uint64_t size)
  {
    if ( this->_odata_size != size ) {
      alignedAllocator<vobj> alloc;

      dealloc();
      
      this->_odata_size = size;
      if ( size )
	this->_odata      = alloc.allocate(this->_odata_size);
      else 
	this->_odata      = nullptr;
    }
  }
public:

  /////////////////////////////////////////////////////////////////////////////////
  // Can use to make accelerator dirty without copy from host ; useful for temporaries "dont care" prev contents
  /////////////////////////////////////////////////////////////////////////////////
  void SetViewMode(ViewMode mode) {
    LatticeView<vobj> accessor(*( (LatticeAccelerator<vobj> *) this),mode);
    accessor.ViewClose();
  }

  // Helper function to print the state of this object in the AccCache
  void PrintCacheState(void)
  {
    MemoryManager::PrintState(this->_odata);
  }

  /////////////////////////////////////////////////////////////////////////////////
  // Return a view object that may be dereferenced in site loops.
  // The view is trivially copy constructible and may be copied to an accelerator device
  // in device lambdas
  /////////////////////////////////////////////////////////////////////////////////

  LatticeView<vobj> View (ViewMode mode) const 
  {
    LatticeView<vobj> accessor(*( (LatticeAccelerator<vobj> *) this),mode);
    return accessor;
  }

  ~Lattice() { 
    if ( this->_odata_size ) {
      dealloc();
    }
   }
  ////////////////////////////////////////////////////////////////////////////////
  // Expression Template closure support
  ////////////////////////////////////////////////////////////////////////////////
  template <typename Op, typename T1> inline Lattice<vobj> & operator=(const LatticeUnaryExpression<Op,T1> &expr)
  {
    GRID_TRACE("ExpressionTemplateEval");
    GridBase *egrid(nullptr);
    GridFromExpression(egrid,expr);
    assert(egrid!=nullptr);
    conformable(this->_grid,egrid);

    int cb=-1;
    CBFromExpression(cb,expr);
    assert( (cb==Odd) || (cb==Even));
    this->checkerboard=cb;
    
    auto exprCopy = expr;
    ExpressionViewOpen(exprCopy);
    auto me  = View(AcceleratorWriteDiscard);
    accelerator_for(ss,me.size(),vobj::Nsimd(),{
      auto tmp = eval(ss,exprCopy);
      coalescedWrite(me[ss],tmp);
    });
    me.ViewClose();
    ExpressionViewClose(exprCopy);
    return *this;
  }
  template <typename Op, typename T1,typename T2> inline Lattice<vobj> & operator=(const LatticeBinaryExpression<Op,T1,T2> &expr)
  {
    GRID_TRACE("ExpressionTemplateEval");
    GridBase *egrid(nullptr);
    GridFromExpression(egrid,expr);
    assert(egrid!=nullptr);
    conformable(this->_grid,egrid);

    int cb=-1;
    CBFromExpression(cb,expr);
    assert( (cb==Odd) || (cb==Even));
    this->checkerboard=cb;

    auto exprCopy = expr;
    ExpressionViewOpen(exprCopy);
    auto me  = View(AcceleratorWriteDiscard);
    accelerator_for(ss,me.size(),vobj::Nsimd(),{
      auto tmp = eval(ss,exprCopy);
      coalescedWrite(me[ss],tmp);
    });
    me.ViewClose();
    ExpressionViewClose(exprCopy);
    return *this;
  }
  template <typename Op, typename T1,typename T2,typename T3> inline Lattice<vobj> & operator=(const LatticeTrinaryExpression<Op,T1,T2,T3> &expr)
  {
    GRID_TRACE("ExpressionTemplateEval");
    GridBase *egrid(nullptr);
    GridFromExpression(egrid,expr);
    assert(egrid!=nullptr);
    conformable(this->_grid,egrid);

    int cb=-1;
    CBFromExpression(cb,expr);
    assert( (cb==Odd) || (cb==Even));
    this->checkerboard=cb;
    auto exprCopy = expr;
    ExpressionViewOpen(exprCopy);
    auto me  = View(AcceleratorWriteDiscard);
    accelerator_for(ss,me.size(),vobj::Nsimd(),{
      auto tmp = eval(ss,exprCopy);
      coalescedWrite(me[ss],tmp);
    });
    me.ViewClose();
    ExpressionViewClose(exprCopy);
    return *this;
  }
  //GridFromExpression is tricky to do
  template<class Op,class T1>
  Lattice(const LatticeUnaryExpression<Op,T1> & expr) {
    this->_grid = nullptr;
    GridFromExpression(this->_grid,expr);
    assert(this->_grid!=nullptr);

    int cb=-1;
    CBFromExpression(cb,expr);
    assert( (cb==Odd) || (cb==Even));
    this->checkerboard=cb;

    resize(this->_grid->oSites());

    *this = expr;
  }
  template<class Op,class T1, class T2>
  Lattice(const LatticeBinaryExpression<Op,T1,T2> & expr) {
    this->_grid = nullptr;
    GridFromExpression(this->_grid,expr);
    assert(this->_grid!=nullptr);

    int cb=-1;
    CBFromExpression(cb,expr);
    assert( (cb==Odd) || (cb==Even));
    this->checkerboard=cb;

    resize(this->_grid->oSites());

    *this = expr;
  }
  template<class Op,class T1, class T2, class T3>
  Lattice(const LatticeTrinaryExpression<Op,T1,T2,T3> & expr) {
    this->_grid = nullptr;
    GridFromExpression(this->_grid,expr);
    assert(this->_grid!=nullptr);

    int cb=-1;
    CBFromExpression(cb,expr);
    assert( (cb==Odd) || (cb==Even));
    this->checkerboard=cb;

    resize(this->_grid->oSites());

    *this = expr;
  }

  template<class sobj> inline Lattice<vobj> & operator = (const sobj & r){
    vobj vtmp;
    vtmp = r;
#if 1
    auto me  = View(CpuWrite);
    thread_for(ss,me.size(),{
       me[ss]= r;
      });
#else    
    auto me  = View(AcceleratorWrite);
    accelerator_for(ss,me.size(),vobj::Nsimd(),{
	auto stmp=coalescedRead(vtmp);
	coalescedWrite(me[ss],stmp);
    });
#endif    
    me.ViewClose();
    return *this;
  }

  //////////////////////////////////////////////////////////////////
  // Follow rule of five, with Constructor requires "grid" passed
  // to user defined constructor
  ///////////////////////////////////////////
  // user defined constructor
  ///////////////////////////////////////////
  Lattice(GridBase *grid,ViewMode mode=AcceleratorWriteDiscard) { 
    this->_grid = grid;
    resize(this->_grid->oSites());
    assert((((uint64_t)&this->_odata[0])&0xF) ==0);
    this->checkerboard=0;
    SetViewMode(mode);
  }
  
  //  virtual ~Lattice(void) = default;
    
  void reset(GridBase* grid) {
    if (this->_grid != grid) {
      this->_grid = grid;
      this->resize(grid->oSites());
      this->checkerboard = 0;
    }
  }
  ///////////////////////////////////////////
  // copy constructor
  ///////////////////////////////////////////
  Lattice(const Lattice& r){ 
    this->_grid = r.Grid();
    resize(this->_grid->oSites());
    *this = r;
  }
  ///////////////////////////////////////////
  // move constructor
  ///////////////////////////////////////////
  Lattice(Lattice && r){ 
    this->_grid = r.Grid();
    this->_odata      = r._odata;
    this->_odata_size = r._odata_size;
    this->checkerboard= r.Checkerboard();
    r._odata      = nullptr;
    r._odata_size = 0;
  }
  ///////////////////////////////////////////
  // assignment template
  ///////////////////////////////////////////
  template<class robj> inline Lattice<vobj> & operator = (const Lattice<robj> & r){
    typename std::enable_if<!std::is_same<robj,vobj>::value,int>::type i=0;
    conformable(*this,r);
    this->checkerboard = r.Checkerboard();
    auto him= r.View(AcceleratorRead);
    auto me =   View(AcceleratorWriteDiscard);
    accelerator_for(ss,me.size(),vobj::Nsimd(),{
      coalescedWrite(me[ss],him(ss));
    });
    me.ViewClose();    him.ViewClose();
    return *this;
  }

  ///////////////////////////////////////////
  // Copy assignment 
  ///////////////////////////////////////////
  inline Lattice<vobj> & operator = (const Lattice<vobj> & r){
    this->checkerboard = r.Checkerboard();
    conformable(*this,r);
    auto him= r.View(AcceleratorRead);
    auto me =   View(AcceleratorWriteDiscard);
    accelerator_for(ss,me.size(),vobj::Nsimd(),{
      coalescedWrite(me[ss],him(ss));
    });
    me.ViewClose();    him.ViewClose();
    return *this;
  }
  ///////////////////////////////////////////
  // Move assignment possible if same type
  ///////////////////////////////////////////
  inline Lattice<vobj> & operator = (Lattice<vobj> && r){

    resize(0); // deletes if appropriate
    this->_grid       = r.Grid();
    this->_odata      = r._odata;
    this->_odata_size = r._odata_size;
    this->checkerboard= r.Checkerboard();

    r._odata      = nullptr;
    r._odata_size = 0;
    
    return *this;
  }

  /////////////////////////////////////////////////////////////////////////////
  // *=,+=,-= operators inherit behvour from correspond */+/- operation
  /////////////////////////////////////////////////////////////////////////////
  template<class T> inline Lattice<vobj> &operator *=(const T &r) {
    *this = (*this)*r;
    return *this;
  }
  
  template<class T> inline Lattice<vobj> &operator -=(const T &r) {
    *this = (*this)-r;
    return *this;
  }
  template<class T> inline Lattice<vobj> &operator +=(const T &r) {
    *this = (*this)+r;
    return *this;
  }

  friend inline void swap(Lattice &l, Lattice &r) { 
    conformable(l,r);
    LatticeAccelerator<vobj> tmp;
    LatticeAccelerator<vobj> *lp = (LatticeAccelerator<vobj> *)&l;
    LatticeAccelerator<vobj> *rp = (LatticeAccelerator<vobj> *)&r;
    tmp = *lp;    *lp=*rp;    *rp=tmp;
  }

}; // class Lattice

template<class vobj> std::ostream& operator<< (std::ostream& stream, const Lattice<vobj> &o){
  typedef typename vobj::scalar_object sobj;
  for(int64_t g=0;g<o.Grid()->_gsites;g++){

    Coordinate gcoor;
    o.Grid()->GlobalIndexToGlobalCoor(g,gcoor);

    sobj ss;
    peekSite(ss,o,gcoor);
    stream<<"[";
    for(int d=0;d<gcoor.size();d++){
      stream<<gcoor[d];
      if(d!=gcoor.size()-1) stream<<",";
    }
    stream<<"]\t";
    stream<<ss<<std::endl;
  }
  return stream;
}
  
NAMESPACE_END(Grid);

