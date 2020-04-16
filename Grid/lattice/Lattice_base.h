/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/lattice/Lattice_base.h

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

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
			   /*  END LEGAL */
#pragma once 

#define STREAMING_STORES

NAMESPACE_BEGIN(Grid);

extern int GridCshiftPermuteMap[4][16];

///////////////////////////////////////////////////////////////////
// Base class which can be used by traits to pick up behaviour
///////////////////////////////////////////////////////////////////
class LatticeBase {};

/////////////////////////////////////////////////////////////////////////////////////////
// Conformable checks; same instance of Grid required
/////////////////////////////////////////////////////////////////////////////////////////
void accelerator_inline conformable(GridBase *lhs,GridBase *rhs)
{
  assert(lhs == rhs);
}

////////////////////////////////////////////////////////////////////////////
// Minimal base class containing only data valid to access from accelerator
// _odata will be a managed pointer in CUDA
////////////////////////////////////////////////////////////////////////////
// Force access to lattice through a view object.
// prevents writing of code that will not offload to GPU, but perhaps annoyingly
// strict since host could could in principle direct access through the lattice object
// Need to decide programming model.
#define LATTICE_VIEW_STRICT
template<class vobj> class LatticeAccelerator : public LatticeBase
{
protected:
  GridBase *_grid;
  int checkerboard;
  vobj     *_odata;    // A managed pointer
  uint64_t _odata_size;    
public:
  accelerator_inline LatticeAccelerator() : checkerboard(0), _odata(nullptr), _odata_size(0), _grid(nullptr) { }; 
  accelerator_inline uint64_t oSites(void) const { return _odata_size; };
  accelerator_inline int  Checkerboard(void) const { return checkerboard; };
  accelerator_inline int &Checkerboard(void) { return this->checkerboard; }; // can assign checkerboard on a container, not a view
  accelerator_inline void Conformable(GridBase * &grid) const
  { 
    if (grid) conformable(grid, _grid);
    else      grid = _grid;
  };
};

/////////////////////////////////////////////////////////////////////////////////////////
// A View class which provides accessor to the data.
// This will be safe to call from accelerator_for and is trivially copy constructible
// The copy constructor for this will need to be used by device lambda functions
/////////////////////////////////////////////////////////////////////////////////////////
template<class vobj> 
class LatticeView : public LatticeAccelerator<vobj>
{
public:


  // Rvalue
#ifdef __CUDA_ARCH__
  accelerator_inline const typename vobj::scalar_object operator()(size_t i) const { return coalescedRead(this->_odata[i]); }
#else 
  accelerator_inline const vobj & operator()(size_t i) const { return this->_odata[i]; }
#endif

  accelerator_inline const vobj & operator[](size_t i) const { return this->_odata[i]; };
  accelerator_inline vobj       & operator[](size_t i)       { return this->_odata[i]; };

  accelerator_inline uint64_t begin(void) const { return 0;};
  accelerator_inline uint64_t end(void)   const { return this->_odata_size; };
  accelerator_inline uint64_t size(void)  const { return this->_odata_size; };

  LatticeView(const LatticeAccelerator<vobj> &refer_to_me) : LatticeAccelerator<vobj> (refer_to_me)
  {
  }
};

/////////////////////////////////////////////////////////////////////////////////////////
// Lattice expression types used by ET to assemble the AST
// 
// Need to be able to detect code paths according to the whether a lattice object or not
// so introduce some trait type things
/////////////////////////////////////////////////////////////////////////////////////////

class LatticeExpressionBase {};

template <typename T> using is_lattice = std::is_base_of<LatticeBase, T>;
template <typename T> using is_lattice_expr = std::is_base_of<LatticeExpressionBase,T >;

template<class T, bool isLattice> struct ViewMapBase { typedef T Type; };
template<class T>                 struct ViewMapBase<T,true> { typedef LatticeView<typename T::vector_object> Type; };
template<class T> using ViewMap = ViewMapBase<T,std::is_base_of<LatticeBase, T>::value >;

template <typename Op, typename _T1>                           
class LatticeUnaryExpression : public  LatticeExpressionBase 
{
public:
  typedef typename ViewMap<_T1>::Type T1;
  Op op;
  T1 arg1;
  LatticeUnaryExpression(Op _op,const _T1 &_arg1) : op(_op), arg1(_arg1) {};
};

template <typename Op, typename _T1, typename _T2>              
class LatticeBinaryExpression : public LatticeExpressionBase 
{
public:
  typedef typename ViewMap<_T1>::Type T1;
  typedef typename ViewMap<_T2>::Type T2;
  Op op;
  T1 arg1;
  T2 arg2;
  LatticeBinaryExpression(Op _op,const _T1 &_arg1,const _T2 &_arg2) : op(_op), arg1(_arg1), arg2(_arg2) {};
};

template <typename Op, typename _T1, typename _T2, typename _T3> 
class LatticeTrinaryExpression : public LatticeExpressionBase 
{
public:
  typedef typename ViewMap<_T1>::Type T1;
  typedef typename ViewMap<_T2>::Type T2;
  typedef typename ViewMap<_T3>::Type T3;
  Op op;
  T1 arg1;
  T2 arg2;
  T3 arg3;
  LatticeTrinaryExpression(Op _op,const _T1 &_arg1,const _T2 &_arg2,const _T3 &_arg3) : op(_op), arg1(_arg1), arg2(_arg2), arg3(_arg3) {};
};

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
  // Return a view object that may be dereferenced in site loops.
  // The view is trivially copy constructible and may be copied to an accelerator device
  // in device lambdas
  /////////////////////////////////////////////////////////////////////////////////
  LatticeView<vobj> View (void) const 
  {
    LatticeView<vobj> accessor(*( (LatticeAccelerator<vobj> *) this));
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
    GridBase *egrid(nullptr);
    GridFromExpression(egrid,expr);
    assert(egrid!=nullptr);
    conformable(this->_grid,egrid);

    int cb=-1;
    CBFromExpression(cb,expr);
    assert( (cb==Odd) || (cb==Even));
    this->checkerboard=cb;

    auto me  = View();
    accelerator_for(ss,me.size(),1,{
      auto tmp = eval(ss,expr);
      vstream(me[ss],tmp);
    });
    return *this;
  }
  template <typename Op, typename T1,typename T2> inline Lattice<vobj> & operator=(const LatticeBinaryExpression<Op,T1,T2> &expr)
  {
    GridBase *egrid(nullptr);
    GridFromExpression(egrid,expr);
    assert(egrid!=nullptr);
    conformable(this->_grid,egrid);

    int cb=-1;
    CBFromExpression(cb,expr);
    assert( (cb==Odd) || (cb==Even));
    this->checkerboard=cb;

    auto me  = View();
    accelerator_for(ss,me.size(),1,{
      auto tmp = eval(ss,expr);
      vstream(me[ss],tmp);
    });
    return *this;
  }
  template <typename Op, typename T1,typename T2,typename T3> inline Lattice<vobj> & operator=(const LatticeTrinaryExpression<Op,T1,T2,T3> &expr)
  {
    GridBase *egrid(nullptr);
    GridFromExpression(egrid,expr);
    assert(egrid!=nullptr);
    conformable(this->_grid,egrid);

    int cb=-1;
    CBFromExpression(cb,expr);
    assert( (cb==Odd) || (cb==Even));
    this->checkerboard=cb;
    auto me  = View();
    accelerator_for(ss,me.size(),1,{
      auto tmp = eval(ss,expr);
      vstream(me[ss],tmp);
    });
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
    auto me  = View();
    thread_for(ss,me.size(),{
      me[ss] = r;
    });
    return *this;
  }

  //////////////////////////////////////////////////////////////////
  // Follow rule of five, with Constructor requires "grid" passed
  // to user defined constructor
  ///////////////////////////////////////////
  // user defined constructor
  ///////////////////////////////////////////
  Lattice(GridBase *grid) { 
    this->_grid = grid;
    resize(this->_grid->oSites());
    assert((((uint64_t)&this->_odata[0])&0xF) ==0);
    this->checkerboard=0;
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
    //    std::cout << "Lattice constructor(const Lattice &) "<<this<<std::endl; 
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
    auto me =   View();
    auto him= r.View();
    accelerator_for(ss,me.size(),vobj::Nsimd(),{
      coalescedWrite(me[ss],him(ss));
    });
    return *this;
  }

  ///////////////////////////////////////////
  // Copy assignment 
  ///////////////////////////////////////////
  inline Lattice<vobj> & operator = (const Lattice<vobj> & r){
    this->checkerboard = r.Checkerboard();
    conformable(*this,r);
    auto me =   View();
    auto him= r.View();
    accelerator_for(ss,me.size(),vobj::Nsimd(),{
      coalescedWrite(me[ss],him(ss));
    });
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
  for(int g=0;g<o.Grid()->_gsites;g++){

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

