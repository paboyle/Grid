#pragma once
NAMESPACE_BEGIN(Grid);
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
class LatticeExprView : public LatticeAccelerator<vobj>
{
public:
  // Rvalue
  ViewMode mode;
  void * cpu_ptr;
#ifdef GRID_SIMT
  accelerator_inline const typename vobj::scalar_object operator()(size_t i) const { 
    return coalescedRead(this->_odata[i]); 
  }
#else 
  accelerator_inline const vobj & operator()(size_t i) const { return this->_odata[i]; }
#endif

  accelerator_inline const vobj & operator[](size_t i) const { return this->_odata[i]; };
  accelerator_inline vobj       & operator[](size_t i)       { return this->_odata[i]; };

  accelerator_inline uint64_t begin(void) const { return 0;};
  accelerator_inline uint64_t end(void)   const { return this->_odata_size; };
  accelerator_inline uint64_t size(void)  const { return this->_odata_size; };

  LatticeExprView(const LatticeAccelerator<vobj> &refer_to_me) : LatticeAccelerator<vobj> (refer_to_me){}

  // Host functions
  void ViewOpen(ViewMode mode)
  { // Translate the pointer, could save a copy. Could use a "Handle" and not save _odata originally in base
    //    std::cout << "View Open"<<std::hex<<this->_odata<<std::dec <<std::endl;
    this->cpu_ptr = (void *)this->_odata;
    this->mode    = mode;
    this->_odata  =(vobj *)
      AllocationCache::ViewOpen(this->cpu_ptr,
				this->_odata_size*sizeof(vobj),
				mode,
				AdviseDefault);    
  }
  void ViewClose(void)
  { // Inform the manager
    //    std::cout << "View Close"<<std::hex<<this->cpu_ptr<<std::dec <<std::endl;
    AllocationCache::ViewClose(this->cpu_ptr,this->mode);    
  }

};


///////////////////////////////////////////////////////////////////////
// An object to be stored in a shared_ptr to clean up after last view.
// UserView constructor,destructor updates view manager
// Non-copyable object??? Second base with copy/= deleted?
///////////////////////////////////////////////////////////////////////
class MemViewDeleter {
 public:
  void *cpu_ptr;
  ViewMode mode;
  ~MemViewDeleter(){
    AllocationCache::ViewClose(cpu_ptr,mode);    
  }
};
template<class vobj> 
class LatticeView : public LatticeExprView<vobj>
{
  std::shared_ptr<MemViewDeleter> Deleter;
public:
  LatticeView(const LatticeView<vobj> &orig) : LatticeExprView<vobj>(orig) {  }
  LatticeView(const LatticeAccelerator<vobj> &refer_to_me,ViewMode mode) : 
    LatticeExprView<vobj> (refer_to_me), Deleter(new MemViewDeleter)
  {
    //    std::cout << "FIXME - copy shared pointer? View Open in LatticeView"<<std::hex<<this->_odata<<std::dec<<" mode "<<mode <<std::endl;
    this->ViewOpen(mode);
    Deleter->cpu_ptr = this->cpu_ptr;
    Deleter->mode    = mode;
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
template<class T>                 struct ViewMapBase<T,true> { typedef LatticeExprView<typename T::vector_object> Type; };
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
NAMESPACE_END(Grid);
