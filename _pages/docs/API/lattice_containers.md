---
title : "API Documentation"
author_profile: false
excerpt: "Lattice containers"
header:
  overlay_color: "#5DADE2"
permalink: /docs/API/lattice_containers.html
sidebar:
  nav : docs
---

{% include toc icon="gears" title="Contents" %}

Lattice objects may be constructed to contain the local portion of a distribued array of any tensor type.
For performance reasons the tensor type uses a vector `Real` or `Complex` as the fundamental datum.

Every lattice requires a `GridBase` object pointer to be provided in its constructor. Memory is allocated
at construction time. If a Lattice is passed a RedBlack grid, it allocates
half the storage of the full grid, and may either store the red or black checkerboard. The Lattice object
will automatically track through assignments which checkerboard it refers to.
For example, shifting a Even checkerboard by an odd distance produces an Odd result field.

Struct of array objects are defined, and used in the template parameters to the lattice class.

**Example** (`lib/qcd/QCD.h`):

```c++
       template<typename vtype> using iSpinMatrix = iScalar<iMatrix<iScalar<vtype>, Ns> >;
       typedef iSpinMatrix<ComplexF>           SpinMatrixF; //scalar
       typedef iSpinMatrix<vComplexF>          vSpinMatrixF;//vectorised
       typedef Lattice<vSpinMatrixF>           LatticeSpinMatrixF;
```

The full range of QCD relevant lattice objects is given below.

|-----------|------------|----------|-----------|---------------|---------------------------------|--------------------------|
|   Lattice |    Lorentz |    Spin  |    Colour | scalar_type   |   Field                         |      Synonym             |
|-----------|------------|----------|-----------|---------------|---------------------------------|--------------------------|
|`Vector`   | `Scalar`   | `Scalar` |  `Scalar` |  `RealD`      | `LatticeRealD`                  |   N/A                    |
|`Vector`   | `Scalar`   | `Scalar` |  `Scalar` |  `ComplexD`   | `LatticeComplexD`               |   N/A                    |
|`Vector`   | `Scalar`   | `Scalar` |  `Matrix` |  `ComplexD`   | `LatticeColourMatrixD`          |   `LatticeGaugeLink`     |
|`Vector`   | `Vector`   | `Scalar` |  `Matrix` |  `ComplexD`   | `LatticeLorentzColourMatrixD`   |   `LatticeGaugeFieldD`   |
|`Vector`   | `Scalar`   | `Vector` |  `Vector` |  `ComplexD`   | `LatticeSpinColourVectorD`      |   `LatticeFermionD`      |
|`Vector`   | `Scalar`   | `Vector` |  `Vector` |  `ComplexD`   | `LatticeHalfSpinColourVectorD`  |   `LatticeHalfFermionD`  |
|`Vector`   | `Scalar`   | `Matrix` |  `Matrix` |  `ComplexD`   | `LatticeSpinColourMatrixD`      |   `LatticePropagatorD`   |
|-----------|------------|----------|-----------|---------------|---------------------------------|--------------------------|

Additional single precison variants are defined with the suffix `F`.
Other lattice objects can be defined using the sort of typedef's shown above if needed.

### Opaque containers

The layout within the container is complicated to enable maximum opportunity for vectorisation, and 
is opaque from the point of view of the API definition. The key implementation observation is that
so long as data parallel operations are performed and adjacent SIMD lanes correspond to well separated
lattice sites, then identical operations are performed on all SIMD lanes and enable good vectorisation.

Because the layout is opaque, import and export routines from naturally ordered x,y,z,t arrays
are provided (`lib/lattice/Lattice_transfer.h`):

```c++
    unvectorizeToLexOrdArray(std::vector<sobj> &out, const Lattice<vobj> &in);
    vectorizeFromLexOrdArray(std::vector<sobj> &in , Lattice<vobj> &out);
```

The Lexicographic order of data in the external vector fields is defined by (`lib/util/Lexicographic.h`):

```c++
    Lexicographic::IndexFromCoor(const Coordinate &lcoor, int &lex,Coordinate *local_dims);
```

This ordering is $$x + L_x * y + L_x*L_y*z + L_x*L_y*L_z *t$$

Peek and poke routines are provided to perform single site operations. These operations are
extremely low performance and are not intended for algorithm development or performance critical code.

The following are "collective" operations and involve communication between nodes. All nodes receive the same
result by broadcast from the owning node:

```c++
    void peekSite(sobj &s,const Lattice<vobj> &l,const Coordinate &site);
    void pokeSite(const sobj &s,Lattice<vobj> &l,const Coordinate &site);
```

The following are executed independently by each node:

```c++
    void peekLocalSite(sobj &s,const Lattice<vobj> &l,Coordinate &site);
    void pokeLocalSite(const sobj &s,Lattice<vobj> &l,Coordinate &site);
```

Lattices of one tensor type may be transformed into lattices of another tensor type by
peeking and poking specific indices in a data parallel manner:

```c++
    template<int Index,class vobj> // Vector data parallel index peek
    auto PeekIndex(const Lattice<vobj> &lhs,int i);

    template<int Index,class vobj> // Matrix data parallel index peek
    auto PeekIndex(const Lattice<vobj> &lhs,int i,int j);

    template<int Index,class vobj>   // Vector poke
    void PokeIndex(Lattice<vobj> &lhs,const Lattice<> & rhs,int i)

    template<int Index,class vobj>   // Matrix poke
    void PokeIndex(Lattice<vobj> &lhs,const Lattice<> & rhs,int i,int j)
```

The inconsistent capitalisation on the letter P is due to an obscure bug in g++ that has not to
our knowledge been fixed in any version. The bug was reported in 2016.

### Global Reduction operations

Reduction operations for any lattice field are provided. The result is identical on each computing node
that is part of the relevant Grid communicator:

```c++
  template<class vobj> 
  RealD norm2(const Lattice<vobj> &arg);

  template<class vobj> 
  ComplexD innerProduct(const Lattice<vobj> &left,const Lattice<vobj> &right);

  template<class vobj> 
  vobj sum(const Lattice<vobj> &arg)
```

### Site local reduction operations

Internal indices may be reduced, site by site, using the following routines:

```c++
  template<class vobj>
  auto localNorm2 (const Lattice<vobj> &rhs)

  template<class vobj>
  auto localInnerProduct (const Lattice<vobj> &lhs,const Lattice<vobj> &rhs) 
```

### Outer product

A site local outer product is defined:

```c++
  template<class ll,class rr>
  auto outerProduct (const Lattice<ll> &lhs,const Lattice<rr> &rhs)
```

### Slice operations

Slice operations are defined to operate on one lower dimension than the full lattice. The omitted dimension
is the parameter orthogdim:

```c++
  template<class vobj> 
  void sliceSum(const Lattice<vobj> &Data,
                std::vector<typename vobj::scalar_object> &result,
 	        int orthogdim);

  template<class vobj> 
  void sliceInnerProductVector( std::vector<ComplexD> & result, 
	                        const Lattice<vobj> &lhs,
	 			const Lattice<vobj> &rhs,
				int orthogdim); 

  template<class vobj> 
  void sliceNorm (std::vector<RealD> &sn,
       		  const Lattice<vobj> &rhs,
		  int orthogdim);
```

### Data parallel expression template engine

The standard arithmetic operators and some data parallel library functions are implemented site by site
on lattice types. 

Operations may only ever combine lattice objects that have been constructed from the **same** grid pointer.

**Example**:

```c++
    LatticeFermionD A(&grid);
    LatticeFermionD B(&grid);
    LatticeFermionD C(&grid);
    
    A = B - C;
```

Such operations are said to be **conformable** and are the lattice are guaranteed to have the same dimensions
and both MPI and SIMD decomposition because they are based on the same grid object. The conformability check
is lightweight and simply requires the same grid pointers be passed to the lattice objects. The data members
of the grid objects are not compared.

Conformable lattice fields may be combined with appropriate scalar types in expressions. The implemented
rules follow those already documented for the tensor types. 




### Unary operators and functions

The following sitewise unary operations are defined:

|-----------------------|---------------------------------------------|
| Operation             |    Description                              |
|-----------------------|---------------------------------------------|
|`operator-`            |  negate                                     |
|`adj`                  |  Hermitian conjugate                        |
|`conjugate`            |  complex conjugate                          |
|`trace`                |  sitewise trace                             |
|`transpose`            |  sitewise transpose                         |
|`Ta`                   |  take traceles anti Hermitian part          |
|`ProjectOnGroup`       |  reunitarise or orthogonalise               |
|`real`                 |  take the real part                         |
|`imag`                 |  take the imaginary part                    |
|`toReal`               |  demote complex to real                     |
|`toComplex`            |  promote real to complex                    |
|`timesI`               |  elementwise +i mult (0 is not multiplied)  |
|`timesMinusI`          |  elementwise -i mult (0 is not multiplied)  |
|`abs`                  |  elementwise absolute value                 |
|`sqrt`                 |     elementwise square root                 |
|`rsqrt`                |     elementwise reciprocal square root      |
|`sin`                  |     elementwise sine                        |
|`cos`                  |     elementwise cosine                      |
|`asin`                 |     elementwise inverse sine                |
|`acos`                 |     elementwise inverse cosine              |
|`log`                  |     elementwise logarithm                   |
|`exp`                  |     elementwise exponentiation              |
|`operator!`            |     Logical negation of integer field       |
|`Not`                  |     Logical negation of integer field       |
|-----------------------|---------------------------------------------|

The following sitewise applied functions with additional parameters are:

```c++
  template<class obj> Lattice<obj> pow(const Lattice<obj> &rhs_i,RealD y);

  template<class obj> Lattice<obj> mod(const Lattice<obj> &rhs_i,Integer y);

  template<class obj> Lattice<obj> div(const Lattice<obj> &rhs_i,Integer y);

  template<class obj> Lattice<obj> 
  expMat(const Lattice<obj> &rhs_i, RealD alpha, Integer Nexp = DEFAULT_MAT_EXP);
```

### Binary operators

The following binary operators are defined:

```
  operator+
  operator-
  operator*
  operator/
```

Logical are defined on LatticeInteger types:

```
  operator&
  operator|
  operator&&
  operator||
```

### Ternary operator, logical operatons and where

Within the data parallel level of the API the only way to perform operations
that are differentiated between sites is use predicated execution.

The predicate takes the form of a `LatticeInteger` which is confromable with both
the `iftrue` and `iffalse` argument:

```c++
  template<class vobj,class iobj> void where(const Lattice<iobj> &pred,
                                                   Lattice<vobj> &iftrue,
                                                   Lattice<vobj> &iffalse);
```
This plays the data parallel analogue of the C++ ternary operator:

```c++
     a = b ? c : d;
```

In order to create the predicate in a coordinate dependent fashion it is often useful
to use the lattice coordinates. 

The `LatticeCoordinate` function:

```c++
    template<class iobj> LatticeCoordinate(Lattice<iobj> &coor,int dir);
```

fills an `Integer` field with the coordinate in the N-th dimension.
A usage example is given

**Example**:

```c++
        int dir  =3;
        int block=4;
        LatticeInteger coor(FineGrid);

	LatticeCoordinate(coor,dir);
	
	result = where(mod(coor,block)==(block-1),x,z);
```

(Other usage cases of LatticeCoordinate include the generation of plane wave momentum phases.)

### Site local fused operations

The biggest limitation of expression template engines is that the optimisation 
visibility is a single assignment statement in the original source code.

There is no scope for loop fusion between multiple statements.
Multi-loop fusion gives scope for greater cache locality.

Two primitives for hardware aware parallel loops are provided.
These will operate directly on the site objects which are expanded by a factor
of the vector length (in our struct of array datatypes). 

Since the mapping of sites
to data lanes is opaque, these vectorised loops
are *only* appropriate for optimisation of site local operations.

### View objects

Due to an obscure aspect of the way that Nvidia handle device C++11 lambda functions,
it is necessary to disable the indexing of a Lattice object.

Rather, a reference to a lattice object must be first obtained.

The reference is copyable to a GPU, and is able to be indexed on either accelerator code,
or by host code.

In order to prevent people developing code that dereferences Lattice objects in a way that
works on CPU compilation, but fails on GPU compilation, we have decided to remove the ability
to index a lattice object on CPU code. 

As a result of Nvidia's constraints, all accesses to lattice objects are required to be made
through a View object. 

In the following, the type is `LatticeView<vobj>`, however it is wise to use the C++11 auto keyword
to avoid naming the type. See code examples below.


### thread_loops

The first parallel primitive is the thread_loop

**Example**:

```c++
  LatticeField r(grid); 
  LatticeField x(grid);
  LatticeField p(grid); 
  LatticeField mmp(grid);
  auto r_v = r.View();  
  auto x_v = x.View();
  auto p_v = p.View(); 
  auto mmp_v = mmp.View();
  thread_loop(s , r_v, {
    r_v[s] = r_v[s]   - a * mmp_v[s];
    x_v[s] = x_v[s]   + a*p_v[s];
    p_v[s] = p_v[s]*b + r_v[s];
  });
```

### accelerator_loops

The second parallel primitive is an accelerated_loop

**Example**:

```c++
  LatticeField r(grid); 
  LatticeField x(grid);
  LatticeField p(grid); 
  LatticeField mmp(grid);
  auto r_v = r.View();  
  auto x_v = x.View();
  auto p_v = p.View(); 
  auto mmp_v = mmp.View();
  accelerator_loop(s , r_v, {
    r_v[s] = r_v[s]   - a * mmp_v[s];
    x_v[s] = x_v[s]   + a*p_v[s];
    p_v[s] = p_v[s]*b + r_v[s];
  });
```


### Cshift 

Site shifting operations are provided using the Cshift function:

```c++
  template<class vobj> 
  Lattice<vobj> Cshift(const Lattice<vobj> &rhs,int dimension,int shift)
```

This shifts the whole vector by any distance shift in the appropriate dimension.

For the avoidance of doubt on direction conventions,a positive shift moves the 
lattice site $$x_mu = 1$$ in the rhs to $$x_mu = 0$$ in the result.

**Example** (`benchmarks/Benchmark_wilson.cc`):

```c++
  { // Naive wilson implementation
    ref = Zero();
    for(int mu=0;mu<Nd;mu++){

      tmp = U[mu]*Cshift(src,mu,1);
      {
	auto ref_v = ref.View();
	auto tmp_v = tmp.View();
	for(int i=0;i<ref_v.size();i++){
	  ref_v[i]+= tmp_v[i] - Gamma(Gmu[mu])*tmp_v[i]; ;
	}
      }

      tmp =adj(U[mu])*src;
      tmp =Cshift(tmp,mu,-1);
      {
	auto ref_v = ref.View();
	auto tmp_v = tmp.View();
	for(int i=0;i<ref_v.size();i++){
	  ref_v[i]+= tmp_v[i] + Gamma(Gmu[mu])*tmp_v[i]; ;
	}
      }
    }
  }
```

### CovariantCshift 

Covariant Cshift operations are provided for common cases of boundary condition. These may be further optimised
in future:

```c++
  template<class covariant,class gauge> 
  Lattice<covariant> CovShiftForward(const Lattice<gauge> &Link, int mu,
			   	     const Lattice<covariant> &field);

  template<class covariant,class gauge> 
  Lattice<covariant> CovShiftBackward(const Lattice<gauge> &Link, int mu,
			              const Lattice<covariant> &field);
```

### Boundary conditions

The covariant shift routines occur in namespaces PeriodicBC and ConjugateBC. The correct covariant shift
for the boundary condition is passed into the gauge actions and wilson loops via an
"Impl" template policy class.

The relevant staples, plaquettes, and loops are formed by using the provided method:

```c++
    Impl::CovShiftForward
    Impl::CovShiftBackward
```

etc... This makes physics code transform appropriately with externally supplied rules about
treating the boundary.

**Example** (`lib/qcd/util/WilsonLoops.h`):

```c++
  static void dirPlaquette(GaugeMat &plaq, const std::vector<GaugeMat> &U,
                           const int mu, const int nu) {
    // ___
    //|   |
    //|<__|
    plaq = Gimpl::CovShiftForward(U[mu],mu,
                    Gimpl::CovShiftForward(U[nu],nu,
			Gimpl::CovShiftBackward(U[mu],mu,
			   Gimpl::CovShiftIdentityBackward(U[nu], nu))));
  }
```

### Inter-grid transfer operations

Transferring between different checkerboards of the same global lattice:

```c++
  template<class vobj> void pickCheckerboard(int cb,Lattice<vobj> &half,const Lattice<vobj> &full);
  template<class vobj> void setCheckerboard(Lattice<vobj> &full,const Lattice<vobj> &half);
```

These are used to set up Schur red-black decomposed solvers, for example.

Multi-grid projection between a fine and coarse grid:

```c++
 template<class vobj,class CComplex,int nbasis>
 void blockProject(Lattice<iVector<CComplex,nbasis > > &coarseData,
                   const             Lattice<vobj>   &fineData,
                   const std::vector<Lattice<vobj> > &Basis);
```

Multi-grid promotion to a finer grid:

```c++
  template<class vobj,class CComplex,int nbasis>
  void blockPromote(const Lattice<iVector<CComplex,nbasis > > &coarseData,
                    Lattice<vobj>   &fineData,
                    const std::vector<Lattice<vobj> > &Basis)
```

Support for sub-block Linear algebra:

```c++
  template<class vobj,class CComplex>
  void blockZAXPY(Lattice<vobj> &fineZ,
                  const Lattice<CComplex> &coarseA,
                  const Lattice<vobj> &fineX,
                  const Lattice<vobj> &fineY)

  template<class vobj,class CComplex>
  void blockInnerProduct(Lattice<CComplex> &CoarseInner,
                         const Lattice<vobj> &fineX,
                         const Lattice<vobj> &fineY)

  template<class vobj,class CComplex>
  void blockNormalise(Lattice<CComplex> &ip,Lattice<vobj> &fineX)

  template<class vobj>
  void blockSum(Lattice<vobj> &coarseData,const Lattice<vobj> &fineData)

  template<class vobj,class CComplex>
  void blockOrthogonalise(Lattice<CComplex> &ip,std::vector<Lattice<vobj> > &Basis)
```

Conversion between different SIMD layouts:

```c++
  template<class vobj,class vvobj>
  void localConvert(const Lattice<vobj> &in,Lattice<vvobj> &out)
```

Slices between grid of dimension N and grid of dimentions N+1:

```c++
  template<class vobj>
  void InsertSlice(const Lattice<vobj> &lowDim,Lattice<vobj> & higherDim,int slice, int orthog)

  template<class vobj>
  void ExtractSlice(Lattice<vobj> &lowDim,const Lattice<vobj> & higherDim,int slice, int orthog)
```

Growing a lattice by a multiple factor, with periodic replication:

```c++
  template<class vobj>
  void Replicate(Lattice<vobj> &coarse,Lattice<vobj> & fine)
```

That latter is useful to, for example, pre-thermalise a smaller volume and then grow the volume in HMC.
It was written while debugging G-parity boundary conditions.
