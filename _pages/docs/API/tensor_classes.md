---
title : "API Documentation"
author_profile: false
excerpt: "Tensor classes"
header:
  overlay_color: "#5DADE2"
permalink: /docs/API/tensor_classes.html
sidebar:
  nav : docs
---

The Tensor data structures are built up from fundamental
scalar matrix and vector classes:

```c++
    template<class vobj      > class iScalar { private: vobj _internal ; }
    template<class vobj,int N> class iVector { private: vobj _internal[N] ; }
    template<class vobj,int N> class iMatrix { private: vobj _internal[N] ; }
```

These are template classes and can be passed a fundamental scalar or vector type, or
nested to form arbitrarily complicated tensor products of indices. All mathematical expressions
are defined to operate recursively, index by index.

Presently the constants

* `Nc`
* `Nd`

are globally predefined. However, this is planned for changed in future and policy classes
for different theories (e.g. QCD, QED, SU2 etc...) will contain these constants and enable multiple
theories to coexist more naturally.

Arbitrary tensor products of fundamental scalar, vector
and matrix objects may be formed in principle by the basic Grid code.

For Lattice field theory, we define types according to the following tensor
product structure ordering. The suffix "D" indicates either double types, and
replacing with "F" gives the corresponding single precision type.

|Lattice  | Lorentz  |  Spin     |   Colour | scalar_type |  Field                   |
|---------|----------|-----------|----------|-------------|--------------------------|
|Scalar   | Scalar   |  Scalar   | Scalar   |  RealD      |   RealD                  |
|Scalar   | Scalar   |  Scalar   | Scalar   |  ComplexD   |   ComplexD               |
|Scalar   | Scalar   |  Scalar   | Matrix   |  ComplexD   |   ColourMatrixD          |
|Scalar   | Vector   |  Scalar   | Matrix   |  ComplexD   |   LorentzColourMatrixD   |
|Scalar   | Scalar   |  Vector   | Vector   |  ComplexD   |   SpinColourVectorD      |
|Scalar   | Scalar   |  Vector   | Vector   |  ComplexD   |   HalfSpinColourVectorD  |
|Scalar   | Scalar   |  Matrix   | Matrix   |  ComplexD   |   SpinColourMatrixD      |
|---------|----------|-----------|----------|-------------|--------------------------|

The types are implemented via a recursive tensor nesting system.

**Example** we declare:

```c++
  template<typename vtype>
  using iLorentzColourMatrix = iVector<iScalar<iMatrix<vtype, Nc> >, Nd > ;

  typedef iLorentzColourMatrix<ComplexD > LorentzColourMatrixD;
```

**Example** we declare:

```c++
  template<typename vtype>
  using iLorentzColourMatrix = iVector<iScalar<iMatrix<vtype, Nc> >, Nd > ;

  typedef iLorentzColourMatrix<ComplexD > LorentzColourMatrixD;
```

Arbitrarily deep tensor nests may be formed. Grid uses a positional and numerical rule to associate indices for contraction
in the Einstein summation sense.

|----------------|---------|------------|
| Symbolic name  | Number  | Position   |
|----------------|---------|------------|
|`LorentzIndex`  | 0       | left       |
|`SpinIndex`     | 1       | middle     |
|`ColourIndex`   | 2       | right      |
|----------------|---------|------------|

The conventions are that the index ordering left to right are: Lorentz, Spin, Colour. A scalar type (either real
or complex, single or double precision) is be provided to the innermost structure.


### Tensor arithmetic rules (`lib/tensors/Tensor_arith.h`)

Arithmetic rules are defined on these types

The multiplication operator follows the natural multiplication
table for each index, index level by index level.

`Operator *`

|----|----|----|----|
| x  |  S |  V |  M |
|----|----|----|----|
| S  | S  | V  | M  |
| V  | S  | S  | V  |
| M  | M  | V  | M  |
|----|----|----|----|

The addition and subtraction rules disallow a scalar to be added to a vector,
and vector to be added to matrix. A scalar adds to a matrix on the diagonal.

`Operator +` and `Operator -`

|----|----|----|----|
| +/-|  S |  V |  M |
|----|----|----|----|
| S  | S  | -  | M  |
| V  | -  | V  | -  |
| M  | M  | -  | M  |
|----|----|----|----|


The rules for a nested objects are recursively inferred level by level from basic rules of multiplication
addition and subtraction for scalar/vector/matrix. Legal expressions can only be formed between objects
with the same number of nested internal indices. All the Grid QCD datatypes have precisely three internal
indices, some of which may be trivial scalar to enable expressions to be formed.

Arithmetic operations are possible where the left or right operand is a scalar type.

**Example**:

```c++
    LatticeColourMatrixD U(grid);
    LatticeColourMatrixD Udag(grid);

    Udag = adj(U);

    RealD unitary_err = norm2(U*adj(U) - 1.0);
```

Will provide a measure of how discrepant from unitarity the matrix U is.

### Internal index manipulation (`lib/tensors/Tensor_index.h`)

General code can access any specific index by number with a peek/poke semantic:

```c++
   // peek index number "Level" of a vector index
   template<int Level,class vtype>  auto peekIndex (const vtype &arg,int i);

   // peek index number "Level" of a vector index
   template<int Level,class vtype>  auto peekIndex (const vtype &arg,int i,int j);

   // poke index number "Level" of a vector index
   template<int Level,class vtype>
   void pokeIndex (vtype &pokeme,arg,int i)

   // poke index number "Level" of a matrix index
   template<int Level,class vtype>
   void pokeIndex (vtype &pokeme,arg,int i,int j)
```

**Example**:

```c++
    for (int mu = 0; mu < Nd; mu++) {
      U[mu] = PeekIndex<LorentzIndex>(Umu, mu);
    }
```

Similar to the QDP++ package convenience routines are provided to access specific elements of
vector and matrix internal index types by physics name or meaning aliases for the above routines
with the appropriate index constant.

* `peekColour`
* `peekSpin`
* `peekLorentz`

and

* `pokeColour`
* `pokeSpin`
* `pokeLorentz`

For example, we often store Gauge Fields with a Lorentz index, but can split them into
polarisations in relevant pieces of code.

**Example**:

```c++
    for (int mu = 0; mu < Nd; mu++) {
      U[mu] = peekLorentz(Umu, mu);
    }
```

For convenience, direct access as both an l-value and an r-value is provided by the parenthesis operator () on each of the Scalar, Vector and Matrix classes.
For example one may write

**Example**:

```c++
  ColourMatrix A, B;
  A()()(i,j) = B()()(j,i);
```

bearing in mind that empty parentheses are need to address a scalar entry in the tensor index nest.

The first (left) empty parentheses move past the (scalar) Lorentz level in the tensor nest, and the second
(middle) empty parantheses move past the (scalar) spin level. The (i,j) index the colour matrix.

Other examples are easy to form for the many cases, and should be obvious to the reader.
This form of addressing is convenient and saves peek, modifying, poke
multiple temporary objects when both spin and colour indices are being accessed.
There are many cases where multiple lines of code are required with a peek/poke semantic which are
easier with direct l-value and r-value addressing.


### Matrix operations

Transposition and tracing specific internal indices are possible using:

```c++
  template<int Level,class vtype>  
  auto traceIndex (const vtype &arg)

  template<int Level,class vtype>  
  auto transposeIndex (const vtype &arg)
```

These may be used as

**Example**:

```c++
  LatticeColourMatrixD Link(grid);
  ComplexD link_trace = traceIndex<ColourIndex> (Link);
```

Again, convenience aliases for QCD naming schemes are provided via       

* `traceColour`
* `traceSpin`

* `transposeColour`
* `transposeSpin`

**Example**:

```c++
  ComplexD link_trace = traceColour (Link);
```

The operations only makes sense for matrix and scalar internal indices.

The trace and transpose over all indices is also defined for matrix and scalar types:

```c++
   template<class vtype,int N> 
   auto trace(const iMatrix<vtype,N> &arg) -> iScalar

   template<class vtype,int N> 
   auto transpose(const iMatrix<vtype,N> &arg  ) -> iMatrix
```

Similar functions are:

* `conjugate`
* `adjoint`

The traceless anti-Hermitian part is taken with:

```c++
    template<class vtype,int N> iMatrix<vtype,N> 
    Ta(const iMatrix<vtype,N> &arg)
```

Reunitarisation (or reorthogonalisation) is enabled by:

```c++
    template<class vtype,int N> iMatrix<vtype,N> 
    ProjectOnGroup(const iMatrix<vtype,N> &arg)
```

**Example**:

```c++
  LatticeColourMatrixD Mom(grid);
  LatticeColourMatrixD TaMom(grid);
  TaMom  = Ta(Mom);
```

### Querying internal index structure

Templated code may find it useful to use query functions on the Grid datatypes they are provided.
For example general Serialisation and I/O code can inspect the nature of a type a routine has been
asked to read from disk, or even generate descriptive type strings:

```c++
      ////////////////////////////////////////////////////
      // Support type queries on template params:
      ////////////////////////////////////////////////////
      // int _ColourScalar  =  isScalar<ColourIndex,vobj>();
      // int _ColourVector  =  isVector<ColourIndex,vobj>();
      // int _ColourMatrix  =  isMatrix<ColourIndex,vobj>();
      template<int Level,class vtype>  int isScalar(void)
      template<int Level,class vtype>  int isVector(void)
      template<int Level,class vtype>  int isMatrix(void)
```

**Example** (`lib/parallelIO/IldgIO.h`):

```c++
  template<class vobj> std::string ScidacRecordTypeString(int &colors, int &spins, int & typesize,int &datacount) { 

  /////////////////////////////////////////
  // Encode a generic tensor as a string
  /////////////////////////////////////////

  typedef typename getPrecision<vobj>::real_scalar_type stype;

  int _ColourN       = indexRank<ColourIndex,vobj>();
  int _ColourScalar  =  isScalar<ColourIndex,vobj>();
  int _ColourVector  =  isVector<ColourIndex,vobj>();
  int _ColourMatrix  =  isMatrix<ColourIndex,vobj>();

  int _SpinN       = indexRank<SpinIndex,vobj>();
  int _SpinScalar  =  isScalar<SpinIndex,vobj>();
  int _SpinVector  =  isVector<SpinIndex,vobj>();
  int _SpinMatrix  =  isMatrix<SpinIndex,vobj>();

  int _LorentzN       = indexRank<LorentzIndex,vobj>();
  int _LorentzScalar  =  isScalar<LorentzIndex,vobj>();
  int _LorentzVector  =  isVector<LorentzIndex,vobj>();
  int _LorentzMatrix  =  isMatrix<LorentzIndex,vobj>();

  std::stringstream stream;

  stream << "GRID_";
  stream << ScidacWordMnemonic<stype>();

  if ( _LorentzVector )   stream << "_LorentzVector"<<_LorentzN;
  if ( _LorentzMatrix )   stream << "_LorentzMatrix"<<_LorentzN;

  if ( _SpinVector )   stream << "_SpinVector"<<_SpinN;
  if ( _SpinMatrix )   stream << "_SpinMatrix"<<_SpinN;

  if ( _ColourVector )   stream << "_ColourVector"<<_ColourN;
  if ( _ColourMatrix )   stream << "_ColourMatrix"<<_ColourN;

  if ( _ColourScalar && _LorentzScalar && _SpinScalar )   stream << "_Complex";

  typesize = sizeof(typename vobj::scalar_type);

  if ( _ColourMatrix ) typesize*= _ColourN*_ColourN;
  else                 typesize*= _ColourN;

  if ( _SpinMatrix )   typesize*= _SpinN*_SpinN;
  else                 typesize*= _SpinN;

  };
```

### Inner and outer products

We recursively define (`tensors/Tensor_inner.h`), ultimately returning scalar in all indices:

```c++
  /////////////////////////////////////////////////////////////////////////
  // innerProduct Scalar x Scalar -> Scalar
  // innerProduct Vector x Vector -> Scalar
  // innerProduct Matrix x Matrix -> Scalar
  /////////////////////////////////////////////////////////////////////////
  template<class l,class r>       
  auto innerProductD (const iScalar<l>& lhs,const iScalar<r>& rhs)

  template<class l,class r,int N> 
  auto innerProductD (const iVector<l,N>& lhs,const iVector<r,N>& rhs)

  template<class l,class r,int N> 
  auto innerProductD (const iMatrix<l,N>& lhs,const iMatrix<r,N>& rhs)

  template<class l,class r>       
  auto innerProduct (const iScalar<l>& lhs,const iScalar<r>& rhs)

  template<class l,class r,int N> 
  auto innerProduct (const iVector<l,N>& lhs,const iVector<r,N>& rhs)

  template<class l,class r,int N> 
  auto innerProduct (const iMatrix<l,N>& lhs,const iMatrix<r,N>& rhs)
```

The sum is always performed in double precision for the `innerProductD` variant.

We recursively define (`tensors/Tensor_outer.h`):

```c++
  /////////////////////////////////////////////////////////////////////////
  // outerProduct Scalar x Scalar -> Scalar
  //              Vector x Vector -> Matrix
  /////////////////////////////////////////////////////////////////////////
  template<class l,class r> 
  auto outerProduct (const iScalar<l>& lhs,const iScalar<r>& rhs)

  template<class l,class r,int N> 
  auto outerProduct (const iVector<l,N>& lhs,const iVector<r,N>& rhs)
```

### Functions of Tensor

The following unary functions are defined, which operate element by element on a tensor 
data structure:

```c++
  sqrt();
  rsqrt();
  sin();
  cos();
  asin();
  acos();
  log();
  exp();
  abs();
  Not();
  toReal();
  toComplex();
```

Element wise functions are defined for::

```c++
  div(tensor,Integer);
  mod(tensor,Integer);
  pow(tensor,RealD);
```

Matrix exponentiation (as opposed to element wise exponentiation is implemented via power series in::

```c++
    Exponentiate(const Tensor &r  ,RealD alpha, Integer Nexp = DEFAULT_MAT_EXP)
```

the exponentiation is distributive across vector indices (i.e. proceeds component by component for a `LorentzColourMatrix`).

Determinant is similar::

```c++
    iScalar Determinant(const Tensor &r )
```