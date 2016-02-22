    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_multishift_sqrt.cc

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
#include <Grid.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

template<class Field> class DumbOperator  : public LinearOperatorBase<Field> {
public:
  LatticeComplex scale;
  LatticeComplex sqrtscale;
  DumbOperator(GridBase *grid)
    : scale(grid),
      sqrtscale(grid)
  {
    GridParallelRNG  pRNG(grid);  
    std::vector<int> seeds({5,6,7,8});
    pRNG.SeedFixedIntegers(seeds);

    random(pRNG,sqrtscale);
    sqrtscale = 0.5*(sqrtscale + conjugate(sqrtscale));
    sqrtscale = sqrtscale*3.0+0.5;// force real pos def
    scale = sqrtscale *sqrtscale; //scale should be bounded by 12.25

    // 
    //    sqrtscale = 2.0;
    //    scale = 4.0;
  }
  // Support for coarsening to a multigrid
  void OpDiag (const Field &in, Field &out) {};
  void OpDir  (const Field &in, Field &out,int dir,int disp){};

  void Op     (const Field &in, Field &out){
    out = scale * in;
  }
  void AdjOp  (const Field &in, Field &out){
    out = scale * in;
  }
  void HermOp(const Field &in, Field &out){
    double n1, n2;
    HermOpAndNorm(in,out,n1,n2);
  }
  void HermOpAndNorm(const Field &in, Field &out,double &n1,double &n2){
    ComplexD dot;

    out = scale * in;

    dot= innerProduct(in,out);
    n1=real(dot);

    dot = innerProduct(out,out);
    n2=real(dot);
  }
  void ApplySqrt(const Field &in, Field &out){
    out = sqrtscale * in;
  }
  void ApplyInverse(const Field &in, Field &out){
    out = pow(scale,-1.0) * in;
  }
};

RealD InverseApproximation(RealD x){
  return 1.0/x;
}
RealD SqrtApproximation(RealD x){
  return std::sqrt(x);
}


int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  GridCartesian *grid = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), 
						       GridDefaultSimd(Nd,vComplex::Nsimd()),
						       GridDefaultMpi());

  double     lo=0.001;
  double     hi=1.0;
  int precision=64;
  int    degree=10;
  AlgRemez remez(lo,hi,precision);

  ////////////////////////////////////////
  // sqrt and inverse sqrt
  ////////////////////////////////////////
  std::cout<<GridLogMessage << "Generating degree "<<degree<<" for x^(1/2)"<<std::endl;
  remez.generateApprox(degree,1,2);

  MultiShiftFunction Sqrt(remez,1.0e-6,false);

  GridParallelRNG  pRNG(grid);  
  std::vector<int> seeds({1,2,3,4});
  pRNG.SeedFixedIntegers(seeds);

  LatticeFermion    src(grid); random(pRNG,src);
  LatticeFermion    combined(grid);
  LatticeFermion    reference(grid);
  LatticeFermion    error(grid);
  std::vector<LatticeFermion> result(degree,grid);
  LatticeFermion    summed(grid);

  ConjugateGradientMultiShift<LatticeFermion> MSCG(10000,Sqrt);

  DumbOperator<LatticeFermion> Diagonal(grid);

  MSCG(Diagonal,src,result);


  //  double a = norm;
  //  for(int n=0;n<poles.size();n++){
  //    a = a + residues[n]/(x+poles[n]);
  //  }
  assert(Sqrt.order==degree);

  combined = Sqrt.norm*src;
  for(int i=0;i<degree;i++){
    combined = combined + Sqrt.residues[i]*result[i];
  }
  
  Diagonal.ApplySqrt(src,reference);

  error = reference - combined;

  std::cout<<GridLogMessage << " Reference "<<norm2(reference)<<std::endl;
  std::cout<<GridLogMessage << " combined  "<<norm2(combined) <<std::endl;
  std::cout<<GridLogMessage << " error     "<<norm2(error)    <<std::endl;

  MSCG(Diagonal,src,summed);
  error = summed - combined;
  std::cout<<GridLogMessage << " summed-combined "<<norm2(error)    <<std::endl;


  src=1.0;
  Chebyshev<LatticeFermion> Cheby(0.1,40.0,200,InverseApproximation);

  std::cout<<GridLogMessage<<"Chebuy approx vector "<<std::endl;
  Cheby(Diagonal,src,combined);
  std::ofstream of("cheby");
  Cheby.csv(of);

  Diagonal.ApplyInverse(src,reference);
  error = reference - combined;

  std::cout<<GridLogMessage << "Chebyshev inverse test "<<std::endl;
  std::cout<<GridLogMessage << " Reference "<<norm2(reference)<<std::endl;
  std::cout<<GridLogMessage << " combined  "<<norm2(combined) <<std::endl;
  std::cout<<GridLogMessage << " error     "<<norm2(error)    <<std::endl;

  Grid_finalize();
}
