/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./benchmarks/Benchmark_dwf_lex.cc

    Copyright (C) 2026

Author: Peter Boyle <pboyle@bnl.gov>

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

//
// Single core domain wall Dhop: vectorised chart against lexicographic
// (Nsimd()==1), identical gauge field and source, in one process.
//
// Run at both precisions.  Where the vectorised chart also has Nsimd()==1 the
// difference is instruction quality and register pressure alone; where it has
// lanes, the excess over that is the layout.
//
#include <Grid/Grid.h>

using namespace Grid;

template<class vobjOut,class vobjIn>
void transfer(Lattice<vobjOut> &out,const Lattice<vobjIn> &in)
{
  typedef typename vobjIn::scalar_object sobj;
  std::vector<sobj> buf;
  unvectorizeToLexOrdArray(buf,in);
  vectorizeFromLexOrdArray(buf,out);
}

template<class Op,class Field>
double TimeDhop(Op &D,Field &src,Field &res,int ncall)
{
  D.Dhop(src,res,DaggerNo);
  double t0 = usecond();
  for(int i=0;i<ncall;i++){
    D.Dhop(src,res,DaggerNo);
  }
  double t1 = usecond();
  return (t1-t0)/ncall;
}

template<class vImpl,class lexImpl>
void Sweep(std::string name,int Ls,Coordinate mpi)
{
  typedef DomainWallFermion<vImpl>   vDwf;
  typedef DomainWallFermion<lexImpl> lexDwf;
  typedef typename vImpl::Simd       vSimd;

  const long unsigned int single_site_flops = 8*Nc*(7+16*Nc);

  Coordinate vsimd = GridDefaultSimd(Nd,vSimd::Nsimd());
  Coordinate lsimd({1,1,1,1});

  std::cout << GridLogMessage << "==== " << name
	    << "  Nsimd(simd chart) = " << vSimd::Nsimd()
	    << "  Nsimd(lex chart) = 1  Ls = " << Ls << std::endl;
  std::cout << GridLogMessage
	    << "  L    localvol    MB/field    simd us     lex us    simd Gflop/s  lex Gflop/s  ratio"
	    << std::endl;

  // L divisible by 4: red-black needs an even reduced dimension, and the simd
  // layout can take 2 in one direction.
  for(auto L : std::vector<int>({4,8,12,16,20})){

    Coordinate latt({L,L,L,L});

    GridCartesian         *vU  = new GridCartesian(latt,vsimd,mpi);
    GridRedBlackCartesian *vUr = SpaceTimeGrid::makeFourDimRedBlackGrid(vU);
    GridCartesian         *vF  = SpaceTimeGrid::makeFiveDimGrid(Ls,vU);
    GridRedBlackCartesian *vFr = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,vU);

    GridCartesian         *lU  = new GridCartesian(latt,lsimd,mpi);
    GridRedBlackCartesian *lUr = SpaceTimeGrid::makeFourDimRedBlackGrid(lU);
    GridCartesian         *lF  = SpaceTimeGrid::makeFiveDimGrid(Ls,lU);
    GridRedBlackCartesian *lFr = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,lU);

    GridParallelRNG pRNG4(vU); pRNG4.SeedFixedIntegers(std::vector<int>({1,2,3,4}));
    GridParallelRNG pRNG5(vF); pRNG5.SeedFixedIntegers(std::vector<int>({5,6,7,8}));

    typename vDwf::GaugeField   Umu(vU);  SU<Nc>::HotConfiguration(pRNG4,Umu);
    typename lexDwf::GaugeField Ulex(lU); transfer(Ulex,Umu);

    typename vDwf::FermionField   vsrc(vF),vres(vF);  random(pRNG5,vsrc);
    typename lexDwf::FermionField lsrc(lF),lres(lF);  transfer(lsrc,vsrc);

    vDwf   Dv(Umu ,*vF,*vFr,*vU,*vUr,0.05,1.8);
    lexDwf Dl(Ulex,*lF,*lFr,*lU,*lUr,0.05,1.8);

    double volume=Ls; for(int mu=0;mu<Nd;mu++) volume=volume*latt[mu];

    int ncall = (L<=8) ? 200 : 50;
    double tv=1e30,tl=1e30;
    for(int rep=0;rep<3;rep++){
      tv = std::min(tv,TimeDhop(Dv,vsrc,vres,ncall));
      tl = std::min(tl,TimeDhop(Dl,lsrc,lres,ncall));
    }

    double gf_v = single_site_flops*volume/tv/1.0e3;
    double gf_l = single_site_flops*volume/tl/1.0e3;
    double mb   = volume*sizeof(typename vDwf::FermionField::scalar_object)/1024./1024.;

    std::cout << GridLogMessage
	      << std::setw(3) << L
	      << std::setw(12) << (long)volume
	      << std::setw(11) << std::fixed << std::setprecision(2) << mb
	      << std::setw(11) << std::setprecision(1) << tv
	      << std::setw(11) << tl
	      << std::setw(13) << std::setprecision(3) << gf_v
	      << std::setw(13) << gf_l
	      << std::setw(8)  << std::setprecision(2) << gf_v/gf_l
	      << std::endl;

    delete vFr; delete vF; delete vUr; delete vU;
    delete lFr; delete lF; delete lUr; delete lU;
  }
}

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  const int Ls = 8;
  Coordinate mpi = GridDefaultMpi();

  Sweep<WilsonImplD,lexWilsonImplD>("DOUBLE",Ls,mpi);
  Sweep<WilsonImplF,lexWilsonImplF>("SINGLE",Ls,mpi);

  Grid_finalize();
}
