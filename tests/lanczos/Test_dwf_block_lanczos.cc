    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_dwf_block_lanczos.cc

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>

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
#include <Grid/Grid.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

typedef typename GparityDomainWallFermionR::FermionField FermionField;

RealD AllZero(RealD x){ return 0.;}

class CmdJobParams 
{
  public:
    std::string gaugefile;

    int Ls;
    double mass;
    double M5;
    double mob_b;
    
    int Nu;
    int Nk;
    int Np;
    int Nstop;
    int MaxIter;
    double resid;
    
    double low;
    double high;
    int order;

    CmdJobParams()
      : gaugefile("Hot"), 
        Ls(8), mass(0.01), M5(1.8), mob_b(1.5), 
        Nu(4), Nk(200), Np(200), Nstop(100), MaxIter(10), resid(1.0e-8), 
        low(0.2), high(5.5), order(11)
    {};
    
    void Parse(char **argv, int argc);
};


void CmdJobParams::Parse(char **argv,int argc)
{
  std::string arg;
  std::vector<int> vi;
  
  if( GridCmdOptionExists(argv,argv+argc,"--gconf") ){
    gaugefile = GridCmdOptionPayload(argv,argv+argc,"--gconf");
  }
  
  if( GridCmdOptionExists(argv,argv+argc,"--Ls") ){
    arg = GridCmdOptionPayload(argv,argv+argc,"--Ls");
    GridCmdOptionInt(arg,Ls);
  }
  
  if( GridCmdOptionExists(argv,argv+argc,"--mass") ){
    arg = GridCmdOptionPayload(argv,argv+argc,"--mass");
    GridCmdOptionFloat(arg,mass);
  }
  
  if( GridCmdOptionExists(argv,argv+argc,"--M5") ){
    arg = GridCmdOptionPayload(argv,argv+argc,"--M5");
    GridCmdOptionFloat(arg,M5);
  }
  
  if( GridCmdOptionExists(argv,argv+argc,"--mob_b") ){
    arg = GridCmdOptionPayload(argv,argv+argc,"--mob_b");
    GridCmdOptionFloat(arg,mob_b);
  }
  
  if( GridCmdOptionExists(argv,argv+argc,"--irbl") ){
    arg = GridCmdOptionPayload(argv,argv+argc,"--irbl");
    GridCmdOptionIntVector(arg,vi);
    Nu = vi[0];
    Nk = vi[1];
    Np = vi[2];
    Nstop = vi[3];
    MaxIter = vi[4];
  }
  
  if( GridCmdOptionExists(argv,argv+argc,"--resid") ){
    arg = GridCmdOptionPayload(argv,argv+argc,"--resid");
    GridCmdOptionFloat(arg,resid);
  }
  
  if( GridCmdOptionExists(argv,argv+argc,"--cheby_l") ){
    arg = GridCmdOptionPayload(argv,argv+argc,"--cheby_l");
    GridCmdOptionFloat(arg,low);
  }
  
  if( GridCmdOptionExists(argv,argv+argc,"--cheby_u") ){
    arg = GridCmdOptionPayload(argv,argv+argc,"--cheby_u");
    GridCmdOptionFloat(arg,high);
  }
  
  if( GridCmdOptionExists(argv,argv+argc,"--cheby_n") ){
    arg = GridCmdOptionPayload(argv,argv+argc,"--cheby_n");
    GridCmdOptionInt(arg,order);
  }
  
  if ( CartesianCommunicator::RankWorld() == 0 ) {
    clog <<" Gauge Configuration "<< gaugefile << '\n';
    clog <<" Ls "<< Ls << '\n';
    clog <<" mass "<< mass << '\n';
    clog <<" M5 "<< M5 << '\n';
    clog <<" mob_b "<< mob_b << '\n';
    clog <<" Nu "<< Nu << '\n'; 
    clog <<" Nk "<< Nk << '\n'; 
    clog <<" Np "<< Np << '\n'; 
    clog <<" Nstop "<< Nstop << '\n'; 
    clog <<" MaxIter "<< MaxIter << '\n'; 
    clog <<" resid "<< resid << '\n'; 
    clog <<" Cheby Poly "<< low << "," << high << "," << order << std::endl; 
  }
}


int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);
  
  CmdJobParams JP;
  JP.Parse(argv,argc);

  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(JP.Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(JP.Ls,UGrid);
  printf("UGrid=%p UrbGrid=%p FGrid=%p FrbGrid=%p\n",UGrid,UrbGrid,FGrid,FrbGrid);

  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});
  GridParallelRNG          RNG5(FGrid);  RNG5.SeedFixedIntegers(seeds5);
  GridParallelRNG          RNG4(UGrid);  RNG4.SeedFixedIntegers(seeds4);
  GridParallelRNG          RNG5rb(FrbGrid);  RNG5.SeedFixedIntegers(seeds5); 
  // ypj [note] why seed RNG5 again? bug? In this case, run with a default seed().
  //GridParallelRNG          RNG5rb(FrbGrid);  //RNG5rb.SeedFixedIntegers(seeds5);

  LatticeGaugeField Umu(UGrid); 
  std::vector<LatticeColourMatrix> U(4,UGrid);
  
  if ( JP.gaugefile.compare("Hot") == 0 ) {
    SU3::HotConfiguration(RNG4, Umu);
  } else {
    FieldMetaData header;
    NerscIO::readConfiguration(Umu,header,JP.gaugefile);
  }
  
  for(int mu=0;mu<Nd;mu++){
    U[mu] = PeekIndex<LorentzIndex>(Umu,mu);
  }
  
  RealD mass = JP.mass;
  RealD M5 = JP.M5;
  RealD mob_b = JP.mob_b;
//  DomainWallFermionR Ddwf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5);
  GparityMobiusFermionD ::ImplParams params;
  std::vector<int> twists({1,1,1,0});
  params.twists = twists;
  GparityMobiusFermionR  Ddwf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,mob_b,mob_b-1.,params);

//  MdagMLinearOperator<DomainWallFermionR,LatticeFermion> HermOp(Ddwf);
//  SchurDiagTwoOperator<DomainWallFermionR,LatticeFermion> HermOp(Ddwf);
  SchurDiagTwoOperator<GparityMobiusFermionR,FermionField> HermOp(Ddwf);
//  SchurDiagMooeeOperator<DomainWallFermionR,LatticeFermion> HermOp(Ddwf);

  int Nu = JP.Nu;
  int Nk = JP.Nk;
  int Nm = Nk+JP.Np;

  //std::vector<double> Coeffs { 0.,-1.}; 
  // ypj [note] this may not be supported by some compilers
  std::vector<double> Coeffs({ 0.,-1.}); 
  Polynomial<FermionField> PolyX(Coeffs);
  //Chebyshev<FermionField> Cheb(0.2,5.5,11);
  Chebyshev<FermionField> Cheb(JP.low,JP.high,JP.order);
//  Cheb.csv(std::cout);
  ImplicitlyRestartedBlockLanczos<FermionField> IRBL(HermOp,
                                                     Cheb,
                                                     JP.Nstop,
                                                     Nu,Nk,Nm,
                                                     JP.resid,
                                                     JP.MaxIter);
  
  std::vector<RealD> eval(Nm);
  
  std::vector<FermionField> src(Nu,FrbGrid);
  for ( int i=0; i<Nu; ++i ) gaussian(RNG5rb,src[i]);
  
  std::vector<FermionField> evec(Nm,FrbGrid);
  for(int i=0;i<1;++i){
    clog << i <<" / "<< Nm <<" grid pointer "<< evec[i]._grid << std::endl;
  };

  int Nconv;
  IRBL.calc(eval,evec,src,Nconv);


  Grid_finalize();
}
