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

//typedef typename GparityDomainWallFermionR::FermionField FermionField;
typedef typename ZMobiusFermionR::FermionField FermionField;

RealD AllZero(RealD x){ return 0.;}

class CmdJobParams 
{
  public:
    std::string gaugefile;

    int Ls;
    double mass;
    double M5;
    double mob_b;
    std::vector<ComplexD> omega;
    std::vector<ComplexD> boundary_phase;
    
    LanczosType Impl;
    int Nu;
    int Nk;
    int Np;
    int Nm;
    int Nstop;
    int Ntest;
    int MaxIter;
    double resid;
    
    double low;
    double high;
    int order;

    CmdJobParams()
      : gaugefile("Hot"),
        Ls(8), mass(0.01), M5(1.8), mob_b(1.5),
        Impl(LanczosType::irbl),
        Nu(4), Nk(200), Np(200), Nstop(100), Ntest(1), MaxIter(10), resid(1.0e-8), 
        low(0.2), high(5.5), order(11)
    {Nm=Nk+Np;};
    
    void Parse(char **argv, int argc);
};


void CmdJobParams::Parse(char **argv,int argc)
{
  std::string arg;
  std::vector<int> vi;
  double re,im;
  int expect, idx;
  std::string vstr;
  std::ifstream pfile;
  
  if( GridCmdOptionExists(argv,argv+argc,"--gconf") ){
    gaugefile = GridCmdOptionPayload(argv,argv+argc,"--gconf");
  }
  
  if( GridCmdOptionExists(argv,argv+argc,"--phase") ){
    arg = GridCmdOptionPayload(argv,argv+argc,"--phase");
    pfile.open(arg);
    assert(pfile);
    expect = 0;
    while( pfile >> vstr ) {
      if ( vstr.compare("boundary_phase") == 0 ) {
        pfile >> vstr;
        GridCmdOptionInt(vstr,idx);
        assert(expect==idx);
        pfile >> vstr;
        GridCmdOptionFloat(vstr,re);
        pfile >> vstr;
        GridCmdOptionFloat(vstr,im);
        boundary_phase.push_back({re,im});
        expect++;
      }
    }
    pfile.close();
  } else {
    for (int i=0; i<4; ++i) boundary_phase.push_back({1.,0.});
  }
  
  if( GridCmdOptionExists(argv,argv+argc,"--omega") ){
    arg = GridCmdOptionPayload(argv,argv+argc,"--omega");
    pfile.open(arg);
    assert(pfile);
    Ls = 0;
    while( pfile >> vstr ) {
      if ( vstr.compare("omega") == 0 ) {
        pfile >> vstr;
        GridCmdOptionInt(vstr,idx);
        assert(Ls==idx);
        pfile >> vstr;
        GridCmdOptionFloat(vstr,re);
        pfile >> vstr;
        GridCmdOptionFloat(vstr,im);
        omega.push_back({re,im});
        Ls++;
      }
    }
    pfile.close();
  } else {
    if( GridCmdOptionExists(argv,argv+argc,"--Ls") ){
      arg = GridCmdOptionPayload(argv,argv+argc,"--Ls");
      GridCmdOptionInt(arg,Ls);
    }
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
    // ypj[fixme] mode overriding message is needed.
    Impl = LanczosType::irbl;
    Nm = Nk+Np;
  }
  
  // block Lanczos with explicit extension of its dimensions
  if( GridCmdOptionExists(argv,argv+argc,"--rbl") ){
    arg = GridCmdOptionPayload(argv,argv+argc,"--rbl");
    GridCmdOptionIntVector(arg,vi);
    Nu = vi[0];
    Nk = vi[1];
    Np = vi[2]; // vector space is enlarged by adding Np vectors
    Nstop = vi[3];
    MaxIter = vi[4];
    // ypj[fixme] mode overriding message is needed.
    Impl = LanczosType::rbl;
    Nm = Nk+Np*MaxIter;
  }
  
  if( GridCmdOptionExists(argv,argv+argc,"--check_int") ){
    arg = GridCmdOptionPayload(argv,argv+argc,"--check_int");
    GridCmdOptionInt(arg,Ntest);
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
    std::streamsize ss = std::cout.precision();
    std::cout << GridLogMessage <<" Gauge Configuration "<< gaugefile << '\n';
    std::cout.precision(15);
    for ( int i=0; i<4; ++i ) std::cout << GridLogMessage <<" boundary_phase["<< i << "] = " << boundary_phase[i] << '\n';
    std::cout.precision(ss);
    std::cout << GridLogMessage <<" Ls "<< Ls << '\n';
    std::cout << GridLogMessage <<" mass "<< mass << '\n';
    std::cout << GridLogMessage <<" M5 "<< M5 << '\n';
    std::cout << GridLogMessage <<" mob_b "<< mob_b << '\n';
    std::cout.precision(15);
    for ( int i=0; i<Ls; ++i ) std::cout << GridLogMessage <<" omega["<< i << "] = " << omega[i] << '\n';
    std::cout.precision(ss);
    std::cout << GridLogMessage <<" Nu "<< Nu << '\n'; 
    std::cout << GridLogMessage <<" Nk "<< Nk << '\n'; 
    std::cout << GridLogMessage <<" Np "<< Np << '\n'; 
    std::cout << GridLogMessage <<" Nm "<< Nm << '\n'; 
    std::cout << GridLogMessage <<" Nstop "<< Nstop << '\n'; 
    std::cout << GridLogMessage <<" Ntest "<< Ntest << '\n'; 
    std::cout << GridLogMessage <<" MaxIter "<< MaxIter << '\n'; 
    std::cout << GridLogMessage <<" resid "<< resid << '\n'; 
    std::cout << GridLogMessage <<" Cheby Poly "<< low << "," << high << "," << order << std::endl; 
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
    // ypj [fixme] additional checks for the loaded configuration?
  }
  
  for(int mu=0;mu<Nd;mu++){
    U[mu] = PeekIndex<LorentzIndex>(Umu,mu);
  }
  
  RealD mass = JP.mass;
  RealD M5 = JP.M5;

// ypj [fixme] flexible support for a various Fermions
//  RealD mob_b = JP.mob_b;      // Gparity
//  std::vector<ComplexD> omega; // ZMobius
  
//  GparityMobiusFermionD ::ImplParams params;
//  std::vector<int> twists({1,1,1,0});
//  params.twists = twists;
//  GparityMobiusFermionR  Ddwf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,mob_b,mob_b-1.,params);
//  SchurDiagTwoOperator<GparityMobiusFermionR,FermionField> HermOp(Ddwf);

  //WilsonFermionR::ImplParams params;
  ZMobiusFermionR::ImplParams params;
  params.overlapCommsCompute = true;
  params.boundary_phases = JP.boundary_phase;
  ZMobiusFermionR  Ddwf(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,JP.omega,1.,0.,params);
  SchurDiagTwoOperator<ZMobiusFermionR,FermionField> HermOp(Ddwf);

  //std::vector<double> Coeffs { 0.,-1.}; 
  // ypj [note] this may not be supported by some compilers
  std::vector<double> Coeffs({ 0.,-1.}); 
  Polynomial<FermionField> PolyX(Coeffs);
  //Chebyshev<FermionField> Cheb(0.2,5.5,11);
  Chebyshev<FermionField> Cheb(JP.low,JP.high,JP.order);
//  Cheb.csv(std::cout);
  ImplicitlyRestartedBlockLanczos<FermionField> IRBL(HermOp,
                                                     Cheb,
                                                     JP.Nstop, JP.Ntest,
                                                     JP.Nu, JP.Nk, JP.Nm,
                                                     JP.resid,
                                                     JP.MaxIter);
  
  std::vector<RealD> eval(JP.Nm);
  
  std::vector<FermionField> src(JP.Nu,FrbGrid);
  for ( int i=0; i<JP.Nu; ++i ) gaussian(RNG5rb,src[i]);
  
  std::vector<FermionField> evec(JP.Nm,FrbGrid);
  for(int i=0;i<1;++i){
    std::cout << GridLogMessage << i <<" / "<< JP.Nm <<" grid pointer "<< evec[i]._grid << std::endl;
  };

  int Nconv;
  IRBL.calc(eval,evec,src,Nconv,JP.Impl);


  Grid_finalize();
}
