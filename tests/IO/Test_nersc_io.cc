    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_nersc_io.cc

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
#include <Grid/Grid.h>

using namespace std;
using namespace Grid;
 ;


int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  std::cout <<GridLogMessage<< " main "<<std::endl;

  Coordinate simd_layout = GridDefaultSimd(4,vComplex::Nsimd());
  Coordinate mpi_layout  = GridDefaultMpi();
  //std::vector<int> latt_size  ({48,48,48,96});
  //std::vector<int> latt_size  ({32,32,32,32});
  Coordinate latt_size  ({16,16,16,32});
  Coordinate clatt_size  ({4,4,4,8});
  int orthodir=3;
  int orthosz =latt_size[orthodir];
    
  GridCartesian     Fine(latt_size,simd_layout,mpi_layout);
  GridCartesian     Coarse(clatt_size,simd_layout,mpi_layout);

  GridParallelRNG   pRNGa(&Fine);
  GridParallelRNG   pRNGb(&Fine);
  GridSerialRNG     sRNGa;
  GridSerialRNG     sRNGb;

  std::cout <<GridLogMessage<< " seeding... "<<std::endl;
  pRNGa.SeedFixedIntegers(std::vector<int>({45,12,81,9}));
  sRNGa.SeedFixedIntegers(std::vector<int>({45,12,81,9}));
  std::cout <<GridLogMessage<< " ...done "<<std::endl;

  std::string rfile("./ckpoint_rng.4000");
  FieldMetaData rngheader;
  NerscIO::writeRNGState(sRNGa,pRNGa,rfile);
  NerscIO::readRNGState (sRNGb,pRNGb,rngheader,rfile);

  LatticeComplex tmpa(&Fine); random(pRNGa,tmpa);
  LatticeComplex tmpb(&Fine); random(pRNGb,tmpb);
  tmpa = tmpa - tmpb;
  std::cout <<GridLogMessage<< " difference between restored randoms and orig "<<norm2( tmpa ) <<" / "<< norm2(tmpb)<<std::endl;

  ComplexD a,b;

  random(sRNGa,a);
  random(sRNGb,b);
  std::cout <<GridLogMessage<< " serial RNG numbers "<<a<<" "<<b<<std::endl;

  LatticeGaugeField Umu(&Fine);
  LatticeGaugeField Umu_diff(&Fine);
  LatticeGaugeField Umu_saved(&Fine);

  std::vector<LatticeColourMatrix> U(4,&Fine);
  
  SU3::HotConfiguration(pRNGa,Umu);

  FieldMetaData header;
  std::string file("./ckpoint_lat.4000");

  int precision32 = 0;
  int tworow      = 0;
  NerscIO::writeConfiguration(Umu,file,tworow,precision32);
  Umu_saved = Umu;
  NerscIO::readConfiguration(Umu,header,file);
  Umu_diff = Umu - Umu_saved;
  //std::cout << "Umu_save "<<Umu_saved[0]<<std::endl;
  //std::cout << "Umu_read "<<Umu[0]<<std::endl;
  std::cout <<GridLogMessage<< "norm2 Gauge Diff = "<<norm2(Umu_diff)<<std::endl;

  for(int mu=0;mu<Nd;mu++){
    U[mu] = PeekIndex<LorentzIndex>(Umu,mu);
  }

  // Painful ; fix syntactical niceness
  LatticeComplex LinkTrace(&Fine);
  LinkTrace=Zero();
  for(int mu=0;mu<Nd;mu++){
    LinkTrace = LinkTrace + trace(U[mu]);
  }

  // (1+2+3)=6 = N(N-1)/2 terms
  LatticeComplex Plaq(&Fine);
  LatticeComplex cPlaq(&Coarse);

  Plaq = Zero();
#if 1
  for(int mu=1;mu<Nd;mu++){
    for(int nu=0;nu<mu;nu++){
      Plaq = Plaq + trace(U[mu]*Cshift(U[nu],mu,1)*adj(Cshift(U[mu],nu,1))*adj(U[nu]));
    }
  }
#endif
  double vol = Fine.gSites();
  Complex PlaqScale(1.0/vol/6.0/3.0);

  std::vector<TComplex> Plaq_T(orthosz);
  sliceSum(Plaq,Plaq_T,Nd-1);
  int Nt = Plaq_T.size();

  TComplex Plaq_T_sum; 
  Plaq_T_sum=Zero();
  for(int t=0;t<Nt;t++){
    Plaq_T_sum = Plaq_T_sum+Plaq_T[t];
    Complex Pt=TensorRemove(Plaq_T[t]);
    std::cout<<GridLogMessage << "sliced ["<<t<<"]" <<Pt*PlaqScale*Real(Nt)<<std::endl;
  }

  {
    Complex Pt = TensorRemove(Plaq_T_sum);
    std::cout<<GridLogMessage << "total " <<Pt*PlaqScale<<std::endl;
  }  


  TComplex Tp = sum(Plaq);
  Complex p  = TensorRemove(Tp);
  std::cout<<GridLogMessage << "calculated plaquettes " <<p*PlaqScale<<std::endl;

  Complex LinkTraceScale(1.0/vol/4.0/3.0);
  TComplex Tl = sum(LinkTrace);
  Complex l  = TensorRemove(Tl);
  std::cout<<GridLogMessage << "calculated link trace " <<l*LinkTraceScale<<std::endl;

  blockSum(cPlaq,Plaq);
  TComplex TcP = sum(cPlaq);
  Complex ll= TensorRemove(TcP);
  std::cout<<GridLogMessage << "coarsened plaquettes sum to " <<ll*PlaqScale<<std::endl;

  std::string clone2x3("./ckpoint_clone2x3.4000");
  std::string clone3x3("./ckpoint_clone3x3.4000");

  NerscIO::writeConfiguration(Umu,clone3x3,0,precision32);
  NerscIO::writeConfiguration(Umu,clone2x3,1,precision32);
  
  Grid_finalize();
}
