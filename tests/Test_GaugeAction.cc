    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_GaugeAction.cc

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

#include <qcd/utils/CovariantCshift.h>
#include <qcd/utils/WilsonLoops.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;


/* For Metropolis */
class Metropolis {
public:
  GridSerialRNG & sRNG;
  Metropolis(GridSerialRNG & _sRNG) : sRNG(_sRNG) {};
  bool AcceptReject(const RealD Delta)
  {
    RealD rand;

    if(Delta <=0.0) return true;

    random(sRNG,rand);
    if(rand <= exp(-Delta))
      return true;
    else 
      return false;
  }
};


int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);


  std::vector<int> simd_layout = GridDefaultSimd(4,vComplex::Nsimd());
  std::vector<int> mpi_layout  = GridDefaultMpi();
  std::vector<int> latt_size  ({16,16,16,32});
  std::vector<int> clatt_size  ({4,4,4,8});
  int orthodir=3;
  int orthosz =latt_size[orthodir];
    
  GridCartesian     Fine(latt_size,simd_layout,mpi_layout);
  GridCartesian     Coarse(clatt_size,simd_layout,mpi_layout);

  LatticeGaugeField Umu(&Fine);

  std::vector<LatticeColourMatrix> U(4,&Fine);
  
  NerscField header;
  
  std::string file("./ckpoint_lat.4000");
  NerscIO::readConfiguration(Umu,header,file);

  for(int mu=0;mu<Nd;mu++){
    U[mu] = PeekIndex<LorentzIndex>(Umu,mu);
  }

  // Painful ; fix syntactical niceness : to check reader
  LatticeComplex LinkTrace(&Fine);
  LinkTrace=zero;
  for(int mu=0;mu<Nd;mu++){
    LinkTrace = LinkTrace + trace(U[mu]);
  }

  // (1+2+3)=6 = N(N-1)/2 terms // equals to Double S Chorma
  // in LatticeGaugeField out Plaq
  // class WilsonLoop {
  // RealD  plaquette(LatticeGaugeField &Umu);
  // void   staple(LatticeSomethingorOther,LatticeGaugeField &Umu);
  // RealD  rectangle(LatticeGaugeField &Umu);
  // LatticeComplex sitePlaquette()
  // } 
  // covariantCshift ???
  // GaugeActionBase
  // GaugeActionPlaquette
  // GaugeActionPlaquettePlusRectangle
  // GaugeActionIwasaki
  // GaugeActionSymanzik
  // GaugeActionWilson
  // Heatbath and quenched update.
  //  
  LatticeColourMatrix tmpU(&Fine);

  LatticeComplex Plaq(&Fine);
  LatticeComplex cPlaq(&Coarse);
  Plaq = zero;
  for(int mu=1;mu<Nd;mu++){
    for(int nu=0;nu<mu;nu++){
      Plaq = Plaq + trace(PeriodicBC::CovShiftForward(U[mu],mu,U[nu])*adj(PeriodicBC::CovShiftForward(U[nu],nu,U[mu])));
    }
  }
  
  double vol = Fine.gSites();
  Complex PlaqScale(1.0/vol/6.0/3.0);
  RealD   StapScale(1.0/vol/6.0/3.0);
  std::cout<<GridLogMessage <<"PlaqScale" << PlaqScale<<std::endl;
  std::vector<TComplex> Plaq_T(orthosz);
  sliceSum(Plaq,Plaq_T,Nd-1);
  int Nt = Plaq_T.size();


  TComplex Plaq_T_sum; 
  Plaq_T_sum=zero;
  for(int t=0;t<Nt;t++){
    Plaq_T_sum = Plaq_T_sum+Plaq_T[t];
    Complex Pt=TensorRemove(Plaq_T[t]);
    std::cout<<GridLogMessage << "sliced ["<<t<<"]" <<Pt*PlaqScale*Real(Nt) << std::endl;
  }

  {
    Complex Pt = TensorRemove(Plaq_T_sum);
    std::cout<<GridLogMessage << "total " <<Pt*PlaqScale<<std::endl;
  }  

  TComplex Tp = sum(Plaq);
  Complex p  = TensorRemove(Tp);
  std::cout<<GridLogMessage << "calculated plaquettes " <<p*PlaqScale<<std::endl;

  RealD avg_plaq = ColourWilsonLoops::avgPlaquette(Umu);
  std::cout<<GridLogMessage << "NEW : calculated real plaquettes " <<avg_plaq<<std::endl;

  RealD stap_plaq=0.0;
  LatticeColourMatrix stap(&Fine);
  LatticeComplex stap_tr(&Fine);
  for(int mu=0;mu<Nd;mu++){
    ColourWilsonLoops::Staple(stap,Umu,mu);
    stap_tr = trace(stap*U[mu]);
    TComplex Ts = sum(stap_tr);
    Complex s  = TensorRemove(Ts);
    stap_plaq+=real(s);
  }
  std::cout<<GridLogMessage << "NEW : plaquette via staples"<< stap_plaq*StapScale*0.25<< std::endl;
  Complex LinkTraceScale(1.0/vol/4.0/3.0);
  TComplex Tl = sum(LinkTrace);
  Complex l  = TensorRemove(Tl);
  std::cout<<GridLogMessage << "calculated link trace " <<l*LinkTraceScale<<std::endl;

  blockSum(cPlaq,Plaq);
  TComplex TcP = sum(cPlaq);
  Complex ll= TensorRemove(TcP);
  std::cout<<GridLogMessage << "coarsened plaquettes sum to " <<ll*PlaqScale<<std::endl;


  Grid_finalize();
  }



