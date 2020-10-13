/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/qdpxx/Test_qdpxx_wilson.cc

    Copyright (C) 2017

    Author: Felix Erben <felix.erben@ed.ac.uk>

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

#include <chroma.h>
#include <Grid/Grid.h>
#include <Grid/qcd/utils/BaryonUtils.h>

typedef Grid::LatticeGaugeField GaugeField;

namespace Chroma
{

class ChromaWrapper
{
public:
  typedef multi1d<LatticeColorMatrix> U;
  typedef LatticeFermion T4;

  static void ImportGauge(GaugeField &gr,
                          QDP::multi1d<QDP::LatticeColorMatrix> &ch)
  {
    Grid::LorentzColourMatrix LCM;
    Grid::Complex cc;
    QDP::ColorMatrix cm;
    QDP::Complex c;

    std::vector<int> x(4);
    QDP::multi1d<int> cx(4);
    Grid::Coordinate gd = gr.Grid()->GlobalDimensions();

    for (x[0] = 0; x[0] < gd[0]; x[0]++)
    {
      for (x[1] = 0; x[1] < gd[1]; x[1]++)
      {
        for (x[2] = 0; x[2] < gd[2]; x[2]++)
        {
          for (x[3] = 0; x[3] < gd[3]; x[3]++)
          {
            cx[0] = x[0];
            cx[1] = x[1];
            cx[2] = x[2];
            cx[3] = x[3];
            Grid::peekSite(LCM, gr, x);

            for (int mu = 0; mu < 4; mu++)
            {
              for (int i = 0; i < 3; i++)
              {
                for (int j = 0; j < 3; j++)
                {
                  cc = LCM(mu)()(i, j);
                  c = QDP::cmplx(QDP::Real(real(cc)), QDP::Real(imag(cc)));
                  QDP::pokeColor(cm, c, i, j);
                }
              }
              QDP::pokeSite(ch[mu], cm, cx);
            }
          }
        }
      }
    }
  }

  static void ExportGauge(GaugeField &gr,
                          QDP::multi1d<QDP::LatticeColorMatrix> &ch)
  {
    Grid::LorentzColourMatrix LCM;
    Grid::Complex cc;
    QDP::ColorMatrix cm;
    QDP::Complex c;

    std::vector<int> x(4);
    QDP::multi1d<int> cx(4);
    Grid::Coordinate gd = gr.Grid()->GlobalDimensions();

    for (x[0] = 0; x[0] < gd[0]; x[0]++)
    {
      for (x[1] = 0; x[1] < gd[1]; x[1]++)
      {
        for (x[2] = 0; x[2] < gd[2]; x[2]++)
        {
          for (x[3] = 0; x[3] < gd[3]; x[3]++)
          {
            cx[0] = x[0];
            cx[1] = x[1];
            cx[2] = x[2];
            cx[3] = x[3];

            for (int mu = 0; mu < 4; mu++)
            {
              for (int i = 0; i < 3; i++)
              {
                for (int j = 0; j < 3; j++)
                {
                  cm = QDP::peekSite(ch[mu], cx);
                  c = QDP::peekColor(cm, i, j);
                  cc = Grid::Complex(toDouble(real(c)), toDouble(imag(c)));
                  LCM(mu)
                  ()(i, j) = cc;
                }
              }
            }
            Grid::pokeSite(LCM, gr, x);
          }
        }
      }
    }
  }

  // Specific for Wilson Fermions
  static void ImportPropagator(Grid::LatticePropagator &gr,
                            QDP::LatticePropagator &ch)
  {
    Grid::LatticeSpinColourVector LF(gr.Grid());
    QDP::LatticeFermion cLF;

    int Nspin=4;
    int Ncolour=3;

    for (int is = 0; is < Nspin; is++){
      for (int ic = 0; ic < Ncolour; ic++){
	  Grid::PropToFerm<Grid::WilsonImplR>(LF,gr,is,ic);
          ImportFermion(LF,cLF);
          Chroma::FermToProp(cLF,ch,ic,is);
      }
    }
  }
  
  static void ExportPropagator(Grid::LatticePropagator &gr,
                            QDP::LatticePropagator &ch)
  {
    Grid::LatticeSpinColourVector LF(gr.Grid());
    QDP::LatticeFermion cLF;

    int Nspin=4;
    int Ncolour=3;

    for (int is = 0; is < Nspin; is++){
      for (int ic = 0; ic < Ncolour; ic++){
          Chroma::PropToFerm(ch,cLF,ic,is);
          ExportFermion(LF,cLF);
	  Grid::FermToProp<Grid::WilsonImplR>(gr,LF,is,ic);
      }
    }
  }

  // Specific for Wilson Fermions
  static void ImportFermion(Grid::LatticeFermion &gr,
                            QDP::LatticeFermion &ch)
  {
    Grid::SpinColourVector F;
    Grid::Complex c;

    QDP::Fermion cF;
    QDP::SpinVector cS;
    QDP::Complex cc;

    std::vector<int> x(4); // explicit 4d fermions in Grid
    QDP::multi1d<int> cx(4);
    Grid::Coordinate gd = gr.Grid()->GlobalDimensions();

    for (x[0] = 0; x[0] < gd[0]; x[0]++)
    {
      for (x[1] = 0; x[1] < gd[1]; x[1]++)
      {
        for (x[2] = 0; x[2] < gd[2]; x[2]++)
        {
          for (x[3] = 0; x[3] < gd[3]; x[3]++)
          {
            cx[0] = x[0];
            cx[1] = x[1];
            cx[2] = x[2];
            cx[3] = x[3];

            Grid::peekSite(F, gr, x);

            for (int j = 0; j < 3; j++)
            {
              for (int sp = 0; sp < 4; sp++)
              {

                c = F()(sp)(j);

                cc = QDP::cmplx(QDP::Real(real(c)), QDP::Real(imag(c)));

                QDP::pokeSpin(cS, cc, sp);
              }
              QDP::pokeColor(cF, cS, j);
            }
            QDP::pokeSite(ch, cF, cx);
          }
        }
      }
    }
  }

  // Specific for 4d Wilson fermions
  static void ExportFermion(Grid::LatticeFermion &gr,
                            QDP::LatticeFermion &ch)
  {
    Grid::SpinColourVector F;
    Grid::Complex c;

    QDP::Fermion cF;
    QDP::SpinVector cS;
    QDP::Complex cc;

    std::vector<int> x(4); // 4d fermions
    QDP::multi1d<int> cx(4);
    Grid::Coordinate gd = gr.Grid()->GlobalDimensions();

    for (x[0] = 0; x[0] < gd[0]; x[0]++)
    {
      for (x[1] = 0; x[1] < gd[1]; x[1]++)
      {
        for (x[2] = 0; x[2] < gd[2]; x[2]++)
        {
          for (x[3] = 0; x[3] < gd[3]; x[3]++)
          {
            cx[0] = x[0];
            cx[1] = x[1];
            cx[2] = x[2];
            cx[3] = x[3];

            cF = QDP::peekSite(ch, cx);
            for (int sp = 0; sp < 4; sp++)
            {
              for (int j = 0; j < 3; j++)
              {
                cS = QDP::peekColor(cF, j);
                cc = QDP::peekSpin(cS, sp);
                c = Grid::Complex(QDP::toDouble(QDP::real(cc)),
                                  QDP::toDouble(QDP::imag(cc)));
                F()
                (sp)(j) = c;
              }
            }
            Grid::pokeSite(F, gr, x);
          }
        }
      }
    }
  }

};
} // namespace Chroma

void make_gauge(GaugeField &Umu, Grid::LatticePropagator &q1,Grid::LatticePropagator &q2,Grid::LatticePropagator &q3)
{
  using namespace Grid;
  using namespace Grid::QCD;

  std::vector<int> seeds4({1, 2, 3, 4});

  Grid::GridCartesian *UGrid = (Grid::GridCartesian *)Umu.Grid();
  Grid::GridParallelRNG RNG4(UGrid);
  RNG4.SeedFixedIntegers(seeds4);
  Grid::SU<Nc>::HotConfiguration(RNG4, Umu);

  // Propagator
  Grid::gaussian(RNG4, q1);
  Grid::gaussian(RNG4, q2);
  Grid::gaussian(RNG4, q3);
}

void calc_chroma(GaugeField &lat, Grid::LatticePropagator &qU,Grid::LatticePropagator &qD,Grid::LatticePropagator &qS, std::vector<QDP::Complex> &res, std::string baryon)
{
  QDP::multi1d<QDP::LatticeColorMatrix> u(4);
  Chroma::ChromaWrapper::ImportGauge(lat, u);

  QDP::LatticePropagator check;
  QDP::LatticePropagator result;
  QDP::LatticePropagator  psiU;
  QDP::LatticePropagator  psiD;
  QDP::LatticePropagator  psiS;


  Chroma::ChromaWrapper::ImportPropagator(qU, psiU);
  Chroma::ChromaWrapper::ImportPropagator(qD, psiD);
  Chroma::ChromaWrapper::ImportPropagator(qS, psiS);

  if(0){
    std::cout << "Testing ImportPropagator(): "  << std::endl;
    Grid::GridCartesian *UGrid = (Grid::GridCartesian *)lat.Grid();
    std::vector<Grid::TComplex> buf;
    Grid::LatticeComplex tmp(UGrid);
    tmp = Grid::trace(qU);  
    Grid::sliceSum(tmp,buf,Grid::Nd-1);
    for (unsigned int t = 0; t < buf.size(); ++t)
    {
  	  std::cout << "Grid qU " << t << " " << Grid::TensorRemove(buf[t]) << std::endl;
    }
    QDP::LatticeComplex ctmp;
    ctmp = QDP::trace(psiU);
    Chroma::SftMom phases0(0,true,3); //How do I circumvent this? sliceSum equivalent?
    QDP::multi2d<DComplex> hsum0;
    hsum0 = phases0.sft(ctmp);
    for(int t = 0; t < phases0.numSubsets(); ++t){
  	  std::cout << "Chroma qU " << t << " " << hsum0[0][t] << std::endl;
    }
  }

  SpinMatrix C;
  SpinMatrix C_5;
  SpinMatrix C_4_5;
  SpinMatrix CG_1;
  SpinMatrix CG_2;
  SpinMatrix CG_3;
  SpinMatrix CG_4;

     
  SpinMatrix g_one = 1.0;
  //C = \gamma_2\gamma_4
  C = (Gamma(10)*g_one);
  
  //C_5 = C*gamma_5
  C_5 = (Gamma(5)*g_one);
           
  //C_4_5 = C*gamma_4*gamma_5
  C_4_5 = (Gamma(13)*g_one);
  
  //CG_1 = C*gamma_1
  CG_1 = (Gamma(11)*g_one);
                     
  //CG_2 = C*gamma_2
  CG_2 = (Gamma(8)*g_one);
                           
  //CG_3 = C*gamma_3
  CG_3 = (Gamma(14)*g_one);
  
  //CG_4 = C*gamma_4
  CG_4 = (Gamma(2)*g_one);
                                     
  // S_proj_unpol = (1/2)(1 + gamma_4)
  SpinMatrix S_proj_unpol = 0.5 * (g_one + (g_one * Gamma(8)));
  
  QDP::LatticeComplex b_prop;
  QDP::LatticePropagator di_quark;

  if(! baryon.compare("OmegaX")){
    // Omega_x - this esentially is degenerate (s C\gamma_1 s)s
    // C gamma_1 = Gamma(10) * Gamma(1) = Gamma(11)
    di_quark = QDP::quarkContract13(psiS * CG_1, CG_1 * psiS);
    b_prop = QDP::trace(S_proj_unpol * QDP::traceColor(psiS * QDP::traceSpin(di_quark)))
		      	     + 2.0 * QDP::trace(S_proj_unpol * QDP::traceColor(psiS * di_quark));
  } else if (! baryon.compare("OmegaY")){
    // Omega_x - this esentially is degenerate (s C\gamma_3 s)s
    // C gamma_1 = Gamma(10) * Gamma(2) = Gamma(8)
    di_quark = QDP::quarkContract13(psiS * CG_2, CG_2 * psiS);
    b_prop = QDP::trace(S_proj_unpol * QDP::traceColor(psiS * QDP::traceSpin(di_quark)))
		      	     + 2.0 * QDP::trace(S_proj_unpol * QDP::traceColor(psiS * di_quark));
  } else if (! baryon.compare("OmegaZ")){
    // Omega_x - this esentially is degenerate (s C\gamma_3 s)s
    // C gamma_1 = Gamma(10) * Gamma(4) = Gamma(14)
    di_quark = QDP::quarkContract13(psiS * CG_3, CG_3 * psiS);
    b_prop = QDP::trace(S_proj_unpol * QDP::traceColor(psiS * QDP::traceSpin(di_quark)))
		      	     + 2.0 * QDP::trace(S_proj_unpol * QDP::traceColor(psiS * di_quark));
  } else if (! baryon.compare("Proton")){
    // Proton - this esentially is degenerate (d C\gamma_5 u)u
// This is how the UKHadron code is written - diquarks are swapped when compared to coment above code.	  
    //di_quark = QDP::quarkContract13(psiU * C_5, C_5 * psiD);
    di_quark = QDP::quarkContract13(psiD * C_5, C_5 * psiU);
    b_prop = QDP::trace(S_proj_unpol * QDP::traceColor(psiU * QDP::traceSpin(di_quark)))
		      	     + QDP::trace(S_proj_unpol * QDP::traceColor(psiU * di_quark));
  } else if (! baryon.compare("Lambda")){
    // Lambda (octet) - This is the totally antisymmetric 
    // one from the middle of the octet
    // Lambda - (d C\gamma_5 s)u - (u C\gamma_5 s)d
    // This is given by: 
    // 1/3[ <us>d + <ds>u + 4<ud>s - (usd) - (dsu) + 2(sud) + 2(sdu) + 2(uds) + 2(dus) ]
   
/* This is how the UKHadron code is written - diquarks are swapped when compared to coments above code.	  
    // This gives <us>d - (usd) -- yes
    di_quark = QDP::quarkContract13(psiU * C_5, C_5 * psiS);
    b_prop = QDP::trace(S_proj_unpol * QDP::traceColor(psiD * QDP::traceSpin(di_quark)))
      	     - QDP::trace(S_proj_unpol * QDP::traceColor(psiD * di_quark));

    // This gives <ds>u - (dsu) -- yes
    di_quark = quarkContract13(psiD * C_5,C_5 * psiS);
    b_prop += QDP::trace(S_proj_unpol * QDP::traceColor(psiU * QDP::traceSpin(di_quark)))
       	      - QDP::trace(S_proj_unpol * QDP::traceColor(psiU * di_quark));

    // This gives 4<ud>s  -- yes
    di_quark = quarkContract13(psiU * C_5,C_5 * psiD);
    b_prop += 4.0 * QDP::trace(S_proj_unpol * QDP::traceColor(psiS * QDP::traceSpin(di_quark)));

    //This gives 2(sud) -- yes
    di_quark = quarkContract13(psiS * C_5,C_5 * psiU);
    b_prop += 2.0 * QDP::trace(S_proj_unpol * QDP::traceColor(psiD * di_quark));

    // This gives 2(sdu) -- yes
    di_quark = quarkContract13(psiS * C_5,C_5 * psiD);
    b_prop += 2.0 * QDP::trace(S_proj_unpol * QDP::traceColor(psiU * di_quark));

    // This gives 2(uds) -- yes
    di_quark = quarkContract13(psiU * C_5,C_5 * psiD);
    b_prop += 2.0 * QDP::trace(S_proj_unpol * QDP::traceColor(psiS * di_quark));

    // This gives 2(dus) -- yes
    di_quark = quarkContract13(psiD * C_5,C_5 * psiU);
    b_prop += 2.0 * QDP::trace(S_proj_unpol * QDP::traceColor(psiS * di_quark));*/

    // This gives <us>d - (usd) -- yes
    di_quark = QDP::quarkContract13(psiS * C_5, C_5 * psiU);
    b_prop = QDP::trace(S_proj_unpol * QDP::traceColor(psiD * QDP::traceSpin(di_quark)))
      	     - QDP::trace(S_proj_unpol * QDP::traceColor(psiD * di_quark));

    // This gives <ds>u - (dsu) -- yes
    di_quark = quarkContract13(psiS * C_5,C_5 * psiD);
    b_prop += QDP::trace(S_proj_unpol * QDP::traceColor(psiU * QDP::traceSpin(di_quark)))
       	      - QDP::trace(S_proj_unpol * QDP::traceColor(psiU * di_quark));

    // This gives 4<ud>s  -- yes
    di_quark = quarkContract13(psiD * C_5,C_5 * psiU);
    b_prop += 4.0 * QDP::trace(S_proj_unpol * QDP::traceColor(psiS * QDP::traceSpin(di_quark)));

    //This gives 2(sud) -- yes
    di_quark = quarkContract13(psiU * C_5,C_5 * psiS);
    b_prop += 2.0 * QDP::trace(S_proj_unpol * QDP::traceColor(psiD * di_quark));

    // This gives 2(sdu) -- yes
    di_quark = quarkContract13(psiD * C_5,C_5 * psiS);
    b_prop += 2.0 * QDP::trace(S_proj_unpol * QDP::traceColor(psiU * di_quark));

    // This gives 2(uds) -- yes
    di_quark = quarkContract13(psiD * C_5,C_5 * psiU);
    b_prop += 2.0 * QDP::trace(S_proj_unpol * QDP::traceColor(psiS * di_quark));

    // This gives 2(dus) -- yes
    di_quark = quarkContract13(psiU * C_5,C_5 * psiD);
    b_prop += 2.0 * QDP::trace(S_proj_unpol * QDP::traceColor(psiS * di_quark));
  } else {
    std::cout << "baryon not part of test " << std::endl;
    return;
  }
  std::cout<< "Chroma computing " << baryon << std::endl;

  Chroma::SftMom phases(0,true,3); //How do I circumvent this? sliceSum equivalent?
  QDP::multi2d<DComplex> hsum;
  hsum = phases.sft(b_prop);
  int length = phases.numSubsets();
  res.resize(length);
  for(int t = 0; t < length; ++t){
    res[t] = hsum[0][t]; //Should I test momentum?
  }

}


void calc_grid(Grid::LatticeGaugeField &Umu, Grid::LatticePropagator &qU, Grid::LatticePropagator &qD, Grid::LatticePropagator &qS, std::vector<Grid::Complex> &res, std::string baryon)
{
  using namespace Grid;
  using namespace Grid::QCD;

  Grid::GridCartesian *UGrid = (Grid::GridCartesian *)Umu.Grid();

  Grid::Gamma G_A = Grid::Gamma(Grid::Gamma::Algebra::Identity);
  Grid::Gamma G_B = Grid::Gamma(Grid::Gamma::Algebra::GammaZGamma5); // OmegaX: C*GammaX = i* GammaZ*Gamma5

  Grid::LatticeComplex c(UGrid);
  Grid::LatticeComplex c1(UGrid);

  if(! baryon.compare("OmegaX")){
    BaryonUtils<Grid::WilsonImplR>::ContractBaryons(qS,qS,qS,G_A,G_B,G_A,G_B,"sss","sss",1,c);  
    c*=0.5;
    std::cout << "Grid-Omega factor 2 larger than Chroma-Omega!!!" << std::endl; 
  } else if (! baryon.compare("OmegaY")){
    G_B = Grid::Gamma(Grid::Gamma::Algebra::GammaT);
    BaryonUtils<Grid::WilsonImplR>::ContractBaryons(qS,qS,qS,G_A,G_B,G_A,G_B,"sss","sss",1,c);  
    c*=0.5;
    std::cout << "Grid-Omega factor 2 larger than Chroma-Omega!!!" << std::endl; 
  } else if (! baryon.compare("OmegaZ")){
    G_B = Grid::Gamma(Grid::Gamma::Algebra::GammaXGamma5);
    BaryonUtils<Grid::WilsonImplR>::ContractBaryons(qS,qS,qS,G_A,G_B,G_A,G_B,"sss","sss",1,c);  
    c*=0.5;
    std::cout << "Grid-Omega factor 2 larger than Chroma-Omega!!!" << std::endl; 
  } else if (! baryon.compare("Proton")){
    G_B = Grid::Gamma(Grid::Gamma::Algebra::SigmaXZ);
    BaryonUtils<Grid::WilsonImplR>::ContractBaryons(qU,qD,qU,G_A,G_B,G_A,G_B,"udu","udu",1,c); 
    std::cout << "UKHadron-Proton has flipped diquarks in original code." << std::endl; 
  } else if (! baryon.compare("Lambda")){
    G_B = Grid::Gamma(Grid::Gamma::Algebra::SigmaXZ);
    BaryonUtils<Grid::WilsonImplR>::ContractBaryons(qS,qU,qD,G_A,G_B,G_A,G_B,"sud","sud",1,c1); //<ud>s 
    c = 4.*c1; 
    BaryonUtils<Grid::WilsonImplR>::ContractBaryons(qD,qU,qS,G_A,G_B,G_A,G_B,"dus","dus",1,c1); //<us>d 
    c += c1; 
    BaryonUtils<Grid::WilsonImplR>::ContractBaryons(qU,qD,qS,G_A,G_B,G_A,G_B,"uds","uds",1,c1); //<ds>u  
    c += c1; 
    BaryonUtils<Grid::WilsonImplR>::ContractBaryons(qD,qU,qS,G_A,G_B,G_A,G_B,"dus","sud",1,c1); //(sud) 
    c += 2.*c1; 
    BaryonUtils<Grid::WilsonImplR>::ContractBaryons(qU,qD,qS,G_A,G_B,G_A,G_B,"uds","sud",1,c1); //(sdu) 
    c -= 2.*c1; 
    BaryonUtils<Grid::WilsonImplR>::ContractBaryons(qS,qU,qD,G_A,G_B,G_A,G_B,"sud","dus",1,c1); //(dus) 
    c += 2.*c1; 
    BaryonUtils<Grid::WilsonImplR>::ContractBaryons(qU,qD,qS,G_A,G_B,G_A,G_B,"uds","dus",1,c1); //-(dsu) 
    c -= c1; 
    BaryonUtils<Grid::WilsonImplR>::ContractBaryons(qS,qU,qD,G_A,G_B,G_A,G_B,"sud","uds",1,c1); //(uds) 
    c -= 2.*c1; 
    BaryonUtils<Grid::WilsonImplR>::ContractBaryons(qD,qU,qS,G_A,G_B,G_A,G_B,"dus","uds",1,c1); //-(usd) 
    c -= c1; 
    std::cout << "UKHadron-Lambda has flipped diquarks in original code." << std::endl; 
  } else {
    std::cout << "baryon not part of test " << std::endl;
    return;
  }
  std::cout<< "Grid computing " << baryon << std::endl;

  std::vector<Grid::TComplex> buf;
  Grid::sliceSum(c,buf,Grid::Nd-1);
  res.resize(buf.size());
  for (unsigned int t = 0; t < buf.size(); ++t)
  {
    res[t]=Grid::TensorRemove(buf[t]);
  }

}

int main(int argc, char **argv)
{

  /********************************************************
   * Setup QDP
   *********************************************************/
  Chroma::initialize(&argc, &argv);
  Chroma::WilsonTypeFermActs4DEnv::registerAll();

  /********************************************************
   * Setup Grid
   *********************************************************/
  Grid::Grid_init(&argc, &argv);
  Grid::GridCartesian *UGrid = Grid::SpaceTimeGrid::makeFourDimGrid(Grid::GridDefaultLatt(),
                                                                         Grid::GridDefaultSimd(Grid::Nd, Grid::vComplex::Nsimd()),
                                                                         Grid::GridDefaultMpi());

  Grid::Coordinate gd = UGrid->GlobalDimensions();
  QDP::multi1d<int> nrow(QDP::Nd);
  for (int mu = 0; mu < 4; mu++)
    nrow[mu] = gd[mu];

  QDP::Layout::setLattSize(nrow);
  QDP::Layout::create();

  GaugeField Ug(UGrid);
  typedef Grid::LatticePropagator PropagatorField;
  PropagatorField up(UGrid);
  PropagatorField down(UGrid);
  PropagatorField strange(UGrid);
  std::vector<ComplexD> res_chroma;
  std::vector<Grid::Complex> res_grid;
  Grid::Complex res_chroma_g;

  std::vector<std::string> baryons({"OmegaX","OmegaY","OmegaZ","Proton","Lambda"});
  int nBaryon=baryons.size();

  for (int iB = 0; iB < nBaryon; iB++)
  {
      make_gauge(Ug, up, down, strange); // fills the gauge field and the propagator with random numbers

      calc_chroma(Ug, up, down, strange, res_chroma,baryons[iB]);
 
      for(int t=0;t<res_chroma.size();t++){
        std::cout << " Chroma baryon "<<t<<" "<< res_chroma[t] << std::endl;
      }
         
      calc_grid(Ug, up, down, strange, res_grid,baryons[iB]);

      for(int t=0;t<res_chroma.size();t++){
        std::cout << " Grid baryon "<<t<<" "<< res_grid[t] << std::endl;
      }
      for(int t=0;t<res_chroma.size();t++){
	res_chroma_g = Grid::Complex(toDouble(real(res_chroma[t])), toDouble(imag(res_chroma[t])));
        std::cout << " Difference "<<t<<" "<< res_chroma_g - res_grid[t] << std::endl;
      }

      std::cout << "Finished test " << std::endl;

  }
  Chroma::finalize();
}
