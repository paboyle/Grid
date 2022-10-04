    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/action/pseudofermion/GeneralEvenOddRationalRatioMixedPrec.h

    Copyright (C) 2015

    Author: Christopher Kelly <ckelly@bnl.gov>
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
#ifndef QCD_PSEUDOFERMION_GENERAL_EVEN_ODD_RATIONAL_RATIO_MIXED_PREC_H
#define QCD_PSEUDOFERMION_GENERAL_EVEN_ODD_RATIONAL_RATIO_MIXED_PREC_H

NAMESPACE_BEGIN(Grid);

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Generic rational approximation for ratios of operators utilizing the mixed precision multishift algorithm
    // cf. GeneralEvenOddRational.h for details
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
      
    template<class ImplD, class ImplF, class ImplD2>
    class GeneralEvenOddRatioRationalMixedPrecPseudoFermionAction : public GeneralEvenOddRatioRationalPseudoFermionAction<ImplD> {
    private:
      typedef typename ImplD2::FermionField FermionFieldD2;
      typedef typename ImplD::FermionField FermionFieldD;
      typedef typename ImplF::FermionField FermionFieldF;

      FermionOperator<ImplD> & NumOpD;
      FermionOperator<ImplD> & DenOpD;

      FermionOperator<ImplD2> & NumOpD2;
      FermionOperator<ImplD2> & DenOpD2;
     
      FermionOperator<ImplF> & NumOpF;
      FermionOperator<ImplF> & DenOpF;

      Integer ReliableUpdateFreq;
    protected:

      //Allow derived classes to override the multishift CG
      virtual void multiShiftInverse(bool numerator, const MultiShiftFunction &approx, const Integer MaxIter, const FermionFieldD &in, FermionFieldD &out){
#if 0
	SchurDifferentiableOperator<ImplD> schurOp(numerator ? NumOp : DenOp);
	ConjugateGradientMultiShift<FermionFieldD> msCG(MaxIter, approx);
	msCG(schurOp,in, out);
#else
	SchurDifferentiableOperator<ImplD2> schurOpD2(numerator ? NumOpD2 : DenOpD2);
	SchurDifferentiableOperator<ImplF> schurOpF(numerator ? NumOpF : DenOpF);
	FermionFieldD2 inD2(NumOpD2.FermionRedBlackGrid());
	FermionFieldD2 outD2(NumOpD2.FermionRedBlackGrid());
	
	ConjugateGradientMultiShiftMixedPrec<FermionFieldD2, FermionFieldF> msCG(MaxIter, approx, NumOpF.FermionRedBlackGrid(), schurOpF, ReliableUpdateFreq);
	precisionChange(inD2,in);
	std::cout << "msCG single solve "<<norm2(inD2)<<" " <<norm2(in)<<std::endl;
	msCG(schurOpD2, inD2, outD2);
	precisionChange(out,outD2);
#endif
      }
      virtual void multiShiftInverse(bool numerator, const MultiShiftFunction &approx, const Integer MaxIter, const FermionFieldD &in, std::vector<FermionFieldD> &out_elems, FermionFieldD &out){
	SchurDifferentiableOperator<ImplD2> schurOpD2(numerator ? NumOpD2 : DenOpD2);
	SchurDifferentiableOperator<ImplF> schurOpF(numerator ? NumOpF : DenOpF);

	FermionFieldD2 inD2(NumOpD2.FermionRedBlackGrid());
	FermionFieldD2 outD2(NumOpD2.FermionRedBlackGrid());
	std::vector<FermionFieldD2> out_elemsD2(out_elems.size(),NumOpD2.FermionRedBlackGrid());
	ConjugateGradientMultiShiftMixedPrec<FermionFieldD2, FermionFieldF> msCG(MaxIter, approx, NumOpF.FermionRedBlackGrid(), schurOpF, ReliableUpdateFreq);
	precisionChange(inD2,in);
	std::cout << "msCG in "<<norm2(inD2)<<" " <<norm2(in)<<std::endl;
	msCG(schurOpD2, inD2, out_elemsD2, outD2);
	precisionChange(out,outD2);
	for(int i=0;i<out_elems.size();i++){
	  precisionChange(out_elems[i],out_elemsD2[i]);
	}
      }
      //Allow derived classes to override the gauge import
      virtual void ImportGauge(const typename ImplD::GaugeField &Ud){

	typename ImplF::GaugeField Uf(NumOpF.GaugeGrid());
	typename ImplD2::GaugeField Ud2(NumOpD2.GaugeGrid());
	precisionChange(Uf, Ud);
	precisionChange(Ud2, Ud);

	std::cout << "Importing "<<norm2(Ud)<<" "<< norm2(Uf)<<" " << norm2(Ud2)<<std::endl;
	
	NumOpD.ImportGauge(Ud);
	DenOpD.ImportGauge(Ud);

	NumOpF.ImportGauge(Uf);
	DenOpF.ImportGauge(Uf);

	NumOpD2.ImportGauge(Ud2);
	DenOpD2.ImportGauge(Ud2);
      }
      
    public:
      GeneralEvenOddRatioRationalMixedPrecPseudoFermionAction(FermionOperator<ImplD>  &_NumOpD, FermionOperator<ImplD>  &_DenOpD, 
							      FermionOperator<ImplF>  &_NumOpF, FermionOperator<ImplF>  &_DenOpF, 
							      FermionOperator<ImplD2>  &_NumOpD2, FermionOperator<ImplD2>  &_DenOpD2, 
							      const RationalActionParams & p, Integer _ReliableUpdateFreq
							      ) : GeneralEvenOddRatioRationalPseudoFermionAction<ImplD>(_NumOpD, _DenOpD, p),
								  ReliableUpdateFreq(_ReliableUpdateFreq),
								  NumOpD(_NumOpD), DenOpD(_DenOpD),
								  NumOpF(_NumOpF), DenOpF(_DenOpF),
								  NumOpD2(_NumOpD2), DenOpD2(_DenOpD2)
      {}
      
      virtual std::string action_name(){return "GeneralEvenOddRatioRationalMixedPrecPseudoFermionAction";}
    };

NAMESPACE_END(Grid);

#endif
