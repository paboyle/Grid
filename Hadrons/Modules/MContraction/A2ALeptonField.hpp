/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MContraction/A2ALeptonField.cc

Copyright (C) 2015-2019

Author: Vera Guelpers <Vera.Guelpers@ed.ac.uk>

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
#ifndef Hadrons_MContraction_A2ALeptonField_hpp_
#define Hadrons_MContraction_A2ALeptonField_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/A2AMatrix.hpp>

#ifndef MLF_IO_TYPE
#define MLF_IO_TYPE ComplexF
#endif

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         A2ALeptonField                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class A2ALeptonFieldPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2ALeptonFieldPar,
                                    int, cacheBlock,
                                    int, block,
                                    std::string, left,
                                    std::string, right,
                                    std::string, output,
                                    std::vector<std::string>, lepton);
};

class A2ALeptonFieldMetadata: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2ALeptonFieldMetadata,
                                    std::string, leptonName,
				    int, sidx1,
				    int, sidx2);
};

template <typename T, typename FImpl>
class LeptonFieldKernel: public A2AKernel<T, typename FImpl::FermionField>
{
public:
    typedef typename FImpl::FermionField FermionField;
public:
    LeptonFieldKernel(const std::vector<LatticeComplex> &lepton02,
                      const std::vector<LatticeComplex> &lepton03,
                      const std::vector<LatticeComplex> &lepton12,
                      const std::vector<LatticeComplex> &lepton13,
                      GridBase *grid)
    : lepton02_(lepton02), lepton03_(lepton03), lepton12_(lepton12), lepton13_(lepton13), grid_(grid)
    {
        vol_ = 1.;
        for (auto &d: grid_->GlobalDimensions())
        {
            vol_ *= d;
        }
    }

    virtual ~LeptonFieldKernel(void) = default;
    virtual void operator()(A2AMatrixSet<T> &m, const FermionField *left, 
                            const FermionField *right,
                            const unsigned int orthogDim, double &t)
    {
        A2Autils<FImpl>::LeptonField(m, left, right, lepton02_, lepton03_, lepton12_, lepton13_, orthogDim, &t);
    }

    virtual double flops(const unsigned int blockSizei, const unsigned int blockSizej)
    {
        return 0.;
    }

    virtual double bytes(const unsigned int blockSizei, const unsigned int blockSizej)
    {
        return 0.;
    }
private:
    const std::vector<LatticeComplex> &lepton02_, &lepton03_, &lepton12_, &lepton13_;
    GridBase                          *grid_;
    double                            vol_;
};

template <typename FImpl, typename PhotonImpl>
class TA2ALeptonField: public Module<A2ALeptonFieldPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    typedef A2AMatrixBlockComputation<Complex, 
                                      FermionField, 
                                      A2ALeptonFieldMetadata, 
                                      MLF_IO_TYPE> Computation;
    typedef LeptonFieldKernel<Complex, FImpl> Kernel;
public:
    // constructor
    TA2ALeptonField(const std::string name);
    // destructor
    virtual ~TA2ALeptonField(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(A2ALeptonField, ARG(TA2ALeptonField<FIMPL, PhotonR>), MContraction);

/******************************************************************************
 *                 TA2ALeptonField implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename PhotonImpl>
TA2ALeptonField<FImpl, PhotonImpl>::TA2ALeptonField(const std::string name)
: Module<A2ALeptonFieldPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename PhotonImpl>
std::vector<std::string> TA2ALeptonField<FImpl, PhotonImpl>::getInput(void)
{
    std::vector<std::string> in = par().lepton;
    in.push_back(par().left);
    in.push_back(par().right);

    return in;
}

template <typename FImpl, typename PhotonImpl>
std::vector<std::string> TA2ALeptonField<FImpl, PhotonImpl>::getOutput(void)
{
    std::vector<std::string> out = {};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, typename PhotonImpl>
void TA2ALeptonField<FImpl, PhotonImpl>::setup(void)
{
    envTmp(Computation, "computation", 1, envGetGrid(FermionField), 
           env().getNd() - 1, par().lepton.size()*8, 1, par().block, 
           par().cacheBlock, this);
    envTmp(std::vector<LatticeComplex>, "Lmu", 1, 
           4, envGetGrid(LatticeComplex));
    envTmp(std::vector<LatticeComplex>, "L02", 1, 
           par().lepton.size()*8, envGetGrid(LatticeComplex));
    envTmp(std::vector<LatticeComplex>, "L03", 1, 
           par().lepton.size()*8, envGetGrid(LatticeComplex)); 
    envTmp(std::vector<LatticeComplex>, "L12", 1, 
           par().lepton.size()*8, envGetGrid(LatticeComplex));
    envTmp(std::vector<LatticeComplex>, "L13", 1, 
           par().lepton.size()*8, envGetGrid(LatticeComplex));
    envTmpLat(PropagatorField, "prop_buf");

}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, typename PhotonImpl>
void TA2ALeptonField<FImpl, PhotonImpl>::execute(void)
{
    auto &left  = envGet(std::vector<FermionField>, par().left);
    auto &right = envGet(std::vector<FermionField>, par().right);

    int nt         = env().getDim().back();
    int N_i        = left.size();
    int N_j        = right.size();
    int block      = par().block;
    int cacheBlock = par().cacheBlock;

    LOG(Message) << "Computing all-to-all lepton insertion fields" << std::endl;
    LOG(Message) << "Left: '" << par().left << "' Right: '" << par().right << "'" << std::endl;
    LOG(Message) << "leptons: " << std::endl;
    for (auto &name: par().lepton)
    {
        LOG(Message) << "  " << name << std::endl;
    }
    LOG(Message) << "lepton insertion field size: " << nt << "*" << N_i << "*" << N_j 
                 << " (filesize " << sizeString(nt*N_i*N_j*sizeof(MLF_IO_TYPE)) 
                 << "/spin index/lepton)" << std::endl;

    envGetTmp(std::vector<ComplexField>, L02);
    envGetTmp(std::vector<ComplexField>, L03);
    envGetTmp(std::vector<ComplexField>, L12);
    envGetTmp(std::vector<ComplexField>, L13);

    unsigned int ii = 0;
    envGetTmp(PropagatorField, prop_buf);
    envGetTmp(std::vector<LatticeComplex>, Lmu);
    for (unsigned int i = 0; i < par().lepton.size(); ++i)
    {
	auto &lepton = envGet(PropagatorField, par().lepton[i]);
	// need only half of the spin indices (s1s2) of the lepton field
	// the other ones are zero because of the left-handed current
	for (unsigned int s1 = 0; s1 < 2; ++s1) 
	for (unsigned int s2 = 0; s2 < 4; ++s2)
	{  	
	    for (unsigned int mu = 0; mu < 4; ++mu)
	    {
		prop_buf = GammaL(Gamma::gmu[mu]) * lepton;
		//lepon is unit matix in color space, so just pick one diagonal enty
		Lmu[mu] = peekColour(peekSpin(prop_buf,s1,s2),0,0);
	    }
	    // build the required combinations of lepton fields for the a2a kernel
	    // [g^L_mu * L_mu]_{ab} = L_{ab} (with mu = (x,y,z,t))
	    L02[ii] = 2.0*Lmu[3] + 2.0*timesI(Lmu[2]);
	    L03[ii] = 2.0*timesI(Lmu[0]) - 2.0*Lmu[1];
	    L12[ii] = 2.0*timesI(Lmu[0]) + 2.0*Lmu[1];
	    L13[ii] = 2.0*Lmu[3] - 2.0*timesI(Lmu[2]);
	    ii++;
    	}
    }

    auto ionameFn = [this](const unsigned int index, const unsigned int dummy)
    {
	//index = 8*l + 4*sindex1 + sindex2
	unsigned int sindex = index % 8;
	unsigned int sindex2 = sindex % 4;
 	unsigned int sindex1 = (sindex - sindex2)/4;
	unsigned int l = (index - sindex)/8;
	return par().lepton[l] + "_" + std::to_string(sindex1) + std::to_string(sindex2);
    };

    auto filenameFn = [this, &ionameFn](const unsigned int index, const unsigned int dummy)
    {
        return par().output + "." + std::to_string(vm().getTrajectory()) 
               + "/" + ionameFn(index, dummy) + ".h5";
    };

    auto metadataFn = [this](const unsigned int index, const unsigned int dummy)
    {
	A2ALeptonFieldMetadata md;

	unsigned int sindex = index % 8;
	unsigned int sindex2 = sindex % 4;
 	unsigned int sindex1 = (sindex - sindex2)/4;
	unsigned int l = (index - sindex)/8;
        md.leptonName = par().lepton[l];
	md.sidx1 = sindex1;
	md.sidx2 = sindex2;

	return md;
    };

    // executing computation
    Kernel kernel(L02, L03, L12, L13, envGetGrid(FermionField));
    envGetTmp(Computation, computation);
    computation.execute(left, right, kernel, ionameFn, filenameFn, metadataFn);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_A2ALeptonField_hpp_
