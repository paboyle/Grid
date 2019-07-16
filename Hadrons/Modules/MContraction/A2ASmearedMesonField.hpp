/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: A2ASmearedMesonField.hpp

Copyright (C) 2015-2019

Author: Antonin Portelli <antonin.portelli@me.com>
Author: Nils Asmussen <n.asmussen@soton.ac.uk>
Author: paboyle <paboyle@ph.ed.ac.uk>
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
#ifndef Hadrons_MContraction_A2ASmearedMesonField_hpp_
#define Hadrons_MContraction_A2ASmearedMesonField_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/A2AMatrix.hpp>
#include <map>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                     All-to-all meson field creation                        *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class A2ASmearedMesonFieldPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2ASmearedMesonFieldPar,
                                    int, cacheBlock,
                                    int, block,
                                    std::string, left,
                                    std::string, right,
                                    std::string, distributions,
                                    std::string, output,
                                    std::string, gammas,
                                    std::vector<std::string>, smearMom,
                                    std::vector<std::string>, mom);
};

class A2ASmearedMesonFieldMetadata: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2ASmearedMesonFieldMetadata,
                                    std::vector<int>, smearMom,
                                    std::vector<int>, mom,
                                    Gamma::Algebra, gamma);
};

template <typename T, typename FImpl>
class SmearedMesonFieldKernel: public A2AKernel<T, typename FImpl::FermionField>
{
public:
    typedef typename FImpl::FermionField FermionField;
public:
    SmearedMesonFieldKernel(const std::vector<Gamma::Algebra> &gamma,
                     const std::vector<LatticeComplex> &mom,
                     GridBase *grid)
    : gamma_(gamma), mom_(mom), grid_(grid)
    {
        vol_ = 1.;
        for (auto &d: grid_->GlobalDimensions())
        {
            vol_ *= d;
        }
    }

    virtual ~SmearedMesonFieldKernel(void) = default;
    virtual void operator()(A2AMatrixSet<T> &m, const FermionField *left,
                            const FermionField *right,
                            const unsigned int orthogDim, double &t)
    {
        A2Autils<FImpl>::MesonField(m, left, right, gamma_, mom_, orthogDim, &t);
    }

    virtual double flops(const unsigned int blockSizei, const unsigned int blockSizej)
    {
        return vol_*(2*8.0+6.0+8.0*mom_.size())*blockSizei*blockSizej*gamma_.size();
    }

    virtual double bytes(const unsigned int blockSizei, const unsigned int blockSizej)
    {
        return vol_*(12.0*sizeof(T))*blockSizei*blockSizej
               +  vol_*(2.0*sizeof(T)*mom_.size())*blockSizei*blockSizej*gamma_.size();
    }
private:
    const std::vector<Gamma::Algebra> &gamma_;
    const std::vector<LatticeComplex> &mom_;
    GridBase                          *grid_;
    double                            vol_;
};

template <typename FImpl>
class TA2ASmearedMesonField : public Module<A2ASmearedMesonFieldPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    typedef typename FImpl::GaugeLinkField GaugeMat;
    typedef A2AMatrixBlockComputation<Complex,
                                      FermionField,
                                      A2ASmearedMesonFieldMetadata,
                                      HADRONS_A2AM_IO_TYPE> Computation;
    typedef SmearedMesonFieldKernel<Complex, FImpl> Kernel;
    typedef std::pair<std::string, std::string> stringPair;
public:
    // constructor
    TA2ASmearedMesonField(const std::string name);
    // destructor
    virtual ~TA2ASmearedMesonField(void){};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    std::vector<Gamma::Algebra>       gamma_;
    std::vector<std::vector<int>>     smearMom_;
    std::vector<std::vector<int>>     mom_;
    bool hasDistributions_{false};
    std::string distributionsCache_;
    std::vector<stringPair> distributionsNames_;
    std::map<std::string, int> distributionsMap_;
    template<typename Lattice, typename Iterable>
    void multidim_Cshift_inplace(Lattice &lat, const Iterable &shifts);
    void smearing_weight(std::vector<ComplexField> &out,
            const std::vector<stringPair> &distributionsNames_,
            const std::map<std::string, int> &distributionsMap_,
            const std::vector<ComplexField> &distrib,
            const std::vector<std::vector<int>> &smearMom,
            const std::vector<int> &mom);
};

MODULE_REGISTER(A2ASmearedMesonField, ARG(TA2ASmearedMesonField<FIMPL>), MContraction);

/******************************************************************************
*                  TA2ASmearedMesonField implementation                       *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TA2ASmearedMesonField<FImpl>::TA2ASmearedMesonField(const std::string name)
: Module<A2ASmearedMesonFieldPar>(name)
{
}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TA2ASmearedMesonField<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().left, par().right};

    if (par().distributions.find('(') != std::string::npos)
    { //pairs
        using std::get;
        auto dists = strToVec<stringPair>(par().distributions);
        for (const auto &i : dists)
        {
            in.push_back(get<0>(i));
            in.push_back(get<1>(i));
        }
    }
    else
    { //no pairs
        auto dists = strToVec<std::string>(par().distributions);
        for (const auto &i : dists)
        {
            in.push_back(i);
        }
    }

    //remove duplicates
    {
        std::set<std::string> set;
        for(auto &i: in)
        {
            set.insert(i);
        }
        in.clear();
        for(auto &i: set)
        {
            in.push_back(i);
        }
    }

    return in;
}

template <typename FImpl>
std::vector<std::string> TA2ASmearedMesonField<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2ASmearedMesonField<FImpl>::setup(void)
{
    if (par().gammas == "all")
    {
        gamma_ = {
            Gamma::Algebra::Gamma5,
            Gamma::Algebra::Identity,
            Gamma::Algebra::GammaX,
            Gamma::Algebra::GammaY,
            Gamma::Algebra::GammaZ,
            Gamma::Algebra::GammaT,
            Gamma::Algebra::GammaXGamma5,
            Gamma::Algebra::GammaYGamma5,
            Gamma::Algebra::GammaZGamma5,
            Gamma::Algebra::GammaTGamma5,
            Gamma::Algebra::SigmaXY,
            Gamma::Algebra::SigmaXZ,
            Gamma::Algebra::SigmaXT,
            Gamma::Algebra::SigmaYZ,
            Gamma::Algebra::SigmaYT,
            Gamma::Algebra::SigmaZT
        };
    }
    else
    {
        gamma_ = strToVec<Gamma::Algebra>(par().gammas);
    }

    auto parse_momenta = [](const std::vector<std::string> &momenta, int dim,
          const std::string &desc)
    {
        std::vector<std::vector<int>> res;
        for (auto &pstr: momenta)
        {
            auto p = strToVec<int>(pstr);
            if (p.size() != dim)
            {
                HADRONS_ERROR(Size, desc + " has " + std::to_string(p.size())
                                    + " components instead of "
                                    + std::to_string(dim));
            }
            res.push_back(p);
        }
        return res;
    };
    smearMom_ = parse_momenta(par().smearMom, env().getNd()-1,
            "Smearing momentum");
    mom_ = parse_momenta(par().mom, env().getNd()-1, "Momentum");

    if(!hasDistributions_) {
        //distributions cache name
        distributionsCache_=getName()+"_distCache";
        //setup distributionsNames_
        distributionsNames_.clear();
        if (par().distributions.find('(') != std::string::npos)
        { //pairs
            distributionsNames_ = strToVec<stringPair>(par().distributions);
        }
        else
        { //no pairs
            auto dists = strToVec<std::string>(par().distributions);
            for (const auto &i: dists)
            {
                distributionsNames_.push_back(make_pair(i,i));
            }
        }
        //create mapping to unique ids (names may appear repeatedly)
        distributionsMap_.clear();
        {
            int id=0;
            for(const auto &i: distributionsNames_)
            {
                using std::get;
                auto res = distributionsMap_.insert(make_pair(get<0>(i), id));
                if(get<1>(res))
                {
                    id++;
                }
                res = distributionsMap_.insert(make_pair(get<1>(i), id));
                if(get<1>(res))
                {
                    id++;
                }
            }
        }
    }

    //allocate storage for the Fourier transforms of the distributions
    envCache(std::vector<ComplexField>, distributionsCache_,
            1, distributionsMap_.size(), envGetGrid(ComplexField));

    const auto smear_size=distributionsNames_.size()*smearMom_.size();
    envTmp(std::vector<ComplexField>, "smear_weight", 1, smear_size,
            envGetGrid(ComplexField));
    envTmp(Computation, "computation", 1, envGetGrid(FermionField),
            env().getNd() - 1, smear_size, gamma_.size(), par().block,
            par().cacheBlock, this);
    envTmp(FFT, "fft", 1, env().getGrid());

    auto &left_orig=envGet(std::vector<LatticeFermion>, par().left);
    auto &right_orig=envGet(std::vector<LatticeFermion>, par().right);
    const auto size_l=left_orig.size();
    const auto size_r=right_orig.size();
    envTmp(std::vector<FermionField>, "left", 1, size_l,
        envGetGrid(FermionField));
    envTmp(std::vector<FermionField>, "right", 1, size_r,
        envGetGrid(FermionField));
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2ASmearedMesonField<FImpl>::execute(void)
{
    //get copies of left and right
    envGetTmp(std::vector<FermionField>, left);
    envGetTmp(std::vector<FermionField>, right);
    {
        auto &orig_left=envGet(std::vector<FermionField>, par().left);
        auto &orig_right=envGet(std::vector<FermionField>, par().right);
        assert(left.size()==orig_left.size());
        assert(right.size()==orig_right.size());
        for(int i=0; i<left.size(); i++)
        {
            left[i]=orig_left[i];
        }
        for(int i=0; i<right.size(); i++)
        {
            right[i]=orig_right[i];
        }
    }

    auto &distrib=envGet(std::vector<ComplexField>, distributionsCache_);
    envGetTmp(FFT, fft);
    envGetTmp(std::vector<ComplexField>, smear_weight);

    int nt         = env().getDim().back();
    int N_i        = left.size();
    int N_j        = right.size();
    int ngamma     = gamma_.size();
    int nSmearMom  = smearMom_.size();
    int block      = par().block;
    int cacheBlock = par().cacheBlock;

    LOG(Message) << "Computing smeared all-to-all meson fields" << std::endl;
    LOG(Message) << "Left: '" << par().left << "' Right: '" << par().right
        << "'" << std::endl;
    LOG(Message) << "Smearing momenta:" << std::endl;
    for (auto &p: smearMom_)
    {
        LOG(Message) << "  " << p << std::endl;
    }
    LOG(Message) << "Momenta:" << std::endl;
    for (auto &p: mom_)
    {
        LOG(Message) << "  " << p << std::endl;
    }
    LOG(Message) << "Spin bilinears:" << std::endl;
    for (auto &g: gamma_)
    {
        LOG(Message) << "  " << g << std::endl;
    }
    LOG(Message) << "Meson field size: " << nt << "*" << N_i << "*" << N_j
                 << " (filesize "
                 << sizeString(nt*N_i*N_j*sizeof(HADRONS_A2AM_IO_TYPE))
                 << "/momentum/bilinear)" << std::endl;

    //---Fourier transform distributions---
    if(hasDistributions_)
    {
        LOG(Message) << "reusing Fourier transformed distributions"
                     << std::endl;
    }
    else
    {
        LOG(Message) << "Fourier transforming distributions" << std::endl;
        startTimer("Fourier transforming distributions");
        std::vector<int> mask(env().getNd(), 1);
        mask.back()=0; //transform only the spatial dimensions
        for(auto &i: distributionsMap_)
        {
            auto &src=envGet(ComplexField, std::get<0>(i));
            auto &dst=distrib.at(std::get<1>(i));
            fft.FFT_dim_mask(dst, src, mask, FFT::forward);
        }
        stopTimer("Fourier transforming distributions");
        hasDistributions_=true;
    }

    //---Fourier transform A2A vectors---
    LOG(Message) << "Fourier transforming A2A vectors" << std::endl;
    {
        startTimer("Fourier transform A2A vectors");
        std::vector<int> mask(env().getNd(), 1);
        mask.back()=0; //transform only the spatial dimensions
        for(auto &i: left)
        {
            fft.FFT_dim_mask(i, i, mask, FFT::forward);
        }
        for(auto &i: right)
        {
            fft.FFT_dim_mask(i, i, mask, FFT::forward);
        }
        stopTimer("Fourier transform A2A vectors");
    }

    for(const auto &m: mom_)
    {

        auto ionameFn = [this,&m,nSmearMom](const unsigned int ms,
                const unsigned int g)
        {
            std::stringstream ss;
            auto &distNames = distributionsNames_.at(ms/nSmearMom);

            ss << gamma_[g] << "_smear";
            for(auto pmu: smearMom_.at(ms%nSmearMom))
            {
                ss << '_' << pmu;
            }
            ss << "_mom";
            for(auto pmu: m)
            {
                ss << '_' << pmu;
            }
            ss << "_" << distNames.first << "_" << distNames.second;

            return ss.str();
        };

        auto filenameFn = [this, &ionameFn](const unsigned int ms,
                const unsigned int g)
        {
            return par().output + "." + std::to_string(vm().getTrajectory())
                   + "/" + ionameFn(ms, g) + ".h5";
        };

        auto metadataFn = [this,&m,nSmearMom](const unsigned int ms,
              const unsigned int g)
        {
            A2ASmearedMesonFieldMetadata md;

            for (auto pmu: smearMom_[ms%nSmearMom])
            {
                md.smearMom.push_back(pmu);
            }
            for (auto pmu: m)
            {
                md.mom.push_back(pmu);
            }
            md.gamma = gamma_[g];

            return md;
        };

        startTimer("Cshift momentum A2A vector");
        for(auto &field: right)
        {
            multidim_Cshift_inplace(field, m);
        }
        stopTimer("Cshift momentum A2A vector");

        startTimer("compute smearing weight");
        smearing_weight(smear_weight, distributionsNames_, distributionsMap_,
                distrib, smearMom_, m);
        stopTimer("compute smearing weight");
        Kernel kernel(gamma_, smear_weight, envGetGrid(FermionField));

        envGetTmp(Computation, computation);
        computation.execute(left, right, kernel,
                ionameFn, filenameFn, metadataFn);

        //restore A2A vectors
        {
            startTimer("Cshift momentum A2A vector");
            auto mNeg(m);
            for(auto &i: mNeg)
            {
                i*=-1;
            }
            for(auto &field: right)
            {
                multidim_Cshift_inplace(field, mNeg);
            }
            stopTimer("Cshift momentum A2A vector");
        }
    }
}

template<typename FImpl>
template<typename Lattice, typename Iterable>
void TA2ASmearedMesonField<FImpl>::multidim_Cshift_inplace(Lattice &lat,
        const Iterable &shifts)
{
    int dim=0;
    for(auto shift: shifts)
    {
        if(shift!=0)
        {
            lat=Cshift(lat, dim, shift);
        }
        dim++;
    }
}

//compute the smearing weight
//out[n*Nm+m](q)=dist_left[n](q)*dist_right[n](q+p-smearMom[m])*invSpVol
//where
//Nm=smearMom.size()
//q: coordinate of the (momentum space) field
//p=mom
//invSpVol=1/spatialVolume: factor needed to normalize the Fourier transform
template<typename FImpl>
void TA2ASmearedMesonField<FImpl>::smearing_weight(
        std::vector<ComplexField> &out,
        const std::vector<stringPair> &distributionsNames_,
        const std::map<std::string, int> &distributionsMap_,
        const std::vector<ComplexField> &distrib,
        const std::vector<std::vector<int>> &smearMom,
        const std::vector<int> &mom)
{
    using std::get;
    assert(out.size() == distributionsNames_.size()*smearMom.size());
    const int numSpatialDims=env().getNd()-1;
    const int timeDim=env().getNd()-1;
    const Complex invSpVol(env().getDim(timeDim) / env().getVolume());

    std::vector<int> totalShift(numSpatialDims);
    auto out_it = out.begin();
    for(auto &dist: distributionsNames_)
    {
        auto getDistribution=[&](const std::string &name) -> const ComplexField&
            {
                auto tmp_it = distributionsMap_.find(name);
                assert(tmp_it != distributionsMap_.end());
                return distrib.at(get<1>(*tmp_it));
            };
        auto &dist_left  = getDistribution(get<0>(dist));
        auto &dist_right = getDistribution(get<1>(dist));
        for(const auto &msmear: smearMom)
        {
            assert(msmear.size() == numSpatialDims);
            assert(mom.size() == numSpatialDims);
            for(int i=0; i<numSpatialDims; i++)
            {
                totalShift[i]=mom[i]-msmear[i];
            }
            *out_it=dist_right;
            multidim_Cshift_inplace(*out_it, totalShift);
            *out_it=*out_it * dist_left * invSpVol;
            out_it++;
        }
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_A2ASmearedMesonField_hpp_
