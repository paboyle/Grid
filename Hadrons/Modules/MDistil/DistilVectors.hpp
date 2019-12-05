/*************************************************************************************
 
 Grid physics library, www.github.com/paboyle/Grid
 
 Source file: Hadrons/Modules/MDistil/DistilVectors.hpp
 
 Copyright (C) 2019
 
 Author: Felix Erben <ferben@ed.ac.uk>
 Author: Michael Marshall <Michael.Marshall@ed.ac.uk>
 
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

#ifndef Hadrons_MDistil_DistilVectors_hpp_
#define Hadrons_MDistil_DistilVectors_hpp_

#include <Hadrons/Modules/MDistil/Distil.hpp>

BEGIN_HADRONS_NAMESPACE
BEGIN_MODULE_NAMESPACE(MDistil)

/******************************************************************************
 *                         DistilVectors                                      *
 *                (Create rho and/or phi vectors)                             *
 ******************************************************************************/

class DistilVectorsPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(DistilVectorsPar,
                                    std::string, noise,
                                    std::string, perambulator,
                                    std::string, lapevec,
                                    std::string, rho,
                                    std::string, phi,
                                    std::string, DistilParams);
};

template <typename FImpl>
class TDistilVectors: public Module<DistilVectorsPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    // constructor
    TDistilVectors(const std::string name);
    // destructor
    virtual ~TDistilVectors(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
protected:
    std::unique_ptr<GridCartesian> grid3d; // Owned by me, so I must delete it
public:
    // These variables contain parameters
    std::string RhoName;
    std::string PhiName;
};

MODULE_REGISTER_TMP(DistilVectors, TDistilVectors<FIMPL>, MDistil);

/******************************************************************************
 *                 TDistilVectors implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TDistilVectors<FImpl>::TDistilVectors(const std::string name) : Module<DistilVectorsPar>(name) {}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TDistilVectors<FImpl>::getInput(void)
{
    return {par().noise,par().perambulator,par().lapevec,par().DistilParams};
}

template <typename FImpl>
std::vector<std::string> TDistilVectors<FImpl>::getOutput(void)
{
    RhoName = par().rho;
    PhiName = par().phi;
    if (RhoName.empty() && PhiName.empty())
    {
        HADRONS_ERROR(Argument,"No output specified");
    }
    std::vector<std::string> out;
    if (!RhoName.empty())
        out.push_back(RhoName);
    if (!PhiName.empty())
        out.push_back(PhiName);
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDistilVectors<FImpl>::setup(void)
{
    // We expect the perambulator to have been created with these indices
    auto &perambulator = envGet(PerambTensor, par().perambulator);
    if (!perambulator.ValidateIndexNames())
    {
        HADRONS_ERROR(Range,"Perambulator index names bad");
    }

    const DistilParameters &dp{envGet(DistilParameters, par().DistilParams)};
    const int Nt{env().getDim(Tdir)};
    const bool full_tdil{ dp.TI == Nt };
    const int Nt_inv{ full_tdil ? 1 : dp.TI };
    
    if (!RhoName.empty())
        envCreate(std::vector<FermionField>, RhoName, 1, dp.nnoise*dp.LI*dp.SI*Nt_inv, envGetGrid(FermionField));
    if (!PhiName.empty())
        envCreate(std::vector<FermionField>, PhiName, 1, dp.nnoise*dp.LI*dp.SI*Nt_inv, envGetGrid(FermionField));
    
    Coordinate latt_size   = GridDefaultLatt();
    Coordinate mpi_layout  = GridDefaultMpi();
    Coordinate simd_layout_3 = GridDefaultSimd(Nd-1, vComplex::Nsimd());
    latt_size[Nd-1] = 1;
    simd_layout_3.push_back( 1 );
    mpi_layout[Nd-1] = 1;
    GridCartesian * const grid4d{env().getGrid()};
    MakeLowerDimGrid(grid3d, grid4d);
    
    envTmp(LatticeSpinColourVector, "source4d",1,LatticeSpinColourVector(grid4d));
    envTmp(LatticeSpinColourVector, "source3d",1,LatticeSpinColourVector(grid3d.get()));
    envTmp(LatticeColourVector, "source3d_nospin",1,LatticeColourVector(grid3d.get()));
    envTmp(LatticeSpinColourVector, "sink3d",1,LatticeSpinColourVector(grid3d.get()));
    envTmp(LatticeColourVector, "evec3d",1,LatticeColourVector(grid3d.get()));
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDistilVectors<FImpl>::execute(void)
{
    auto &noise        = envGet(NoiseTensor,  par().noise);
    auto &perambulator = envGet(PerambTensor, par().perambulator);
    auto &epack        = envGet(Grid::Hadrons::EigenPack<LatticeColourVector>, par().lapevec);
    const DistilParameters &dp{envGet(DistilParameters, par().DistilParams)};
    
    envGetTmp(LatticeSpinColourVector, source4d);
    envGetTmp(LatticeSpinColourVector, source3d);
    envGetTmp(LatticeColourVector,     source3d_nospin);
    envGetTmp(LatticeSpinColourVector, sink3d);
    envGetTmp(LatticeColourVector,     evec3d);
    
    GridCartesian * const grid4d{env().getGrid()};
    const int Ntlocal{ grid4d->LocalDimensions()[3] };
    const int Ntfirst{ grid4d->LocalStarts()[3] };
    
    const int Nt{env().getDim(Tdir)}; 
    const bool full_tdil{ dp.TI == Nt }; 
    const int Nt_inv{ full_tdil ? 1 : dp.TI };
    
    int vecindex;
    if (!RhoName.empty())
    {
        auto &rho = envGet(std::vector<FermionField>, RhoName);
        for (int inoise = 0; inoise < dp.nnoise; inoise++) 
	{
            for (int dk = 0; dk < dp.LI; dk++) 
	    {
                for (int dt = 0; dt < Nt_inv; dt++) 
		{
                    for (int ds = 0; ds < dp.SI; ds++) 
		    {
                        vecindex = inoise + dp.nnoise * (dk + dp.LI * (ds + dp.SI * dt));
                        rho[vecindex] = 0;
                        for (int it = dt; it < Nt; it += dp.TI)
			{
                            const int t_inv{full_tdil ? dp.tsrc : it};
                            if (t_inv >= Ntfirst && t_inv < Ntfirst + Ntlocal) 
			    {
                                for (int ik = dk; ik < dp.nvec; ik += dp.LI)
				{
                                    for (int is = ds; is < Ns; is += dp.SI)
				    {
                                        ExtractSliceLocal(evec3d,epack.evec[ik],0,t_inv-Ntfirst,Tdir);
                                        source3d_nospin = evec3d * noise.tensor(inoise, t_inv, ik, is);
                                        source3d=0;
                                        pokeSpin(source3d,source3d_nospin,is);
                                        source4d=0;
                                        InsertSliceLocal(source3d,source4d,0,t_inv-Ntfirst,Tdir);
                                        rho[vecindex] += source4d;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    if (!PhiName.empty()) 
    {
        auto &phi = envGet(std::vector<FermionField>, PhiName);
        for (int inoise = 0; inoise < dp.nnoise; inoise++) 
	{
            for (int dk = 0; dk < dp.LI; dk++) 
	    {
                for (int dt = 0; dt < Nt_inv; dt++) 
		{
                    for (int ds = 0; ds < dp.SI; ds++) 
		    {
                        vecindex = inoise + dp.nnoise * (dk + dp.LI * (ds + dp.SI * dt));
                        phi[vecindex] = 0;
                        for (int t = Ntfirst; t < Ntfirst + Ntlocal; t++) 
			{
                            sink3d=0;
                            for (int ivec = 0; ivec < dp.nvec; ivec++) 
			    {
                                ExtractSliceLocal(evec3d,epack.evec[ivec],0,t-Ntfirst,Tdir);
                                sink3d += evec3d * perambulator.tensor(t, ivec, dk, inoise,dt,ds);
                            }
                            InsertSliceLocal(sink3d,phi[vecindex],0,t-Ntfirst,Tdir);
                        }
                    }
                }
            }
        }
    }
}

END_MODULE_NAMESPACE
END_HADRONS_NAMESPACE
#endif // Hadrons_MDistil_DistilVectors_hpp_
