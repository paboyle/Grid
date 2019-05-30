/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MFermion/EMLepton.hpp

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


#ifndef Hadrons_MFermion_EMLepton_hpp_
#define Hadrons_MFermion_EMLepton_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/*******************************************************************************
*
* Calculates a free lepton propagator with a sequential insertion of 
* i*\gamma_mu A_mu with a photon field A_mu 
*
*	L(x) = \sum_y S(x,y) i*\gamma_mu*A_mu S(y,xl) \delta_{(tl-x0),dt} 
*
* with a wall source for the lepton at tl
*
* In addition outputs the propagator without photon vertex
* 
*	L^{free}(x) =  S(x,xl) \delta_{(tl-x0),dt}
*
*
* options:
*  - action: fermion action used for propagator (string)
*  - emField: photon field A_mu (string)
*  - mass: input mass for the lepton propagator
*  - twist: twisted boundary for lepton propagator, e.g. "0.0 0.0 0.0 0.5"
*  - deltat: source-sink separation
*
*******************************************************************************/


/******************************************************************************
 *                         EMLepton                                           *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MFermion)

class EMLeptonPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(EMLeptonPar,
				    std::string,  action,
				    std::string, emField,
				    double, mass,
                                    std::string , boundary,
				    std::string,  twist,
                                    unsigned int, deltat);
};

template <typename FImpl>
class TEMLepton: public Module<EMLeptonPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    typedef PhotonR::GaugeField     EmField;
public:
    // constructor
    TEMLepton(const std::string name);
    // destructor
    virtual ~TEMLepton(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    unsigned int Ls_;
};

MODULE_REGISTER_TMP(EMLepton, TEMLepton<FIMPL>, MFermion);

/******************************************************************************
 *                 TEMLepton implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TEMLepton<FImpl>::TEMLepton(const std::string name)
: Module<EMLeptonPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TEMLepton<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().action, par().emField};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TEMLepton<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName(), getName() + "_free"};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TEMLepton<FImpl>::setup(void)
{
    Ls_ = env().getObjectLs(par().action);
    envCreateLat(PropagatorField, getName());
    envCreateLat(PropagatorField, getName() + "_free");
    envTmpLat(FermionField, "source", Ls_);
    envTmpLat(FermionField, "sol", Ls_);
    envTmpLat(FermionField, "tmp");
    envTmpLat(PropagatorField, "sourcetmp");
    envTmpLat(PropagatorField, "proptmp");
    envTmpLat(PropagatorField, "freetmp");
    envTmp(Lattice<iScalar<vInteger>>, "tlat",1, envGetGrid(LatticeComplex));

}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TEMLepton<FImpl>::execute(void)
{
    LOG(Message) << "Computing free fermion propagator '" << getName() << "'"
                 << std::endl;
    
    auto        &mat = envGet(FMat, par().action);
    RealD mass = par().mass;
    Complex ci(0.0,1.0);

    PropagatorField &Aslashlep = envGet(PropagatorField, getName());
    PropagatorField &lep = envGet(PropagatorField, getName() + "_free");
    
    envGetTmp(FermionField, source);
    envGetTmp(FermionField, sol);
    envGetTmp(FermionField, tmp);
    LOG(Message) << "Calculating a lepton Propagator with sequential Aslash insertion with lepton mass " 
                 << mass << " using the action '" << par().action
                 << "' for fixed source-sink separation of " << par().deltat << std::endl;

    envGetTmp(Lattice<iScalar<vInteger>>, tlat);
    LatticeCoordinate(tlat, Tp);


    std::vector<double> twist = strToVec<double>(par().twist);
    if(twist.size() != Nd)
    {
	HADRONS_ERROR(Size, "number of twist angles does not match number of dimensions");
    }
    std::vector<Complex> boundary = strToVec<Complex>(par().boundary);
    if(boundary.size() != Nd)
    {
	HADRONS_ERROR(Size, "number of boundary conditions does not match number of dimensions");
    }

    auto &stoch_photon = envGet(EmField,  par().emField);
    unsigned int nt = env().getDim(Tp);

    envGetTmp(PropagatorField, proptmp);
    envGetTmp(PropagatorField, freetmp);
    envGetTmp(PropagatorField, sourcetmp);

    std::vector<int> position;
    SitePropagator   id;
    id  = 1.;

    unsigned int tl=0; 

    //wallsource at tl
    sourcetmp = 1.;
    sourcetmp = where((tlat == tl), sourcetmp, 0.*sourcetmp);

    //free propagator from pt source 
    for (unsigned int s = 0; s < Ns; ++s)
    {
        LOG(Message) << "Calculation for spin= " << s << std::endl;
	if (Ls_ == 1)
	{
	    PropToFerm<FImpl>(source, sourcetmp, s, 0);
	}
	else
	{
	    PropToFerm<FImpl>(tmp, sourcetmp, s, 0);
	    // 5D source if action is 5d
	    mat.ImportPhysicalFermionSource(tmp, source);
	}
        sol = zero;
	mat.FreePropagator(source,sol,mass,boundary,twist);
	if (Ls_ == 1)
	{
            FermToProp<FImpl>(freetmp, sol, s, 0);
	}
        // create 4D propagators from 5D one if necessary
        if (Ls_ > 1)
        {
            mat.ExportPhysicalFermionSolution(sol, tmp);
            FermToProp<FImpl>(freetmp, tmp, s, 0);
        }
    }

    for(tl=0;tl<nt;tl++){

	//shift free propagator to different source positions
	//account for possible anti-periodic boundary in time
	proptmp = Cshift(freetmp,Tp, -tl);
	proptmp = where( tlat < tl, boundary[Tp]*proptmp, proptmp);

        // free propagator for fixed source-sink separation 
	lep = where(tlat == (tl-par().deltat+nt)%nt, proptmp, lep);

        // i*A_mu*gamma_mu
        sourcetmp = zero;
        for(unsigned int mu=0;mu<=3;mu++)
        {
	    Gamma gmu(Gamma::gmu[mu]);
	    sourcetmp +=  ci * PeekIndex<LorentzIndex>(stoch_photon, mu) *  (gmu * proptmp );
        }

        proptmp = zero;

        //sequential propagator from i*Aslash*S
        LOG(Message) << "Sequential propagator for t= " << tl << std::endl;
        for (unsigned int s = 0; s < Ns; ++s)
        {
            LOG(Message) << "Calculation for spin= " << s << std::endl;
	    if (Ls_ == 1)
	    {
		PropToFerm<FImpl>(source, sourcetmp, s, 0);
	    }
	    else
	    {
		PropToFerm<FImpl>(tmp, sourcetmp, s, 0);
		// 5D source if action is 5d
		mat.ImportPhysicalFermionSource(tmp, source);
	    }
            sol = zero;
	    mat.FreePropagator(source,sol,mass,boundary,twist);
	    if (Ls_ == 1)
	    {
        	FermToProp<FImpl>(proptmp, sol, s, 0);
	    }
            // create 4D propagators from 5D one if necessary
            if (Ls_ > 1)
            {
                mat.ExportPhysicalFermionSolution(sol, tmp);
                FermToProp<FImpl>(proptmp, tmp, s, 0);
            }
	}
	// keep the result for the desired delta t
	Aslashlep = where(tlat == (tl-par().deltat+nt)%nt, proptmp, Aslashlep);
    }

    //account for possible anti-periodic boundary in time
    Aslashlep = where( tlat >= nt-par().deltat, boundary[Tp]*Aslashlep, Aslashlep);
    lep = where( tlat >= nt-par().deltat, boundary[Tp]*lep, lep);

}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MFermion_EMLepton_hpp_
