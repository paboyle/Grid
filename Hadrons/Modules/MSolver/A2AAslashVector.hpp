/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MSolver/A2AAslashVector.hpp

Copyright (C) 2015-2018

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
#ifndef Hadrons_MSolver_A2AAslashVector_hpp_
#define Hadrons_MSolver_A2AAslashVector_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Solver.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                       Create all-to-all V & W vectors                      *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSolver)

/****************************************************************************
*  Calculate a sequential propagator on an insertion of i*g_mu*A_mu 
*  on an A2A vector
*
*  vv_i(y) = S(y,x) * i * g_mu*A_mu(x) * v_i(x)
*
*  with
*
*  - vector: A2A vector v_i(x)
*  - emField: A_mu(x): electromagnetic photon field
*  - solver: the solver for calculating the sequential propagator
*
*****************************************************************************/


class A2AAslashVectorPar: Serializable
{
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(A2AAslashVectorPar,
                                  std::string, vector,
                                  std::string, emField,
                                  std::string, solver);
};

template <typename FImpl>
class TA2AAslashVector : public Module<A2AAslashVectorPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    SOLVER_TYPE_ALIASES(FImpl,);
public:
    typedef PhotonR::GaugeField     EmField;
public:
    // constructor
    TA2AAslashVector(const std::string name);
    // destructor
    virtual ~TA2AAslashVector(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    unsigned int Ls_;
};

MODULE_REGISTER_TMP(A2AAslashVector,TA2AAslashVector<FIMPL>, MSolver);
MODULE_REGISTER_TMP(ZA2AAslashVector,TA2AAslashVector<ZFIMPL>, MSolver);

/******************************************************************************
 *                       TA2AAslashVector implementation                       *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TA2AAslashVector<FImpl>::TA2AAslashVector(const std::string name)
: Module<A2AAslashVectorPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TA2AAslashVector<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().vector, par().emField, par().solver};

    return in;
}

template <typename FImpl>
std::vector<std::string> TA2AAslashVector<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2AAslashVector<FImpl>::setup(void)
{
    Ls_  = env().getObjectLs(par().solver);
    auto &vvector = envGet(std::vector<FermionField>, par().vector);
    unsigned int Nmodes = vvector.size();
    envCreate(std::vector<FermionField>, getName(), 1, 
              Nmodes, envGetGrid(FermionField));
   
    envTmpLat(FermionField, "v4dtmp");
    envTmpLat(FermionField, "v5dtmp", Ls_);
    envTmpLat(FermionField, "v5dtmp_sol", Ls_);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2AAslashVector<FImpl>::execute(void)
{
    auto &solver = envGet(Solver, par().solver);
    auto &stoch_photon = envGet(EmField,  par().emField);
    auto &vvector = envGet(std::vector<FermionField>, par().vector);
    auto &Aslashv = envGet(std::vector<FermionField>, getName());
    unsigned int Nmodes = vvector.size();
    auto &mat = solver.getFMat();
    envGetTmp(FermionField, v4dtmp);
    envGetTmp(FermionField, v5dtmp);
    envGetTmp(FermionField, v5dtmp_sol);

    Complex ci(0.0,1.0);


    startTimer("Seq Aslash");

    LOG(Message) << "Calculate Sequential propagator on Aslash * v with the A2A vector " << par().vector
                  << " and the photon field " << par().emField << std::endl;


    for(unsigned int i=0; i<Nmodes; i++)
    {
	v4dtmp = zero;
	startTimer("Multiply Aslash");
	for(unsigned int mu=0;mu<=3;mu++)
    	{
	    Gamma gmu(Gamma::gmu[mu]);
		 v4dtmp +=  ci * PeekIndex<LorentzIndex>(stoch_photon, mu) * (gmu * vvector[i]);
	}
	stopTimer("Multiply Aslash");

	if (Ls_ == 1)
	{
	    solver(Aslashv[i], v4dtmp);
	}
	else
	{
	    mat.ImportPhysicalFermionSource(v4dtmp, v5dtmp);
	    solver(v5dtmp_sol, v5dtmp);
	    mat.ExportPhysicalFermionSolution(v5dtmp_sol, v4dtmp);
	    Aslashv[i] = v4dtmp;
	}
    }

    stopTimer("Seq Aslash");
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSolver_A2AAslashVector_hpp_
