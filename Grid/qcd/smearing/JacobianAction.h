/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/action/gauge/JacobianAction.h

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

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
			   /*  END LEGAL */
#pragma once

NAMESPACE_BEGIN(Grid);

////////////////////////////////////////////////////////////////////////
// Jacobian Action .. 
////////////////////////////////////////////////////////////////////////
template <class Gimpl>
class JacobianAction : public Action<typename Gimpl::GaugeField> {
public:  
  INHERIT_GIMPL_TYPES(Gimpl);

  SmearedConfigurationMasked<Gimpl> * smearer;
  /////////////////////////// constructors
  explicit JacobianAction(SmearedConfigurationMasked<Gimpl> * _smearer ) { smearer=_smearer;};

  virtual std::string action_name() {return "JacobianAction";}

  virtual std::string LogParameters(){
    std::stringstream sstream;
    sstream << GridLogMessage << "[JacobianAction] " << std::endl;
    return sstream.str();
  }

  //////////////////////////////////
  // Usual cases are not used
  //////////////////////////////////
  virtual void refresh(const GaugeField &U, GridSerialRNG &sRNG, GridParallelRNG &pRNG){ assert(0);};
  virtual RealD S(const GaugeField &U) { assert(0); }
  virtual void deriv(const GaugeField &U, GaugeField &dSdU) { assert(0);  }

  //////////////////////////////////
  // Functions of smart configs only
  //////////////////////////////////
  virtual void refresh(ConfigurationBase<GaugeField> & U, GridSerialRNG &sRNG, GridParallelRNG& pRNG)
  {
    return;
  }
  virtual RealD S(ConfigurationBase<GaugeField>& U)
  {
    // det M = e^{ - ( - logDetM) }
    assert( &U == smearer );
    return -smearer->logDetJacobian();
  }
  virtual RealD Sinitial(ConfigurationBase<GaugeField>& U) 
  {
    return S(U);
  }
  virtual void deriv(ConfigurationBase<GaugeField>& U, GaugeField& dSdU)
  {
    assert( &U == smearer );
    smearer->logDetJacobianForce(dSdU);
  }

private:
 };

NAMESPACE_END(Grid);

