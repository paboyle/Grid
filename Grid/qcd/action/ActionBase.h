/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/action/ActionBase.h

Copyright (C) 2015-2016

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: neo <cossu@post.kek.jp>
Author: Guido Cossu <guido.cossu@ed.ac.uk>

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

#ifndef ACTION_BASE_H
#define ACTION_BASE_H

NAMESPACE_BEGIN(Grid);

template <class GaugeField >
class Action 
{

public:
  bool is_smeared = false;
  RealD deriv_norm_sum;
  RealD deriv_max_sum;
  RealD Fdt_norm_sum;
  RealD Fdt_max_sum;
  int   deriv_num;
  RealD deriv_us;
  RealD S_us;
  RealD refresh_us;
  void  reset_timer(void)        {
    deriv_us = S_us = refresh_us = 0.0;
    deriv_norm_sum = deriv_max_sum=0.0;
    Fdt_max_sum =  Fdt_norm_sum = 0.0;
    deriv_num=0;
  }
  void  deriv_log(RealD nrm, RealD max,RealD Fdt_nrm,RealD Fdt_max) {
    if ( max > deriv_max_sum ) {
      deriv_max_sum=max;
    }
    deriv_norm_sum+=nrm;
    if ( Fdt_max > Fdt_max_sum ) {
      Fdt_max_sum=Fdt_max;
    }
    Fdt_norm_sum+=Fdt_nrm; deriv_num++;
  }
  RealD deriv_max_average(void)       { return deriv_max_sum; };
  RealD deriv_norm_average(void)      { return deriv_norm_sum/deriv_num; };
  RealD Fdt_max_average(void)         { return Fdt_max_sum; };
  RealD Fdt_norm_average(void)        { return Fdt_norm_sum/deriv_num; };
  RealD deriv_timer(void)        { return deriv_us; };
  RealD S_timer(void)            { return S_us; };
  RealD refresh_timer(void)      { return refresh_us; };
  void deriv_timer_start(void)   { deriv_us-=usecond(); }
  void deriv_timer_stop(void)    { deriv_us+=usecond(); }
  void refresh_timer_start(void) { refresh_us-=usecond(); }
  void refresh_timer_stop(void)  { refresh_us+=usecond(); }
  void S_timer_start(void)       { S_us-=usecond(); }
  void S_timer_stop(void)        { S_us+=usecond(); }
  // Heatbath?
  virtual void refresh(const GaugeField& U, GridSerialRNG &sRNG, GridParallelRNG& pRNG) = 0; // refresh pseudofermions
  virtual RealD S(const GaugeField& U) = 0;                             // evaluate the action
  virtual RealD Sinitial(const GaugeField& U) { return this->S(U); } ;  // if the refresh computes the action, can cache it. Alternately refreshAndAction() ?
  virtual void deriv(const GaugeField& U, GaugeField& dSdU) = 0;        // evaluate the action derivative
  virtual std::string action_name()    = 0;                             // return the action name
  virtual std::string LogParameters()  = 0;                             // prints action parameters
  virtual ~Action(){}
};

NAMESPACE_END(Grid);

#endif // ACTION_BASE_H
