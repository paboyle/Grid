/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/action/gauge/Gauge_aggregate.h

Copyright (C) 2015

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

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
/*  END LEGAL */
#ifndef GRID_QCD_GAUGE_H
#define GRID_QCD_GAUGE_H

#include <Grid/qcd/action/gauge/GaugeImplementations.h>
#include <Grid/qcd/utils/WilsonLoops.h>
#include <Grid/qcd/action/gauge/WilsonGaugeAction.h>
#include <Grid/qcd/action/gauge/PlaqPlusRectangleAction.h>

namespace Grid {
namespace QCD {

typedef WilsonGaugeAction<PeriodicGimplR>          WilsonGaugeActionR;
typedef WilsonGaugeAction<PeriodicGimplF>          WilsonGaugeActionF;
typedef WilsonGaugeAction<PeriodicGimplD>          WilsonGaugeActionD;
typedef PlaqPlusRectangleAction<PeriodicGimplR>    PlaqPlusRectangleActionR;
typedef PlaqPlusRectangleAction<PeriodicGimplF>    PlaqPlusRectangleActionF;
typedef PlaqPlusRectangleAction<PeriodicGimplD>    PlaqPlusRectangleActionD;
typedef IwasakiGaugeAction<PeriodicGimplR>         IwasakiGaugeActionR;
typedef IwasakiGaugeAction<PeriodicGimplF>         IwasakiGaugeActionF;
typedef IwasakiGaugeAction<PeriodicGimplD>         IwasakiGaugeActionD;
typedef SymanzikGaugeAction<PeriodicGimplR>        SymanzikGaugeActionR;
typedef SymanzikGaugeAction<PeriodicGimplF>        SymanzikGaugeActionF;
typedef SymanzikGaugeAction<PeriodicGimplD>        SymanzikGaugeActionD;


typedef WilsonGaugeAction<ConjugateGimplR>          ConjugateWilsonGaugeActionR;
typedef WilsonGaugeAction<ConjugateGimplF>          ConjugateWilsonGaugeActionF;
typedef WilsonGaugeAction<ConjugateGimplD>          ConjugateWilsonGaugeActionD;
typedef PlaqPlusRectangleAction<ConjugateGimplR>    ConjugatePlaqPlusRectangleActionR;
typedef PlaqPlusRectangleAction<ConjugateGimplF>    ConjugatePlaqPlusRectangleActionF;
typedef PlaqPlusRectangleAction<ConjugateGimplD>    ConjugatePlaqPlusRectangleActionD;
typedef IwasakiGaugeAction<ConjugateGimplR>         ConjugateIwasakiGaugeActionR;
typedef IwasakiGaugeAction<ConjugateGimplF>         ConjugateIwasakiGaugeActionF;
typedef IwasakiGaugeAction<ConjugateGimplD>         ConjugateIwasakiGaugeActionD;
typedef SymanzikGaugeAction<ConjugateGimplR>        ConjugateSymanzikGaugeActionR;
typedef SymanzikGaugeAction<ConjugateGimplF>        ConjugateSymanzikGaugeActionF;
typedef SymanzikGaugeAction<ConjugateGimplD>        ConjugateSymanzikGaugeActionD;

}}


#endif
