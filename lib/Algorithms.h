    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/Algorithms.h

    Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
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
#ifndef GRID_ALGORITHMS_H
#define GRID_ALGORITHMS_H

#include <algorithms/SparseMatrix.h>
#include <algorithms/LinearOperator.h>
#include <algorithms/Preconditioner.h>

#include <algorithms/approx/Zolotarev.h>
#include <algorithms/approx/Chebyshev.h>
#include <algorithms/approx/Remez.h>
#include <algorithms/approx/MultiShiftFunction.h>

#include <algorithms/iterative/ConjugateGradient.h>
#include <algorithms/iterative/ConjugateResidual.h>
#include <algorithms/iterative/NormalEquations.h>
#include <algorithms/iterative/SchurRedBlack.h>

#include <algorithms/iterative/ConjugateGradientMultiShift.h>

// Lanczos support
#include <algorithms/iterative/MatrixUtils.h>
#include <algorithms/iterative/ImplicitlyRestartedLanczos.h>

#include <algorithms/CoarsenedMatrix.h>

// Eigen/lanczos
// EigCg
// MCR
// Pcg
// Multishift CG
// Hdcg
// GCR
// etc..

// integrator/Leapfrog
// integrator/Omelyan
// integrator/ForceGradient

// montecarlo/hmc
// montecarlo/rhmc
// montecarlo/metropolis
// etc...


#endif
