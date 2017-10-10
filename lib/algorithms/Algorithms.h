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

#include <Grid/algorithms/SparseMatrix.h>
#include <Grid/algorithms/LinearOperator.h>
#include <Grid/algorithms/Preconditioner.h>

#include <Grid/algorithms/approx/Zolotarev.h>
#include <Grid/algorithms/approx/Chebyshev.h>
#include <Grid/algorithms/approx/Remez.h>
#include <Grid/algorithms/approx/MultiShiftFunction.h>
#include <Grid/algorithms/approx/Forecast.h>

#include <Grid/algorithms/densematrix/DenseMatrix.h>
#include <Grid/algorithms/densematrix/Francis.h>
#include <Grid/algorithms/densematrix/Householder.h>

#include <Grid/algorithms/iterative/ConjugateGradient.h>
#include <Grid/algorithms/iterative/ConjugateResidual.h>
#include <Grid/algorithms/iterative/NormalEquations.h>
#include <Grid/algorithms/iterative/SchurRedBlack.h>
#include <Grid/algorithms/iterative/ConjugateGradientMultiShift.h>
#include <Grid/algorithms/iterative/ConjugateGradientMixedPrec.h>
#include <Grid/algorithms/iterative/BlockConjugateGradient.h>
#include <Grid/algorithms/iterative/ConjugateGradientReliableUpdate.h>
#include <Grid/algorithms/iterative/ImplicitlyRestartedLanczos.h>
#include <Grid/algorithms/iterative/ImplicitlyRestartedLanczosCJ.h>
#include <Grid/algorithms/iterative/SimpleLanczos.h>
#include <Grid/algorithms/CoarsenedMatrix.h>
#include <Grid/algorithms/FFT.h>

// EigCg
// Pcg
// Hdcg
// GCR
// etc..

#endif
