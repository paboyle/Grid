    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: Grid/algorithms/multigrid/MultiGrid.h

    Copyright (C) 2023

Author: Peter Boyle <pboyle@bnl.gov>

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
#pragma once

#include <Grid/algorithms/multigrid/Aggregates.h>
#include <Grid/algorithms/multigrid/Geometry.h>
#include <Grid/algorithms/multigrid/CoarsenedMatrix.h>
#include <Grid/algorithms/multigrid/GeneralCoarsenedMatrixMultiRHSV2.h>
// DEPRECATED: the V1 coarse operators and the Aggregation-based single-RHS
// ADEF-2.  Nothing in the library uses them; the PVdagM and HDCG chains are
// on V2 (PVdagMMultiGrid.h, HDCGMultiGrid.h).  Kept so the pre-2026 drivers
// in tests/debug and examples still build; removing this block is the
// deletion gate for them.
#include <Grid/algorithms/multigrid/deprecated/GeneralCoarsenedMatrix.h>
#include <Grid/algorithms/multigrid/deprecated/GeneralCoarsenedMatrixMultiRHS.h>
#include <Grid/algorithms/multigrid/deprecated/TwoLevelADEF2.h>
#include <Grid/algorithms/multigrid/MrhsPromotedOperator.h>
#include <Grid/algorithms/multigrid/Smoothers.h>
#include <Grid/algorithms/multigrid/PVdagMMultiGridParams.h>
// PVdagMOperators.h / MrhsMultiGrid.h / PVdagMMultiGrid.h /
// DenseCoarseMatrix.h / MultiGridIO.h / HDCGMultiGrid.h are NOT in this
// umbrella: consumers include PVdagMMultiGrid.h or HDCGMultiGrid.h
// explicitly (they pull the dense stack, BLAS, and the scidac I/O, which
// is declared after Algorithms.h in Grid.h).  Keeping MrhsMultiGrid.h out
// of the Algorithms.h chain also keeps its class names away from the
// pre-2026 drivers that define their own copies.
