/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: lib/algorithms/iterative/GeneralisedMinimalResidual.h

Copyright (C) 2015
Copyright (C) 2016


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
#ifndef GRID_GENERALISED_MINIMAL_RESIDUAL_H
#define GRID_GENERALISED_MINIMAL_RESIDUAL_H

///////////////////////////////////////////////////////////////////////////////////////////////////////
// from Y. Saad - Iterative Methods for Sparse Linear Systems, PP 172
// Compute r0 = b − Ax0 , β := ||r0||2 , and v1 := r0 /β
// For j = 1, 2, ..., m Do:
//   Compute wj := Avj
//   For i = 1, ..., j Do:
//     hij := (wj , vi)
//     wj := wj − hij vi
//   EndDo
//   hj+1,j = ||wj||2 . If hj+1,j = 0 set m := j and go to HERE
//   vj+1 = wj /hj+1,j
// EndDo
// Define the (m + 1) × m Hessenberg matrix H̄m = {hij}1≤i≤m+1,1≤j≤m. [HERE]
// Compute ym the minimizer of ||βe1 − H̄m y||2 and xm = x0 + Vm ym.
///////////////////////////////////////////////////////////////////////////////////////////////////////

// want to solve Ax = b -> A = LinOp, psi = x, b = src

namespace Grid
{
template< class Field >
class GeneralisedMinimalResidual : public OperatorFunction< Field >
{
public:
    bool ErrorOnNoConverge; // Throw an assert when GMRES fails to converge,
                            // defaults to True.
    RealD   Tolerance;
    Integer MaxIterations;
    Integer IterationsToComplete; // Number of iterations the GMRES took to
                                  // finish. Filled in upon completion

    GeneralisedMinimalResidual( RealD   tol,
                                Integer maxit,
                                bool    err_on_no_conv = true )
        : Tolerance( tol )
        , MaxIterations( maxit )
        , ErrorOnNoConverge( err_on_no_conv ){};

    // want to solve Ax = b -> A = LinOp, psi = x, b = src

    void operator()( LinearOperatorBase< Field > &LinOp,
                     const Field &                src,
                     Field &                      psi )
    {
        std::cout << GridLogMessage
                  << "GeneralisedMinimalResidual: Start of operator()"
                  << std::endl;
        psi.checkerboard = src.checkerboard;
        conformable( psi, src );

        Field r( src );
        Field mmv( src );

        std::vector< Field > v( MaxIterations + 1, src );

        RealD beta{};
        RealD b{};
        RealD d{};

        Eigen::MatrixXcd H
            = Eigen::MatrixXcd::Zero( MaxIterations + 1, MaxIterations );

        // Compute r0 = b − Ax0 , β := ||r0||2 , and v1 := r0 /β
        LinOp.Op( psi, mmv );

        r      = src - mmv;
        beta   = norm2( r );
        V[ 0 ] = ( 1 / beta ) * r;

        for( auto j = 0; j < MaxIterations; ++j )
        {
            LinOp.Op( V[ j ], mmv );

            for( auto i = 0; i < j; ++i )
            {
                std::cout
                    << GridLogMessage
                    << "GeneralisedMinimalResidual: End of inner iteration "
                    << i << std::endl;
                H( i, j ) = innerProduct( mmv, v[ i ] );
                mmv = mmv - H( i, j ) * V[ i ];
            }

            H( j + 1, j ) = norm2( mmv );

            std::cout << GridLogMessage << "GeneralisedMinimalResidual: H"
                      << j + 1 << "," << j << "= " << H( j + 1, j )
                      << std::endl;
            if( H( j + 1, j ) == 0. )
            {
                IterationsToComplete = j;
                break;
            }

            V[ j + 1 ] = ( 1. / H( j + 1, j ) ) * mmv;
            std::cout << GridLogMessage
                      << "GeneralisedMinimalResidual: End of outer iteration "
                      << j << std::endl;
        }
        std::cout << GridLogMessage
                  << "GeneralisedMinimalResidual: End of operator()"
                  << std::endl;
    }
};
}
#endif

// Note: The DD-αAMG codebase turns around the Hessenberg matrix

void arnoldiStep()
{
    w = D * V[ j ];

    for( auto i = 0; i <= j; ++i )
        H( i, j ) = innerProduct( V[ j + 1 ], w );
    w = w - H( i, j ) * V[ i ];

    H( j + 1, j ) = norm2( w );

    V[ j + 1 ] = w / H( j + 1, j );
}

void qr_update_PRECISION()
{
    // update QR factorization
    // apply previous Givens rotation
    for( auto i = 0; i < j; i++ )
    {
        beta = -s[ i ] * H( i, j ) + c[ i ] * H( i + 1, j );
        H( i, j ) = std::conj( c[ i ] ) * H( i, j )
                    + std::conj( s[ i ] ) * H( i + 1, j );
        H( i + 1, j ) = beta;
    }

    // compute current Givens rotation
    beta   = sqrt( std::norm( H( j, j ) ) + std::norm( H( j, j + 1 ) ) );
    s[ j ] = H( j + 1, j ) / beta;
    c[ j ] = H( j, j ) / beta;

    // update right column
    gamma[ j + 1 ] = -s[ j ] * gamma[ j ];
    gamma[ j ]     = std::conj( c[ j ] ) * gamma[ j ];

    // apply current Givens rotation
    H( j, j )     = beta;
    H( j + 1, j ) = 0;
}

// check
void compute_solution_PRECISION()
{
    for( auto i = j; i >= 0; i-- )
    {
        y[ i ] = gamma[ i ];
        for( auto k = i + 1; k <= j; k++ )
            y[ i ] -= H( i, k ) * y[ k ];
        y[ i ] /= H( i, i );
    }

    if( true ) // TODO ???
    {
        for( i = 0; i <= j; i++ )
            x = x + V[ i ] * y[ i ];
    }
    else
    {
        x = y[ 0 ] * V[ 0 ];
        for( i = 1; i <= j; i++ )
            x = x + V[ i ] * y[ i ];
    }
}
