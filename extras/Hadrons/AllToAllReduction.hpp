#ifndef A2A_Reduction_hpp_
#define A2A_Reduction_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Environment.hpp>
#include <Grid/Hadrons/Solver.hpp>

BEGIN_HADRONS_NAMESPACE

////////////////////////////////////////////
// A2A Meson Field Inner Product
////////////////////////////////////////////

template <class FermionField>
void sliceInnerProductMesonField(std::vector<std::vector<ComplexD>> &mat,
                                 const std::vector<Lattice<FermionField>> &lhs,
                                 const std::vector<Lattice<FermionField>> &rhs,
                                 int orthogdim)
{
    typedef typename FermionField::scalar_type scalar_type;
    typedef typename FermionField::vector_type vector_type;

    int Lblock = lhs.size();
    int Rblock = rhs.size();

    GridBase *grid = lhs[0]._grid;

    const int Nd = grid->_ndimension;
    const int Nsimd = grid->Nsimd();
    int Nt = grid->GlobalDimensions()[orthogdim];

    assert(mat.size() == Lblock * Rblock);
    for (int t = 0; t < mat.size(); t++)
    {
        assert(mat[t].size() == Nt);
    }

    int fd = grid->_fdimensions[orthogdim];
    int ld = grid->_ldimensions[orthogdim];
    int rd = grid->_rdimensions[orthogdim];

    // will locally sum vectors first
    // sum across these down to scalars
    // splitting the SIMD
    std::vector<vector_type, alignedAllocator<vector_type>> lvSum(rd * Lblock * Rblock);
    for(int r=0;r<rd * Lblock * Rblock;r++)
    {
        lvSum[r]=zero;
    }
    std::vector<scalar_type> lsSum(ld * Lblock * Rblock, scalar_type(0.0));

    int e1 = grid->_slice_nblock[orthogdim];
    int e2 = grid->_slice_block[orthogdim];
    int stride = grid->_slice_stride[orthogdim];

    // std::cout << GridLogMessage << " Entering first parallel loop " << std::endl;
    // Parallelise over t-direction doesn't expose as much parallelism as needed for KNL
    parallel_for(int r = 0; r < rd; r++)
    {
        int so = r * grid->_ostride[orthogdim]; // base offset for start of plane
        for (int n = 0; n < e1; n++)
        {
            for (int b = 0; b < e2; b++)
            {
                int ss = so + n * stride + b;
                for (int i = 0; i < Lblock; i++)
                {
                    auto left = conjugate(lhs[i]._odata[ss]);
                    for (int j = 0; j < Rblock; j++)
                    {
                        int idx = i + Lblock * j + Lblock * Rblock * r;
                        auto right = rhs[j]._odata[ss];
                        vector_type vv = left()(0)(0) * right()(0)(0) 
                                       + left()(0)(1) * right()(0)(1) 
                                       + left()(0)(2) * right()(0)(2) 
                                       + left()(1)(0) * right()(1)(0) 
                                       + left()(1)(1) * right()(1)(1) 
                                       + left()(1)(2) * right()(1)(2) 
                                       + left()(2)(0) * right()(2)(0) 
                                       + left()(2)(1) * right()(2)(1) 
                                       + left()(2)(2) * right()(2)(2) 
                                       + left()(3)(0) * right()(3)(0) 
                                       + left()(3)(1) * right()(3)(1) 
                                       + left()(3)(2) * right()(3)(2);

                        lvSum[idx] = lvSum[idx] + vv;
                    }
                }
            }
        }
    }

    // std::cout << GridLogMessage << " Entering second parallel loop " << std::endl;
    // Sum across simd lanes in the plane, breaking out orthog dir.
    parallel_for(int rt = 0; rt < rd; rt++)
    {
        std::vector<int> icoor(Nd);
        for (int i = 0; i < Lblock; i++)
        {
            for (int j = 0; j < Rblock; j++)
            {
                iScalar<vector_type> temp;
                std::vector<iScalar<scalar_type>> extracted(Nsimd);
                temp._internal = lvSum[i + Lblock * j + Lblock * Rblock * rt];
                extract(temp, extracted);
                for (int idx = 0; idx < Nsimd; idx++)
                {
                    grid->iCoorFromIindex(icoor, idx);
                    int ldx = rt + icoor[orthogdim] * rd;
                    int ij_dx = i + Lblock * j + Lblock * Rblock * ldx;
                    lsSum[ij_dx] = lsSum[ij_dx] + extracted[idx]._internal;
                }
            }
        }
    }

    // std::cout << GridLogMessage << " Entering non parallel loop " << std::endl;
    for (int t = 0; t < fd; t++)
    {
        int pt = t/ld; // processor plane
        int lt = t%ld;
        for (int i = 0; i < Lblock; i++)
        {
            for (int j = 0; j < Rblock; j++)
            {
                if (pt == grid->_processor_coor[orthogdim])
                {
                    int ij_dx = i + Lblock * j + Lblock * Rblock * lt;
                    mat[i + j * Lblock][t] = lsSum[ij_dx];
                }
                else
                {
                    mat[i + j * Lblock][t] = scalar_type(0.0);
                }
                
            }
        }
    }
    // std::cout << GridLogMessage << " Done " << std::endl;
    // defer sum over nodes.
    return;
}

END_HADRONS_NAMESPACE

#endif // A2A_Reduction_hpp_