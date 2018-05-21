/*******************************************************************************
 Grid physics library, www.github.com/paboyle/Grid

 Source file: tests/hadrons/Test_hadrons_quark.cc

 Copyright (C) 2017

 Author: Andrew Lawson <andrew.lawson1991@gmail.com>

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
 directory.
 *******************************************************************************/

#include "Test_hadrons.hpp"

using namespace Grid;
using namespace QCD;
using namespace Hadrons;

/*******************************************************************************
 * Unit test functions within Quark module.
 ******************************************************************************/

// Alternative 4D & 5D projections
template<class vobj>
inline void make_4D_with_gammas(Lattice<vobj> &in_5d, Lattice<vobj> &out_4d, int Ls)
{
    GridBase *_grid(out_4d._grid);
    Lattice<vobj> tmp(_grid);
    Gamma G5(Gamma::Algebra::Gamma5);

    ExtractSlice(tmp, in_5d, 0, 0);
    out_4d = 0.5 * (tmp - G5*tmp);
    ExtractSlice(tmp, in_5d, Ls - 1, 0);
    out_4d += 0.5 * (tmp + G5*tmp);
}

template<class vobj>
inline void make_5D_with_gammas(Lattice<vobj> &in_4d, Lattice<vobj> &out_5d, int Ls)
{
    out_5d = zero;
    Gamma G5(Gamma::Algebra::Gamma5);
    GridBase *_grid(in_4d._grid);
    Lattice<vobj> tmp(_grid);

    tmp = 0.5 * (in_4d + G5*in_4d);
    InsertSlice(tmp, out_5d, 0, 0);
    tmp = 0.5 * (in_4d - G5*in_4d);
    InsertSlice(tmp, out_5d, Ls - 1, 0);
}

int main(int argc, char **argv)
{
    /***************************************************************************
     * Initialisation.
     **************************************************************************/
    Grid_init(&argc, &argv);

    std::vector<int> latt_size   = GridDefaultLatt();
    std::vector<int> simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
    std::vector<int> mpi_layout  = GridDefaultMpi();

    const int Ls = 8;

    GridCartesian   UGrid(latt_size,simd_layout,mpi_layout);
    GridCartesian   *FGrid = SpaceTimeGrid::makeFiveDimGrid(Ls, &UGrid);
    GridSerialRNG   sRNG;
    GridParallelRNG pRNG(&UGrid);

    std::vector<int> seeds4({1,2,3,4});
    std::vector<int> seeds5({5,6,7,8});
    GridParallelRNG  rng4(&UGrid);
    GridParallelRNG  rng5(FGrid);
    rng4.SeedFixedIntegers(seeds4);
    rng5.SeedFixedIntegers(seeds5);

    /***************************************************************************
     * Build a 4D random source, and convert it to 5D.
     **************************************************************************/
    LatticeFermion test4(&UGrid);
    LatticeFermion test5(FGrid);
    LatticeFermion check5(FGrid);

    gaussian(rng4, test4);
    make_5D(test4, test5, Ls);
    make_5D_with_gammas(test4, check5, Ls);
    test5 -= check5;
    std::cout << "4D -> 5D comparison, diff = " << Grid::sqrt(norm2(test5)) << std::endl;

    /***************************************************************************
     * Build a 5D random source, and project down to 4D.
     **************************************************************************/
    LatticeFermion check4(&UGrid);
    gaussian(rng5, test5);
    check5 = test5;

    make_4D(test5, test4, Ls);
    make_4D_with_gammas(check5, check4, Ls);
    test4 -= check4;
    std::cout << "5D -> 4D comparison, diff = " << Grid::sqrt(norm2(test4)) << std::endl;

    /***************************************************************************
     * Convert a propagator to a fermion & back.
     **************************************************************************/
    LatticeFermion ferm(&UGrid);
    LatticePropagator prop(&UGrid), ref(&UGrid);
    gaussian(rng4, prop);

    // Define variables for sanity checking a single site.
    typename SpinColourVector::scalar_object fermSite;
    typename SpinColourMatrix::scalar_object propSite;
    std::vector<int> site(Nd, 0);

    for (int s = 0; s < Ns; ++s)
    for (int c = 0; c < Nc; ++c)
    {
        ref = prop;
        PropToFerm<WilsonImplR>(ferm, prop, s, c);
        FermToProp<WilsonImplR>(prop, ferm, s, c);

        std::cout << "Spin = " << s << ", Colour = " << c << std::endl;
        ref -= prop;
        std::cout << "Prop->Ferm->Prop test, diff = " << Grid::sqrt(norm2(ref)) << std::endl;

        peekSite(fermSite, ferm, site);
        peekSite(propSite, prop, site);
        for (int s2 = 0; s2 < Ns; ++s2)
        for (int c2 = 0; c2 < Nc; ++c2)
        {
            if (propSite()(s2, s)(c2, c) != fermSite()(s2)(c2))
            {
                std::cout << propSite()(s2, s)(c2, c) << " != "
                          << fermSite()(s2)(c2) << " for spin = " << s2
                          << ", col = " << c2 << std::endl;
            }
        }
    }

    Grid_finalize();
    return EXIT_SUCCESS;
}
