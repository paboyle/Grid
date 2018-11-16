/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Utilities/EigenPackCast.cc

Copyright (C) 2015-2018

Author: Antonin Portelli <antonin.portelli@me.com>

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
#include <Hadrons/EigenPack.hpp>
#include <Hadrons/Environment.hpp>

using namespace Grid;
using namespace QCD;
using namespace Hadrons;

template <typename FOut, typename FIn>
void convert(const std::string outFilename, const std::string inFilename, 
             const unsigned int Ls, const bool rb, const unsigned int size, 
             const bool multiFile, const bool testRead)
{
    assert(outFilename != inFilename);
    
    typedef EigenPack<FOut>            EPOut;
    typedef EigenPack<FIn>             EPIn;
    typedef typename FOut::vector_type VTypeOut;
    typedef typename FIn::vector_type  VTypeIn;

    std::shared_ptr<GridCartesian>         gInBase, gOutBase, gIn5, gOut5;
    std::shared_ptr<GridRedBlackCartesian> rbgIn, rbgOut;
    GridBase                               *gIn, *gOut;

    auto         dim     = GridDefaultLatt();
    unsigned int nd      = dim.size();
    auto         simdOut = GridDefaultSimd(nd, VTypeOut::Nsimd());
    auto         simdIn  = GridDefaultSimd(nd, VTypeIn::Nsimd());

    gOutBase.reset(SpaceTimeGrid::makeFourDimGrid(dim, simdOut, GridDefaultMpi()));
    gInBase.reset(SpaceTimeGrid::makeFourDimGrid(dim, simdIn, GridDefaultMpi()));
    if (rb)
    {
        if (Ls > 1)
        {
            rbgOut.reset(SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls, gOutBase.get()));
            rbgIn.reset(SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls, gInBase.get()));
        }
        else
        {
            rbgOut.reset(SpaceTimeGrid::makeFourDimRedBlackGrid(gOutBase.get()));
            rbgIn.reset(SpaceTimeGrid::makeFourDimRedBlackGrid(gInBase.get()));
        }
        gOut = rbgOut.get();
        gIn  = rbgIn.get();
    }
    else
    {
        if (Ls > 1)
        {
            gOut5.reset(SpaceTimeGrid::makeFiveDimGrid(Ls, gOutBase.get()));
            gIn5.reset(SpaceTimeGrid::makeFiveDimGrid(Ls, gInBase.get()));
            gOut = gOut5.get();
            gIn  = gIn5.get();
        }
        else
        {
            gOut = gOutBase.get();
            gIn  = gInBase.get();
        }
    }

    FOut         bufOut(gOut);
    FIn          bufIn(gIn), testIn(gIn);
    ScidacWriter binWriter(gOut->IsBoss());
    ScidacReader binReader;
    PackRecord   record;
    RealD        eval;

    LOG(Message) << "==== EIGENPACK CONVERSION" << std::endl;
    LOG(Message) << "Lattice       : " << gIn->GlobalDimensions() << std::endl;
    LOG(Message) << "Checkerboarded: " << (rb ? "yes" : "no") << std::endl;
    LOG(Message) << "In path       : " << inFilename  << std::endl;
    LOG(Message) << "In type       : " << typeName<FIn>() << std::endl;
    LOG(Message) << "Out path      : " << outFilename << std::endl;
    LOG(Message) << "Out type      : " << typeName<FOut>() << std::endl;
    LOG(Message) << "#vectors      : " << size << std::endl;
    LOG(Message) << "Multifile     : " << (multiFile ? "yes" : "no") << std::endl;
    LOG(Message) << "Test read     : " << (testRead ? "yes" : "no") << std::endl;
    if (multiFile)
    {
        for(unsigned int k = 0; k < size; ++k)
        {
            std::string  outV = outFilename + "/v" + std::to_string(k) + ".bin";
            std::string  inV  = inFilename + "/v" + std::to_string(k) + ".bin";

            LOG(Message) << "==== Converting vector " << k << std::endl;
            LOG(Message) << "In : " << inV  << std::endl;
            LOG(Message) << "Out: " << outV << std::endl;
            // conversion
            LOG(Message) << "-- Doing conversion" << std::endl;
            makeFileDir(outV, gOut);
            binWriter.open(outV);
            binReader.open(inV);
            EigenPackIo::readHeader(record, binReader);
            EigenPackIo::writeHeader(binWriter, record);
            EigenPackIo::readElement<FIn>(bufIn, eval, k, binReader);
            EigenPackIo::writeElement<FIn, FOut>(binWriter, bufIn, eval, k, &bufOut, &testIn);
            binWriter.close();
            binReader.close();
            // read test
            if (testRead)
            {
                LOG(Message) << "-- Test read" << std::endl;
                binReader.open(outV);
                EigenPackIo::readElement<FOut>(bufOut, eval, k, binReader);
                binReader.close();
            }
        }
    }
    else
    {
        // conversion
        LOG(Message) << "-- Doing conversion" << std::endl;
        makeFileDir(outFilename, gOut);
        binWriter.open(outFilename);
        binReader.open(inFilename);
        EigenPackIo::readHeader(record, binReader);
        EigenPackIo::writeHeader(binWriter, record);
        for(unsigned int k = 0; k < size; ++k)
        {
            EigenPackIo::readElement<FIn>(bufIn, eval, k, binReader);
            EigenPackIo::writeElement<FIn, FOut>(binWriter, bufIn, eval, k, &bufOut, &testIn);
        }
        binWriter.close();
        binReader.close();
        // read test
        if (testRead)
        {
            LOG(Message) << "-- Test read" << std::endl;
            binReader.open(outFilename);
            EigenPackIo::readHeader(record, binReader);
            for(unsigned int k = 0; k < size; ++k)
            {
                EigenPackIo::readElement<FOut>(bufOut, eval, k, binReader);
            }
            binReader.close();
        }
    }
}

#ifndef FOUT
#warning "FOUT undefined (set to WilsonImplF::FermionField by default)"
#define FOUT WilsonImplF::FermionField
#endif
#ifndef FIN
#warning "FIN undefined (set to WilsonImplD::FermionField by default)"
#define FIN WilsonImplD::FermionField
#endif

int main(int argc, char *argv[])
{
    // parse command line
    std::string  outFilename, inFilename;
    unsigned int size, Ls;
    bool         rb, multiFile, testRead;
    
    if (argc < 8)
    {
        std::cerr << "usage: " << argv[0] << " <out eigenpack> <in eigenpack> <Ls> <red-black {0|1}> <#vector> <multifile {0|1}> <test read {0|1}> [Grid options]";
        std::cerr << std::endl;
        std::exit(EXIT_FAILURE);
    }
    outFilename = argv[1];
    inFilename  = argv[2];
    Ls          = std::stoi(std::string(argv[3]));
    rb          = (std::string(argv[4]) != "0");
    size        = std::stoi(std::string(argv[5]));
    multiFile   = (std::string(argv[6]) != "0");
    testRead    = (std::string(argv[7]) != "0");
    
    // initialization
    Grid_init(&argc, &argv);
    initLogger();

    // execution
    try
    {
        convert<FOUT, FIN>(outFilename, inFilename, Ls, rb, size, multiFile, testRead);
    }
    catch (const std::exception& e)
    {
        Exceptions::abort(e);
    }

    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();
    
    return EXIT_SUCCESS;
}
