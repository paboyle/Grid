#include <Hadrons/EigenPack.hpp>
#include <Hadrons/Environment.hpp>

using namespace Grid;
using namespace QCD;
using namespace Hadrons;

template <typename FOut, typename FIn>
void convert(const std::string outFilename, const std::string inFilename, 
             const unsigned int Ls, const bool rb, const unsigned int size, 
             const bool multiFile)
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

    FOut bufOut(gOut);
    FIn  bufIn(gIn), testIn(gIn);

    LOG(Message) << "==== EIGENPACK CONVERSION" << std::endl;
    LOG(Message) << "Lattice       : " << gIn->GlobalDimensions() << std::endl;
    LOG(Message) << "Checkerboarded: " << (rb ? "yes" : "no") << std::endl;
    LOG(Message) << "In path       : " << inFilename  << std::endl;
    LOG(Message) << "In type       : " << typeName<FIn>() << std::endl;
    LOG(Message) << "Out path      : " << outFilename << std::endl;
    LOG(Message) << "Out type      : " << typeName<FOut>() << std::endl;
    LOG(Message) << "#vectors      : " << size << std::endl;
    LOG(Message) << "Multifile     : " << (multiFile ? "yes" : "no") << std::endl;
    if (multiFile)
    {
        for(unsigned int k = 0; k < size; ++k)
        {
            ScidacWriter binWriter(gOut->IsBoss());
            ScidacReader binReader;
            PackRecord   record;
            VecRecord    vecRecord;
            std::string  outV = outFilename + "/v" + std::to_string(k) + ".bin";
            std::string  inV  = inFilename + "/v" + std::to_string(k) + ".bin";

            LOG(Message) << "==== Converting vector " << k << std::endl;
            LOG(Message) << "In : " << inV  << std::endl;
            LOG(Message) << "Out: " << outV << std::endl;
            makeFileDir(outV, gOut);
            binWriter.open(outV);
            binReader.open(inV);
            EPIn::readHeader(record, binReader);
            EPOut::writeHeader(binWriter, record);
            EPIn::readElement(bufIn, vecRecord, binReader);
            precisionChange(bufOut, bufIn);
            precisionChange(testIn, bufOut);
            testIn -= bufIn;
            LOG(Message) << "Diff norm^2: " << norm2(testIn) << std::endl;
            EPOut::writeElement(binWriter, bufOut, vecRecord);
            binWriter.close();
            binReader.close();
        }
    }
    else
    {
        ScidacWriter binWriter(gOut->IsBoss());
        ScidacReader binReader;
        PackRecord   record;

        makeFileDir(outFilename, gOut);
        binWriter.open(outFilename);
        binReader.open(inFilename);
        EPIn::readHeader(record, binReader);
        EPOut::writeHeader(binWriter, record);
        for(unsigned int k = 0; k < size; ++k)
        {
            VecRecord vecRecord;

            LOG(Message) << "==== Converting vector " << k << std::endl;
            EPIn::readElement(bufIn, vecRecord, binReader);
            precisionChange(bufOut, bufIn);
            precisionChange(testIn, bufOut);
            testIn -= bufIn;
            LOG(Message) << "Diff norm^2: " << norm2(testIn) << std::endl;
            EPOut::writeElement(binWriter, bufOut, vecRecord);
        }
        binWriter.close();
        binReader.close();
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
    bool         rb, multiFile;
    
    if (argc < 7)
    {
        std::cerr << "usage: " << argv[0] << " <out eigenpack> <in eigenpack> <Ls> <red-black (0|1)> <#vector> <multifile (0|1)> [Grid options]";
        std::cerr << std::endl;
        std::exit(EXIT_FAILURE);
    }
    outFilename = argv[1];
    inFilename  = argv[2];
    Ls          = std::stoi(std::string(argv[3]));
    rb          = (std::string(argv[4]) != "0");
    size        = std::stoi(std::string(argv[5]));
    multiFile   = (std::string(argv[6]) != "0");
    
    // initialization
    Grid_init(&argc, &argv);
    initLogger();

    // execution
    try
    {
        convert<FOUT, FIN>(outFilename, inFilename, Ls, rb, size, multiFile);
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
