#include <Hadrons/EigenPack.hpp>
#include <Hadrons/Environment.hpp>

using namespace Grid;
using namespace QCD;
using namespace Hadrons;

template <typename FOut, typename FIn>
void convert(const std::string outFilename, const std::string inFilename, 
             const unsigned int size, const bool multiFile)
{
    assert(outFilename != inFilename);
    
    typedef EigenPack<FOut>            EPOut;
    typedef EigenPack<FIn>             EPIn;
    typedef typename FOut::vector_type VTypeOut;
    typedef typename FIn::vector_type  VTypeIn;

    std::shared_ptr<GridCartesian>         fgIn, fgOut;
    std::shared_ptr<GridRedBlackCartesian> gIn, gOut;

    auto                      dim     = GridDefaultLatt();
    unsigned int              nd      = dim.size();
    auto                      simdOut = GridDefaultSimd(nd, VTypeOut::Nsimd());
    auto                      simdIn  = GridDefaultSimd(nd, VTypeIn::Nsimd());

    fgOut.reset(SpaceTimeGrid::makeFourDimGrid(dim, simdOut, GridDefaultMpi()));
    gOut.reset(SpaceTimeGrid::makeFourDimRedBlackGrid(fgOut.get()));
    fgIn.reset(SpaceTimeGrid::makeFourDimGrid(dim, simdIn, GridDefaultMpi()));
    gIn.reset(SpaceTimeGrid::makeFourDimRedBlackGrid(fgIn.get()));

    FOut bufOut(gOut.get());
    FIn  bufIn(gIn.get()), testIn(gIn.get());

    LOG(Message) << "==== EIGENPACK CONVERSION" << std::endl;
    LOG(Message) << "In path  : " << inFilename  << std::endl;
    LOG(Message) << "In type  : " << typeName<FIn>() << std::endl;
    LOG(Message) << "Out path : " << outFilename << std::endl;
    LOG(Message) << "Out type : " << typeName<FOut>() << std::endl;
    LOG(Message) << "Multifile: " << (multiFile ? "yes" : "no") << std::endl;
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
            makeFileDir(outV, gOut.get());
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

        makeFileDir(outFilename, gOut.get());
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
    unsigned int size;
    bool         multiFile;
    
    if (argc < 5)
    {
        std::cerr << "usage: " << argv[0] << " <out eigenpack> <in eigenpack> <size> <multifile (0|1)> [Grid options]";
        std::cerr << std::endl;
        std::exit(EXIT_FAILURE);
    }
    outFilename = argv[1];
    inFilename  = argv[2];
    size        = std::stoi(std::string(argv[3]));
    multiFile   = (std::string(argv[4]) != "0");
    
    // initialization
    Grid_init(&argc, &argv);

    // execution
    try
    {
        convert<FOUT, FIN>(outFilename, inFilename, size, multiFile);
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
