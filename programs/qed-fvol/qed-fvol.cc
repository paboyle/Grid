#include <Global.hpp>

using namespace Grid;
using namespace QCD;
using namespace QedFVol;

int main(int argc, char *argv[])
{
    // parse command line
    std::string parameterFileName;
    
    if (argc < 2)
    {
        std::cerr << "usage: " << argv[0] << " <parameter file> [Grid options]";
        std::cerr << std::endl;
        std::exit(EXIT_FAILURE);
    }
    parameterFileName = argv[1];
    
    // initialization
    Grid_init(&argc, &argv);
    QedFVolLogError.Active(GridLogError.isActive());
    QedFVolLogWarning.Active(GridLogWarning.isActive());
    QedFVolLogMessage.Active(GridLogMessage.isActive());
    QedFVolLogIterative.Active(GridLogIterative.isActive());
    QedFVolLogDebug.Active(GridLogDebug.isActive());
    LOG(Message) << "Grid initialized" << std::endl;
    

    
    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();
    
    return EXIT_SUCCESS;
}
