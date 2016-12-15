#ifndef QedFVol_Global_hpp_
#define QedFVol_Global_hpp_

#include <Grid/Grid.h>

#define BEGIN_QEDFVOL_NAMESPACE \
namespace Grid {\
using namespace QCD;\
namespace QedFVol {\
using Grid::operator<<;
#define END_QEDFVOL_NAMESPACE }}

/* the 'using Grid::operator<<;' statement prevents a very nasty compilation
 * error with GCC (clang compiles fine without it).
 */

BEGIN_QEDFVOL_NAMESPACE

class QedFVolLogger: public Logger
{
public:
    QedFVolLogger(int on, std::string nm): Logger("QedFVol", on, nm,
                                                  GridLogColours, "BLACK"){};
};

#define LOG(channel) std::cout << QedFVolLog##channel
#define QEDFVOL_ERROR(msg)\
LOG(Error) << msg << " (" << __FUNCTION__ << " at " << __FILE__ << ":"\
           << __LINE__ << ")" << std::endl;\
abort();

#define DEBUG_VAR(var) LOG(Debug) << #var << "= " << (var) << std::endl;

extern QedFVolLogger QedFVolLogError;
extern QedFVolLogger QedFVolLogWarning;
extern QedFVolLogger QedFVolLogMessage;
extern QedFVolLogger QedFVolLogIterative;
extern QedFVolLogger QedFVolLogDebug;

END_QEDFVOL_NAMESPACE

#endif // QedFVol_Global_hpp_
