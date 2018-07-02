/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: extras/Hadrons/Global.hpp

Copyright (C) 2015-2018

Author: Antonin Portelli <antonin.portelli@me.com>
Author: Lanny91 <andrew.lawson@gmail.com>

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

#ifndef Hadrons_Global_hpp_
#define Hadrons_Global_hpp_

#include <set>
#include <stack>
#include <Grid/Grid.h>
#include <cxxabi.h>

#ifndef SITE_SIZE_TYPE
#define SITE_SIZE_TYPE size_t
#endif

#ifndef DEFAULT_ASCII_PREC
#define DEFAULT_ASCII_PREC 16
#endif

/* the 'using Grid::operator<<;' statement prevents a very nasty compilation
 * error with GCC 5 (clang & GCC 6 compile fine without it).
 */

#define BEGIN_HADRONS_NAMESPACE \
namespace Grid {\
using namespace QCD;\
namespace Hadrons {\
using Grid::operator<<;\
using Grid::operator>>;
#define END_HADRONS_NAMESPACE }}

#define BEGIN_MODULE_NAMESPACE(name)\
namespace name {\
using Grid::operator<<;\
using Grid::operator>>;

#define END_MODULE_NAMESPACE }

#ifndef FIMPL
#define FIMPL WilsonImplR
#endif
#ifndef ZFIMPL
#define ZFIMPL ZWilsonImplR
#endif
#ifndef SIMPL
#define SIMPL ScalarImplCR
#endif
#ifndef GIMPL
#define GIMPL PeriodicGimplR
#endif

BEGIN_HADRONS_NAMESPACE

// type aliases
#define FERM_TYPE_ALIASES(FImpl, suffix)\
typedef FermionOperator<FImpl>                        FMat##suffix;            \
typedef typename FImpl::FermionField                  FermionField##suffix;    \
typedef typename FImpl::PropagatorField               PropagatorField##suffix; \
typedef typename FImpl::SitePropagator::scalar_object SitePropagator##suffix;  \
typedef std::vector<SitePropagator##suffix>           SlicedPropagator##suffix;

#define GAUGE_TYPE_ALIASES(FImpl, suffix)\
typedef typename FImpl::DoubledGaugeField DoubledGaugeField##suffix;

#define SCALAR_TYPE_ALIASES(SImpl, suffix)\
typedef typename SImpl::Field ScalarField##suffix;\
typedef typename SImpl::Field PropagatorField##suffix;

#define SOLVER_TYPE_ALIASES(FImpl, suffix)\
typedef Solver<FImpl> Solver##suffix;

#define SINK_TYPE_ALIASES(suffix)\
typedef std::function<SlicedPropagator##suffix\
                      (const PropagatorField##suffix &)> SinkFn##suffix;

#define FG_TYPE_ALIASES(FImpl, suffix)\
FERM_TYPE_ALIASES(FImpl, suffix)\
GAUGE_TYPE_ALIASES(FImpl, suffix)

// logger
class HadronsLogger: public Logger
{
public:
    HadronsLogger(int on, std::string nm): Logger("Hadrons", on, nm,
                                                  GridLogColours, "BLACK"){};
};

#define LOG(channel) std::cout << HadronsLog##channel
#define HADRONS_DEBUG_VAR(var) LOG(Debug) << #var << "= " << (var) << std::endl;

extern HadronsLogger HadronsLogError;
extern HadronsLogger HadronsLogWarning;
extern HadronsLogger HadronsLogMessage;
extern HadronsLogger HadronsLogIterative;
extern HadronsLogger HadronsLogDebug;
extern HadronsLogger HadronsLogIRL;

void initLogger(void);

// singleton pattern
#define SINGLETON(name)\
public:\
    name(const name &e) = delete;\
    void operator=(const name &e) = delete;\
    static name & getInstance(void)\
    {\
        static name e;\
        return e;\
    }\
private:\
    name(void);

#define SINGLETON_DEFCTOR(name)\
public:\
    name(const name &e) = delete;\
    void operator=(const name &e) = delete;\
    static name & getInstance(void)\
    {\
        static name e;\
        return e;\
    }\
private:\
    name(void) = default;

// type utilities
template <typename T>
const std::type_info * typeIdPt(const T &x)
{
    return &typeid(x);
}

std::string typeName(const std::type_info *info);

template <typename T>
const std::type_info * typeIdPt(void)
{
    return &typeid(T);
}

template <typename T>
std::string typeName(const T &x)
{
    return typeName(typeIdPt(x));
}

template <typename T>
std::string typeName(void)
{
    return typeName(typeIdPt<T>());
}

// default writers/readers
extern const std::string resultFileExt;

#ifdef HAVE_HDF5
typedef Hdf5Reader ResultReader;
typedef Hdf5Writer ResultWriter;
#else
typedef XmlReader ResultReader;
typedef XmlWriter ResultWriter;
#endif

#define RESULT_FILE_NAME(name) \
name + "." + std::to_string(vm().getTrajectory()) + "." + resultFileExt

// recursive mkdir
#define MAX_PATH_LENGTH 512u
int         mkdir(const std::string dirName);
std::string basename(const std::string &s);
std::string dirname(const std::string &s);
void        makeFileDir(const std::string filename, GridBase *g);

// default Schur convention
#ifndef HADRONS_DEFAULT_SCHUR 
#define HADRONS_DEFAULT_SCHUR DiagTwo
#endif
#define _HADRONS_SCHUR_OP_(conv) Schur##conv##Operator
#define HADRONS_SCHUR_OP(conv) _HADRONS_SCHUR_OP_(conv)
#define HADRONS_DEFAULT_SCHUR_OP HADRONS_SCHUR_OP(HADRONS_DEFAULT_SCHUR)
#define _HADRONS_SCHUR_SOLVE_(conv) SchurRedBlack##conv##Solve
#define HADRONS_SCHUR_SOLVE(conv) _HADRONS_SCHUR_SOLVE_(conv)
#define HADRONS_DEFAULT_SCHUR_SOLVE HADRONS_SCHUR_SOLVE(HADRONS_DEFAULT_SCHUR)

// stringify macro
#define _HADRONS_STR(x) #x
#define HADRONS_STR(x) _HADRONS_STR(x)

END_HADRONS_NAMESPACE

#include <Grid/Hadrons/Exceptions.hpp>

#endif // Hadrons_Global_hpp_
