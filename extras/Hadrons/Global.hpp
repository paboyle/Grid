/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: extras/Hadrons/Global.hpp

Copyright (C) 2015
Copyright (C) 2016

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

#ifndef Hadrons_Global_hpp_
#define Hadrons_Global_hpp_

#include <set>
#include <stack>
#include <Grid/Grid.h>
#include <cxxabi.h>

#ifndef SITE_SIZE_TYPE
#define SITE_SIZE_TYPE size_t
#endif

#define BEGIN_HADRONS_NAMESPACE \
namespace Grid {\
 \
namespace Hadrons {\
using Grid::operator<<;
#define END_HADRONS_NAMESPACE }}

#define BEGIN_MODULE_NAMESPACE(name)\
namespace name {\
using Grid::operator<<;
#define END_MODULE_NAMESPACE }

/* the 'using Grid::operator<<;' statement prevents a very nasty compilation
 * error with GCC 5 (clang & GCC 6 compile fine without it).
 */

#ifndef FIMPL
#define FIMPL WilsonImplR
#endif
#ifndef SIMPL
#define SIMPL ScalarImplCR
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
typedef std::function<void(FermionField##suffix &,\
                      const FermionField##suffix &)> SolverFn##suffix;

#define SINK_TYPE_ALIASES(suffix)\
typedef std::function<SlicedPropagator##suffix(const PropagatorField##suffix &)> SinkFn##suffix;

#define FGS_TYPE_ALIASES(FImpl, suffix)\
FERM_TYPE_ALIASES(FImpl, suffix)\
GAUGE_TYPE_ALIASES(FImpl, suffix)\
SOLVER_TYPE_ALIASES(FImpl, suffix)

// logger
class HadronsLogger: public Logger
{
public:
    HadronsLogger(int on, std::string nm): Logger("Hadrons", on, nm,
                                                  GridLogColours, "BLACK"){};
};

#define LOG(channel) std::cout << HadronsLog##channel
#define DEBUG_VAR(var) LOG(Debug) << #var << "= " << (var) << std::endl;

extern HadronsLogger HadronsLogError;
extern HadronsLogger HadronsLogWarning;
extern HadronsLogger HadronsLogMessage;
extern HadronsLogger HadronsLogIterative;
extern HadronsLogger HadronsLogDebug;

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
#ifdef HAVE_HDF5
typedef Hdf5Reader CorrReader;
typedef Hdf5Writer CorrWriter;
#else
typedef XmlReader CorrReader;
typedef XmlWriter CorrWriter;
#endif

END_HADRONS_NAMESPACE

#include <Grid/Hadrons/Exceptions.hpp>

#endif // Hadrons_Global_hpp_
