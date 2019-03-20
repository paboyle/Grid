/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Module.hpp

Copyright (C) 2015-2019

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

#ifndef Hadrons_Module_hpp_
#define Hadrons_Module_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/TimerArray.hpp>
#include <Hadrons/VirtualMachine.hpp>

BEGIN_HADRONS_NAMESPACE

// module registration macros
#define MODULE_REGISTER(mod, base, ns)\
class mod: public base\
{\
public:\
    typedef base Base;\
    using Base::Base;\
    virtual std::string getRegisteredName(void)\
    {\
        return std::string(#ns "::" #mod);\
    }\
};\
class ns##mod##ModuleRegistrar\
{\
public:\
    ns##mod##ModuleRegistrar(void)\
    {\
        ModuleFactory &modFac = ModuleFactory::getInstance();\
        modFac.registerBuilder(#ns "::" #mod, [&](const std::string name)\
                              {\
                                  return std::unique_ptr<ns::mod>(new ns::mod(name));\
                              });\
    }\
};\
static ns##mod##ModuleRegistrar ns##mod##ModuleRegistrarInstance;

#define MODULE_REGISTER_TMP(mod, base, ns)\
extern template class base;\
MODULE_REGISTER(mod, ARG(base), ns);

#define HADRONS_MACRO_REDIRECT_12(arg1, arg2, macro, ...) macro
#define HADRONS_MACRO_REDIRECT_23(arg1, arg2, arg3, macro, ...) macro

#define envGetGrid4(latticeType)\
env().template getGrid<typename latticeType::vector_type>()

#define envGetGrid5(latticeType, Ls)\
env().template getGrid<typename latticeType::vector_type>(Ls)

#define envGetGrid(...)\
HADRONS_MACRO_REDIRECT_12(__VA_ARGS__, envGetGrid5, envGetGrid4)(__VA_ARGS__)

#define envGetCoarseGrid4(latticeType, blockSize)\
env().template getCoarseGrid<typename latticeType::vector_type>(blockSize)

#define envGetCoarseGrid5(latticeType, blockSize, Ls)\
env().template getCoarseGrid<typename latticeType::vector_type>(blockSize, Ls)

#define envGetCoarseGrid(...)\
HADRONS_MACRO_REDIRECT_23(__VA_ARGS__, envGetCoarseGrid5, envGetCoarseGrid4)(__VA_ARGS__)

#define envGetRbGrid4(latticeType)\
env().template getRbGrid<typename latticeType::vector_type>()

#define envGetRbGrid5(latticeType, Ls)\
env().template getRbGrid<typename latticeType::vector_type>(Ls)

#define envGetRbGrid(...)\
HADRONS_MACRO_REDIRECT_12(__VA_ARGS__, envGetRbGrid5, envGetRbGrid4)(__VA_ARGS__)

#define envGet(type, name)\
*env().template getObject<type>(name)

#define envGetDerived(base, type, name)\
*env().template getDerivedObject<base, type>(name)

#define envGetTmp(type, var)\
type &var = *env().template getObject<type>(getName() + "_tmp_" + #var)

#define envHasType(type, name)\
env().template isObjectOfType<type>(name)

#define envCreate(type, name, Ls, ...)\
env().template createObject<type>(name, Environment::Storage::object, Ls, __VA_ARGS__)

#define envCreateDerived(base, type, name, Ls, ...)\
env().template createDerivedObject<base, type>(name, Environment::Storage::object, Ls, __VA_ARGS__)

#define envCreateLat4(type, name)\
envCreate(type, name, 1, envGetGrid(type))

#define envCreateLat5(type, name, Ls)\
envCreate(type, name, Ls, envGetGrid(type, Ls))

#define envCreateLat(...)\
HADRONS_MACRO_REDIRECT_23(__VA_ARGS__, envCreateLat5, envCreateLat4)(__VA_ARGS__)

#define envCache(type, name, Ls, ...)\
env().template createObject<type>(name, Environment::Storage::cache, Ls, __VA_ARGS__)

#define envCacheLat4(type, name)\
envCache(type, name, 1, envGetGrid(type))

#define envCacheLat5(type, name, Ls)\
envCache(type, name, Ls, envGetGrid(type, Ls))

#define envCacheLat(...)\
HADRONS_MACRO_REDIRECT_23(__VA_ARGS__, envCacheLat5, envCacheLat4)(__VA_ARGS__)

#define envTmp(type, name, Ls, ...)\
env().template createObject<type>(getName() + "_tmp_" + name,         \
                                  Environment::Storage::temporary, Ls, __VA_ARGS__)

#define envTmpLat4(type, name)\
envTmp(type, name, 1, envGetGrid(type))

#define envTmpLat5(type, name, Ls)\
envTmp(type, name, Ls, envGetGrid(type, Ls))

#define envTmpLat(...)\
HADRONS_MACRO_REDIRECT_23(__VA_ARGS__, envTmpLat5, envTmpLat4)(__VA_ARGS__)

#define saveResult(ioStem, name, result)\
if (env().getGrid()->IsBoss() and !ioStem.empty())\
{\
    makeFileDir(ioStem, env().getGrid());\
    {\
        ResultWriter _writer(RESULT_FILE_NAME(ioStem, vm().getTrajectory()));\
        write(_writer, name, result);\
    }\
}

/******************************************************************************
 *                            Module class                                    *
 ******************************************************************************/
// base class
class ModuleBase: public TimerArray
{
public:
    // constructor
    ModuleBase(const std::string name);
    // destructor
    virtual ~ModuleBase(void) = default;
    // access
    std::string getName(void) const;
    // get factory registration name if available
    virtual std::string getRegisteredName(void);
    // dependencies/products
    virtual std::vector<std::string> getInput(void) = 0;
    virtual std::vector<std::string> getReference(void)
    {
        return std::vector<std::string>(0);
    };
    virtual std::vector<std::string> getOutput(void) = 0;
    // parse parameters
    virtual void parseParameters(XmlReader &reader, const std::string name) = 0;
    virtual void saveParameters(XmlWriter &writer, const std::string name) = 0;
    // parameter string
    virtual std::string parString(void) const = 0;
    // setup
    virtual void setup(void) {};
    virtual void execute(void) = 0;
    // execution
    void operator()(void);
protected:
    // environment shortcut
    DEFINE_ENV_ALIAS;
    // virtual machine shortcut
    DEFINE_VM_ALIAS;
    // RNG seeded from module string
    GridParallelRNG &rng4d(void);
private:
    std::string makeSeedString(void);
private:
    std::string                          name_, currentTimer_, seed_;
    std::map<std::string, GridStopWatch> timer_; 
};

// derived class, templating the parameter class
template <typename P>
class Module: public ModuleBase
{
public:
    typedef P Par;
public:
    // constructor
    Module(const std::string name);
    // destructor
    virtual ~Module(void) = default;
    // parse parameters
    virtual void parseParameters(XmlReader &reader, const std::string name);
    virtual void saveParameters(XmlWriter &writer, const std::string name);
    // parameter string
    virtual std::string parString(void) const;
    // parameter access
    const P &   par(void) const;
    void        setPar(const P &par);
private:
    P par_;
};

// no parameter type
class NoPar {};

template <>
class Module<NoPar>: public ModuleBase
{
public:
    // constructor
    Module(const std::string name): ModuleBase(name) {};
    // destructor
    virtual ~Module(void) = default;
    // parse parameters (do nothing)
    virtual void parseParameters(XmlReader &reader, const std::string name) {};
    virtual void saveParameters(XmlWriter &writer, const std::string name)
    {
        push(writer, "options");
        pop(writer);
    };
    // parameter string (empty)
    virtual std::string parString(void) const {return "";};
};

/******************************************************************************
 *                           Template implementation                          *
 ******************************************************************************/
template <typename P>
Module<P>::Module(const std::string name)
: ModuleBase(name)
{}

template <typename P>
void Module<P>::parseParameters(XmlReader &reader, const std::string name)
{
    read(reader, name, par_);
}

template <typename P>
void Module<P>::saveParameters(XmlWriter &writer, const std::string name)
{
    write(writer, name, par_);
}

template <typename P>
std::string Module<P>::parString(void) const
{
    XmlWriter writer("", "");

    write(writer, par_.SerialisableClassName(), par_);

    return writer.string();
}

template <typename P>
const P & Module<P>::par(void) const
{
    return par_;
}

template <typename P>
void Module<P>::setPar(const P &par)
{
    par_ = par;
}

END_HADRONS_NAMESPACE

#endif // Hadrons_Module_hpp_
