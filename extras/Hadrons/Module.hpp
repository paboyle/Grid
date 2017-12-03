/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: extras/Hadrons/Module.hpp

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

#ifndef Hadrons_Module_hpp_
#define Hadrons_Module_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Environment.hpp>

BEGIN_HADRONS_NAMESPACE

// module registration macros
#define MODULE_REGISTER(mod, base)\
class mod: public base\
{\
public:\
    typedef base Base;\
    using Base::Base;\
    virtual std::string getRegisteredName(void)\
    {\
        return std::string(#mod);\
    }\
};\
class mod##ModuleRegistrar\
{\
public:\
    mod##ModuleRegistrar(void)\
    {\
        ModuleFactory &modFac = ModuleFactory::getInstance();\
        modFac.registerBuilder(#mod, [&](const std::string name)\
                              {\
                                  return std::unique_ptr<mod>(new mod(name));\
                              });\
    }\
};\
static mod##ModuleRegistrar mod##ModuleRegistrarInstance;

#define MODULE_REGISTER_NS(mod, base, ns)\
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

#define ARG(...) __VA_ARGS__
#define MACRO_REDIRECT(arg1, arg2, arg3, macro, ...) macro

#define envGet(type, name)\
*env().template getObject<type>(name)

#define envGetTmp(type, name)\
*env().template getObject<type>(getName() + "_tmp_" + name)

#define envHasType(type, name)\
env().template isObjectOfType<type>(name)

#define envCreate(type, name, Ls, pt)\
env().template createObject<type>(name, Environment::Storage::object, Ls, pt)

#define envCreateLat4(type, name)\
envCreate(type, name, 1, new type(env().getGrid()))

#define envCreateLat5(type, name, Ls)\
envCreate(type, name, Ls, new type(env().getGrid(Ls)))

#define envCreateLat(...)\
MACRO_REDIRECT(__VA_ARGS__, envCreateLat5, envCreateLat4)(__VA_ARGS__)

#define envCache(type, name, Ls, pt)\
env().template createObject<type>(name, Environment::Storage::cache, Ls, pt)

#define envCacheLat4(type, name)\
envCache(type, name, 1, new type(env().getGrid()))

#define envCacheLat5(type, name, Ls)\
envCache(type, name, Ls, new type(env().getGrid(Ls)))

#define envCacheLat(...)\
MACRO_REDIRECT(__VA_ARGS__, envCacheLat5, envCacheLat4)(__VA_ARGS__)

#define envTmp(type, name, Ls, pt)\
env().template createObject<type>(getName() + "_tmp_" + name,         \
                                  Environment::Storage::temporary, Ls, pt)

#define envTmpLat4(type, name)\
envTmp(type, name, 1, new type(env().getGrid()))

#define envTmpLat5(type, name, Ls)\
envTmp(type, name, Ls, new type(env().getGrid(Ls)))

#define envTmpLat(...)\
MACRO_REDIRECT(__VA_ARGS__, envTmpLat5, envTmpLat4)(__VA_ARGS__)

/******************************************************************************
 *                            Module class                                    *
 ******************************************************************************/
// base class
class ModuleBase
{
public:
    // constructor
    ModuleBase(const std::string name);
    // destructor
    virtual ~ModuleBase(void) = default;
    // access
    std::string getName(void) const;
    Environment &env(void) const;
    // get factory registration name if available
    virtual std::string getRegisteredName(void);
    // dependencies/products
    virtual std::vector<std::string> getInput(void) = 0;
    virtual std::vector<std::string> getOutput(void) = 0;
    // parse parameters
    virtual void parseParameters(XmlReader &reader, const std::string name) = 0;
    virtual void saveParameters(XmlWriter &writer, const std::string name) = 0;
    // setup
    virtual void setup(void) {};
    // execution
    void operator()(void);
    virtual void execute(void) = 0;
private:
    std::string name_;
    Environment &env_;
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
    // parameter access
    const P & par(void) const;
    void      setPar(const P &par);
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
