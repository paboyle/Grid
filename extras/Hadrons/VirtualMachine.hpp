/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: extras/Hadrons/VirtualMachine.hpp

Copyright (C) 2017

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

#ifndef Hadrons_VirtualMachine_hpp_
#define Hadrons_VirtualMachine_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Graph.hpp>
#include <Grid/Hadrons/Environment.hpp>

BEGIN_HADRONS_NAMESPACE

#define DEFINE_VM_ALIAS \
inline VirtualMachine & vm(void) const\
{\
    return VirtualMachine::getInstance();\
}

/******************************************************************************
 *                   Virtual machine for module execution                     *
 ******************************************************************************/
// forward declaration of Module
class ModuleBase;

class VirtualMachine
{
    SINGLETON_DEFCTOR(VirtualMachine);
public:
    typedef SITE_SIZE_TYPE                             Size;
    typedef std::unique_ptr<ModuleBase>                ModPt;
    struct MemoryPrint
    {
        Size         size;
        unsigned int module;
    };
    struct MemoryProfile
    {
        std::vector<std::map<unsigned int, Size>> module;
        std::vector<MemoryPrint>                  object;
    };
private:
    struct ModuleInfo
    {
        const std::type_info      *type{nullptr};
        std::string               name;
        ModPt                     data{nullptr};
        std::vector<unsigned int> input;
        size_t                    maxAllocated;
    };
public:
    // dry run
    void                dryRun(const bool isDry);
    bool                isDryRun(void) const;
    void                memoryProfile(const bool doMemoryProfile);
    bool                doMemoryProfile(void) const;
    // trajectory counter
    void                setTrajectory(const unsigned int traj);
    unsigned int        getTrajectory(void) const;
    // module management
    void                pushModule(ModPt &pt);
    template <typename M>
    void                createModule(const std::string name);
    template <typename M>
    void                createModule(const std::string name,
                                         const typename M::Par &par);
    void                createModule(const std::string name,
                                         const std::string type,
                                         XmlReader &reader);
    unsigned int        getNModule(void) const;
    ModuleBase *        getModule(const unsigned int address) const;
    ModuleBase *        getModule(const std::string name) const;
    template <typename M>
    M *                 getModule(const unsigned int address) const;
    template <typename M>
    M *                 getModule(const std::string name) const;
    unsigned int        getModuleAddress(const std::string name) const;
    std::string         getModuleName(const unsigned int address) const;
    std::string         getModuleType(const unsigned int address) const;
    std::string         getModuleType(const std::string name) const;
    std::string         getModuleNamespace(const unsigned int address) const;
    std::string         getModuleNamespace(const std::string name) const;
    bool                hasModule(const unsigned int address) const;
    bool                hasModule(const std::string name) const;
    Graph<unsigned int> makeModuleGraph(void) const;
    void                checkGraph(void) const;
    // print VM content
    void                printContent(void) const;
    // memory profile
    MemoryProfile       memoryProfile(void) const;
    // general execution
    Size                executeProgram(const std::vector<unsigned int> &p);
    Size                executeProgram(const std::vector<std::string> &p);
private:
    // environment shortcut
    DEFINE_ENV_ALIAS;
    // memory profile
    void resizeProfile(MemoryProfile &profile) const;
    void updateProfile(MemoryProfile &profile, const unsigned int address) const;
    void cleanEnvironment(MemoryProfile &profile) const;
    void memoryProfile(MemoryProfile &profile, const std::string name) const;
    void memoryProfile(MemoryProfile &profile, const unsigned int address) const;
private:
    // general
    bool                                dryRun_{false}, memoryProfile_{false};
    unsigned int                        traj_;
    // module and related maps
    std::vector<ModuleInfo>             module_;
    std::map<std::string, unsigned int> moduleAddress_;
    std::string                         currentModule_{""};
};

/******************************************************************************
 *                   VirtualMachine template implementation                   *
 ******************************************************************************/
// module management ///////////////////////////////////////////////////////////
template <typename M>
void VirtualMachine::createModule(const std::string name)
{
    ModPt pt(new M(name));
    
    pushModule(pt);
}

template <typename M>
void VirtualMachine::createModule(const std::string name,
                               const typename M::Par &par)
{
    ModPt pt(new M(name));
    
    static_cast<M *>(pt.get())->setPar(par);
    pushModule(pt);
}

template <typename M>
M * VirtualMachine::getModule(const unsigned int address) const
{
    if (auto *pt = dynamic_cast<M *>(getModule(address)))
    {
        return pt;
    }
    else
    {
        HADRON_ERROR(Definition, "module '" + module_[address].name
                     + "' does not have type " + typeid(M).name()
                     + "(has type: " + getModuleType(address) + ")");
    }
}

template <typename M>
M * VirtualMachine::getModule(const std::string name) const
{
    return getModule<M>(getModuleAddress(name));
}

END_HADRONS_NAMESPACE

#endif // Hadrons_VirtualMachine_hpp_
