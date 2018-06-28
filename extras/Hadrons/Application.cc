/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: extras/Hadrons/Application.cc

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

#include <Grid/Hadrons/Application.hpp>
#include <Grid/Hadrons/GeneticScheduler.hpp>
#include <Grid/Hadrons/Modules.hpp>

using namespace Grid;
using namespace QCD;
using namespace Hadrons;

#define BIG_SEP "==============="
#define SEP     "---------------"

/******************************************************************************
 *                       Application implementation                           *
 ******************************************************************************/
// constructors ////////////////////////////////////////////////////////////////
#define MACOUT(macro)    macro              << " (" << #macro << ")"
#define MACOUTS(macro) HADRONS_STR(macro) << " (" << #macro << ")"

Application::Application(void)
{
    initLogger();
    auto dim = GridDefaultLatt(), mpi = GridDefaultMpi(), loc(dim);
    locVol_ = 1;
    for (unsigned int d = 0; d < dim.size(); ++d)
    {
        loc[d]  /= mpi[d];
        locVol_ *= loc[d];
    }
    LOG(Message) << "====== HADRONS APPLICATION STARTING ======" << std::endl;
    LOG(Message) << "** Dimensions" << std::endl;
    LOG(Message) << "Global lattice       : " << dim << std::endl;
    LOG(Message) << "MPI partition        : " << mpi << std::endl;
    LOG(Message) << "Local lattice        : " << loc << std::endl;
    LOG(Message) << std::endl;
    LOG(Message) << "** Default parameters (and associated C macro)" << std::endl;
    LOG(Message) << "ASCII output precision  : " << MACOUT(DEFAULT_ASCII_PREC) << std::endl;
    LOG(Message) << "Fermion implementation  : " << MACOUTS(FIMPL) << std::endl;
    LOG(Message) << "z-Fermion implementation: " << MACOUTS(ZFIMPL) << std::endl;
    LOG(Message) << "Scalar implementation   : " << MACOUTS(SIMPL) << std::endl;
    LOG(Message) << "Gauge implementation    : " << MACOUTS(GIMPL) << std::endl;
    LOG(Message) << "Eigenvector base size   : " 
                 << MACOUT(HADRONS_DEFAULT_LANCZOS_NBASIS) << std::endl;
    LOG(Message) << "Schur decomposition     : " << MACOUTS(HADRONS_DEFAULT_SCHUR) << std::endl;
    LOG(Message) << std::endl;
}

Application::Application(const Application::GlobalPar &par)
: Application()
{
    setPar(par);
}

Application::Application(const std::string parameterFileName)
: Application()
{
    parameterFileName_ = parameterFileName;
}

// access //////////////////////////////////////////////////////////////////////
void Application::setPar(const Application::GlobalPar &par)
{
    par_ = par;
    env().setSeed(strToVec<int>(par_.seed));
}

const Application::GlobalPar & Application::getPar(void)
{
    return par_;
}

// execute /////////////////////////////////////////////////////////////////////
void Application::run(void)
{
    if (!parameterFileName_.empty() and (vm().getNModule() == 0))
    {
        parseParameterFile(parameterFileName_);
    }
    vm().printContent();
    env().printContent();
    schedule();
    printSchedule();
    configLoop();
}

// parse parameter file ////////////////////////////////////////////////////////
class ObjectId: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(ObjectId,
                                    std::string, name,
                                    std::string, type);
};

void Application::parseParameterFile(const std::string parameterFileName)
{
    XmlReader reader(parameterFileName);
    GlobalPar par;
    ObjectId  id;
    
    LOG(Message) << "Building application from '" << parameterFileName << "'..." << std::endl;
    read(reader, "parameters", par);
    setPar(par);
    if (!push(reader, "modules"))
    {
        HADRONS_ERROR(Parsing, "Cannot open node 'modules' in parameter file '" 
                              + parameterFileName + "'");
    }
    if (!push(reader, "module"))
    {
        HADRONS_ERROR(Parsing, "Cannot open node 'modules/module' in parameter file '" 
                              + parameterFileName + "'");
    }
    do
    {
        read(reader, "id", id);
        vm().createModule(id.name, id.type, reader);
    } while (reader.nextElement("module"));
    pop(reader);
    pop(reader);
}

void Application::saveParameterFile(const std::string parameterFileName)
{
    LOG(Message) << "Saving application to '" << parameterFileName << "'..." << std::endl;
    if (env().getGrid()->IsBoss())
    {
        XmlWriter          writer(parameterFileName);
        ObjectId           id;
        const unsigned int nMod = vm().getNModule();

        write(writer, "parameters", getPar());
        push(writer, "modules");
        for (unsigned int i = 0; i < nMod; ++i)
        {
            push(writer, "module");
            id.name = vm().getModuleName(i);
            id.type = vm().getModule(i)->getRegisteredName();
            write(writer, "id", id);
            vm().getModule(i)->saveParameters(writer, "options");
            pop(writer);
        }
        pop(writer);
        pop(writer);
    }
}

// schedule computation ////////////////////////////////////////////////////////
void Application::schedule(void)
{
    if (!scheduled_ and !loadedSchedule_)
    {
        program_   = vm().schedule(par_.genetic);
        scheduled_ = true;
    }
}

void Application::saveSchedule(const std::string filename)
{
    LOG(Message) << "Saving current schedule to '" << filename << "'..."
                 << std::endl;
    if (env().getGrid()->IsBoss())
    {
        TextWriter               writer(filename);
        std::vector<std::string> program;
        
        if (!scheduled_)
        {
            HADRONS_ERROR(Definition, "Computation not scheduled");
        }

        for (auto address: program_)
        {
            program.push_back(vm().getModuleName(address));
        }
        write(writer, "schedule", program);
    }
}

void Application::loadSchedule(const std::string filename)
{
    TextReader               reader(filename);
    std::vector<std::string> program;
    
    LOG(Message) << "Loading schedule from '" << filename << "'..."
                 << std::endl;
    read(reader, "schedule", program);
    program_.clear();
    for (auto &name: program)
    {
        program_.push_back(vm().getModuleAddress(name));
    }
    loadedSchedule_ = true;
}

void Application::printSchedule(void)
{
    if (!scheduled_)
    {
        HADRONS_ERROR(Definition, "Computation not scheduled");
    }
    auto peak = vm().memoryNeeded(program_);
    LOG(Message) << "Schedule (memory needed: " << sizeString(peak) << "):"
                 << std::endl;
    for (unsigned int i = 0; i < program_.size(); ++i)
    {
        LOG(Message) << std::setw(4) << i + 1 << ": "
                     << vm().getModuleName(program_[i]) << std::endl;
    }
}

// loop on configurations //////////////////////////////////////////////////////
void Application::configLoop(void)
{
    auto range = par_.trajCounter;
    
    for (unsigned int t = range.start; t < range.end; t += range.step)
    {
        LOG(Message) << BIG_SEP << " Starting measurement for trajectory " << t
                     << " " << BIG_SEP << std::endl;
        vm().setTrajectory(t);
        vm().executeProgram(program_);
    }
    LOG(Message) << BIG_SEP << " End of measurement " << BIG_SEP << std::endl;
    env().freeAll();
}
