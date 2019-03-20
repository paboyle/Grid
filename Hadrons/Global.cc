/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Global.cc

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

#include <Hadrons/Global.hpp>

using namespace Grid;
using namespace QCD;
using namespace Hadrons;

HadronsLogger Hadrons::HadronsLogError(1,"Error");
HadronsLogger Hadrons::HadronsLogWarning(1,"Warning");
HadronsLogger Hadrons::HadronsLogMessage(1,"Message");
HadronsLogger Hadrons::HadronsLogIterative(1,"Iterative");
HadronsLogger Hadrons::HadronsLogDebug(1,"Debug");
HadronsLogger Hadrons::HadronsLogIRL(1,"IRL");

void Hadrons::initLogger(void)
{
    auto w  = std::string("Hadrons").length();
    int  cw = 8;


    GridLogError.setTopWidth(w);
    GridLogWarning.setTopWidth(w);
    GridLogMessage.setTopWidth(w);
    GridLogIterative.setTopWidth(w);
    GridLogDebug.setTopWidth(w);
    GridLogIRL.setTopWidth(w);
    GridLogError.setChanWidth(cw);
    GridLogWarning.setChanWidth(cw);
    GridLogMessage.setChanWidth(cw);
    GridLogIterative.setChanWidth(cw);
    GridLogDebug.setChanWidth(cw);
    GridLogIRL.setChanWidth(cw);
    HadronsLogError.Active(true);
    HadronsLogWarning.Active(true);
    HadronsLogMessage.Active(GridLogMessage.isActive());
    HadronsLogIterative.Active(GridLogIterative.isActive());
    HadronsLogDebug.Active(GridLogDebug.isActive());
    HadronsLogIRL.Active(GridLogIRL.isActive());
    HadronsLogError.setChanWidth(cw);
    HadronsLogWarning.setChanWidth(cw);
    HadronsLogMessage.setChanWidth(cw);
    HadronsLogIterative.setChanWidth(cw);
    HadronsLogDebug.setChanWidth(cw);
    HadronsLogIRL.setChanWidth(cw);
}

// type utilities //////////////////////////////////////////////////////////////
size_t Hadrons::typeHash(const std::type_info *info)
{
    return info->hash_code();
}

constexpr unsigned int maxNameSize = 1024u;

std::string Hadrons::typeName(const std::type_info *info)
{
    char        *buf;
    std::string name;
    
    buf  = abi::__cxa_demangle(info->name(), nullptr, nullptr, nullptr);
    name = buf;
    free(buf);
    
    return name;
}

// default writers/readers /////////////////////////////////////////////////////
#ifdef HAVE_HDF5
const std::string Hadrons::resultFileExt = "h5";
#else
const std::string Hadrons::resultFileExt = "xml";
#endif

// recursive mkdir /////////////////////////////////////////////////////////////
int Hadrons::mkdir(const std::string dirName)
{
    if (!dirName.empty() and access(dirName.c_str(), R_OK|W_OK|X_OK))
    {
        mode_t mode755;
        char   tmp[MAX_PATH_LENGTH];
        char   *p = NULL;
        size_t len;

        mode755 = S_IRWXU|S_IRGRP|S_IXGRP|S_IROTH|S_IXOTH;

        snprintf(tmp, sizeof(tmp), "%s", dirName.c_str());
        len = strlen(tmp);
        if(tmp[len - 1] == '/')
        {
            tmp[len - 1] = 0;
        }
        for(p = tmp + 1; *p; p++)
        {
            if(*p == '/')
            {
                *p = 0;
                ::mkdir(tmp, mode755);
                *p = '/';
            }
        }

        return ::mkdir(tmp, mode755);
    }
    else
    {
        return 0;
    }
}

std::string Hadrons::basename(const std::string &s)
{
    constexpr char sep = '/';
    size_t         i   = s.rfind(sep, s.length());
    
    if (i != std::string::npos)
    {
        return s.substr(i+1, s.length() - i);
    }
    else
    {
        return s;
    }
}

std::string Hadrons::dirname(const std::string &s)
{
    constexpr char sep = '/';
    size_t         i   = s.rfind(sep, s.length());
    
    if (i != std::string::npos)
    {
        return s.substr(0, i);
    }
    else
    {
        return "";
    }
}

void Hadrons::makeFileDir(const std::string filename, GridBase *g)
{
    bool doIt = true;

    if (g)
    {
        doIt = g->IsBoss();
    }
    if (doIt)
    {
        std::string dir    = dirname(filename);
        int         status = mkdir(dir);

        if (status)
        {
            HADRONS_ERROR(Io, "cannot create directory '" + dir
                          + "' ( " + std::strerror(errno) + ")");
        }
    }
}

void Hadrons::printTimeProfile(const std::map<std::string, GridTime> &timing, 
                               GridTime total)
{
    typedef decltype(total.count()) Count;

    std::map<Count, std::string, std::greater<Count>> rtiming;
    const double dtotal = static_cast<double>(total.count());
    auto cf = std::cout.flags();
    auto p  = std::cout.precision();
    unsigned int width = 0;

    for (auto &t: timing)
    {
        width = std::max(width, static_cast<unsigned int>(t.first.length()));
        rtiming[t.second.count()] = t.first;
    }
    for (auto &rt: rtiming)
    {
        LOG(Message) << std::setw(width) << rt.second << ": " 
                     << rt.first << " us (" << std::fixed 
                     << std::setprecision(1) 
                     << static_cast<double>(rt.first)/dtotal*100 << "%)"
                     << std::endl;
    }
    std::cout.flags(cf);
    std::cout.precision(p);
}
