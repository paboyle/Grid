/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./Grid/util/Commandline.h

    Copyright (C) 2015 - 2020

    Author: Daniel Richtmann <daniel.richtmann@gmail.com>

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
#pragma once

NAMESPACE_BEGIN(Grid);

class Commandline {
public:
  static int readInt(int* argc, char*** argv, const std::string& option, int defaultValue) {
    std::string arg;
    int         ret = defaultValue;
    if(GridCmdOptionExists(*argv, *argv + *argc, option)) {
      arg = GridCmdOptionPayload(*argv, *argv + *argc, option);
      GridCmdOptionInt(arg, ret);
    }
    return ret;
  }

  static Coordinate readCoordinate(int*               argc,
                                   char***            argv,
                                   const std::string& option,
                                   const Coordinate&  defaultValues) {
    std::string      arg;
    std::vector<int> ret(defaultValues.toVector());
    if(GridCmdOptionExists(*argv, *argv + *argc, option)) {
      arg = GridCmdOptionPayload(*argv, *argv + *argc, option);
      GridCmdOptionIntVector(arg, ret);
    }
    return ret;
  }

  static bool readToggle(int* argc, char*** argv, const std::string& option) {
    return GridCmdOptionExists(*argv, *argv + *argc, option);
  }

  static std::vector<std::string> readCSL(int*                            argc,
                                          char***                         argv,
                                          const std::string&              option,
                                          const std::vector<std::string>& defaultValues) {
    std::string              arg;
    std::vector<std::string> ret(defaultValues);
    if(GridCmdOptionExists(*argv, *argv + *argc, option)) {
      arg = GridCmdOptionPayload(*argv, *argv + *argc, option);
      GridCmdOptionCSL(arg, ret);
    }
    return ret;
  }
};

NAMESPACE_END(Grid);
