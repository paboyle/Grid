// vim: set ts=2 sw=2 expandtab:

// Copyright (c) 2016 Luchang Jin
// All rights reserved.

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef INCLUDE_SHOW_H
#define INCLUDE_SHOW_H

#include <sstream>
#include <string>
#include <cstdarg>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <sstream>

#ifdef CURRENT_DEFAULT_NAMESPACE_NAME
namespace CURRENT_DEFAULT_NAMESPACE_NAME {
#endif

inline std::string vssprintf(const char* fmt, va_list args)
{
  std::string str;
  char* cstr;
  vasprintf(&cstr, fmt, args);
  str += std::string(cstr);
  std::free(cstr);
  return str;
}

inline std::string ssprintf(const char* fmt, ...)
{
  va_list args;
  va_start(args, fmt);
  return vssprintf(fmt, args);
}

inline std::string show()
{
  return "";
}

inline std::string show(const int& x)
{
  return ssprintf("%d", x);
}

inline std::string show(const unsigned int& x)
{
  return ssprintf("%u", x);
}

inline std::string show(const long& x)
{
  return ssprintf("%ld", x);
}

inline std::string show(const unsigned long& x)
{
  return ssprintf("%lu", x);
}

inline std::string show(const double& x)
{
  return ssprintf("%24.17E", x);
}

inline std::string show(const bool& x)
{
  return x ? "true" : "false";
}

inline std::string show(const std::string& x)
{
  std::ostringstream out;
  out << x;
  return out.str();
}

template <class T>
std::string shows(const T& x)
{
  std::ostringstream out;
  out << x;
  return out.str();
}

template <class T>
T& reads(T& x, const std::string& str)
{
  std::istringstream in(str);
  in >> x;
  return x;
}

inline FILE*& get_output_file()
{
  static FILE* out = stdout;
  return out;
}

inline void display(const std::string& str, FILE* fp = NULL)
{
  if (NULL == fp) {
    fp = get_output_file();
  }
  fprintf(fp, "%s", str.c_str());
}

inline void displayln(const std::string& str, FILE* fp = NULL)
{
  if (NULL == fp) {
    fp = get_output_file();
  }
  fprintf(fp, "%s\n", str.c_str());
}

//////////////////////////////////////////////////////////////////

inline void fdisplay(FILE* fp, const std::string& str)
{
  fprintf(fp, "%s", str.c_str());
}

inline void fdisplayln(FILE* fp, const std::string& str)
{
  fprintf(fp, "%s\n", str.c_str());
}

#ifdef CURRENT_DEFAULT_NAMESPACE_NAME
}
#endif

#endif
