   /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/util/Sha.h

    Copyright (C) 2018

    Author: Peter Boyle

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
extern "C" {
#include <openssl/sha.h>
}

#pragma once

class GridChecksum
{
public:
  static inline uint32_t crc32(void *data,size_t bytes)
  {
    return ::crc32(0L,(unsigned char *)data,bytes);
  }
  template <typename T>
  static inline std::string sha256_string(const std::vector<T> &hash)
  {
    std::stringstream sha;
    std::string       s;

    for(unsigned int i = 0; i < hash.size(); i++) 
    { 
        sha << std::hex << static_cast<unsigned int>(hash[i]);
    }
    s = sha.str();

    return s;
  }
  static inline std::vector<unsigned char> sha256(const void *data,size_t bytes)
  {
    std::vector<unsigned char> hash(SHA256_DIGEST_LENGTH);
    SHA256_CTX sha256;
    SHA256_Init  (&sha256);
    SHA256_Update(&sha256, data,bytes);
    SHA256_Final (&hash[0], &sha256);
    return hash;
  }
  static inline std::vector<int> sha256_seeds(const std::string &s)
  {
    std::vector<int> seeds;
    std::vector<unsigned char> uchars = sha256((void *)s.c_str(),s.size());
    for(int i=0;i<uchars.size();i++) seeds.push_back(uchars[i]);
    return seeds;
  }
};

/*
int main(int argc,char **argv)
{
  std::string s("The quick brown fox jumps over the lazy dog");
  auto csum = GridChecksum::sha256_seeds(s);
  std::cout << "SHA256 sum is 0x";
  for(int i=0;i<csum.size;i++) { 
    std::cout << std::hex << csum[i];
  }
  std::cout << std::endl;
}
*/
