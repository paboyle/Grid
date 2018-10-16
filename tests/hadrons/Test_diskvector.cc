/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Tests/Hadrons/Test_diskvector.cc

Copyright (C) 2015-2018


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
#define DV_DEBUG
#include <Hadrons/DiskVector.hpp>

using namespace Grid;
using namespace Grid::QCD;
using namespace Grid::Hadrons;

GRID_SERIALIZABLE_ENUM(Enum, undef, red, 1, blue, 2, green, 3);

class Object: Serializable {
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(Object,
                                  Enum, e,
                                  SpinColourMatrix, scm);
};

#ifdef HAVE_HDF5
typedef Hdf5Reader TestReader;
typedef Hdf5Writer TestWriter;
#else
typedef BinaryReader TestReader;
typedef BinaryWriter TestWriter;
#endif

int main(int argc, char *argv[])
{
    Grid_init(&argc, &argv);

    GridSerialRNG rng;
    Object        obj, v2w, v2r, v13w, v13r;

    SerializableDiskVector<Object, TestReader, TestWriter> v("diskvector_test", 1000, 4);

    obj.e = Enum::red;
    random(rng, obj.scm);
    v[32] = obj;
    random(rng, obj.scm);
    v[2] = obj;
    v2w  = obj;
    random(rng, obj.scm);
    v[6] = obj;
    random(rng, obj.scm);
    v[7] = obj;
    random(rng, obj.scm);
    v[8] = obj;
    random(rng, obj.scm);
    v[9] = obj;
    random(rng, obj.scm);
    v[10] = obj;
    random(rng, obj.scm);
    v[11] = obj;
    random(rng, obj.scm);
    v[12] = obj;
    random(rng, obj.scm);
    v[13] = obj;
    v13w  = obj;
    random(rng, obj.scm);
    v[14] = obj;
    random(rng, obj.scm);
    v[15] = obj;

    v2r = v[2];
    LOG(Message) << "v[2] correct? " 
                 << ((v2r == v2w) ? "yes" : "no" ) << std::endl;
    v13r = v[13];
    LOG(Message) << "v[13] correct? " 
                 << ((v13r == v13w) ? "yes" : "no" ) << std::endl;
    LOG(Message) << "hit ratio " << v.hitRatio() << std::endl;

    EigenDiskVector<ComplexD>         w("eigendiskvector_test", 1000, 4);
    EigenDiskVector<ComplexD>::Matrix m,n;

    w[2] = EigenDiskVectorMat<ComplexD>::Random(2000, 2000);
    m    = w[2];
    w[3] = EigenDiskVectorMat<ComplexD>::Random(2000, 2000);
    w[4] = EigenDiskVectorMat<ComplexD>::Random(2000, 2000);
    w[5] = EigenDiskVectorMat<ComplexD>::Random(2000, 2000);
    w[6] = EigenDiskVectorMat<ComplexD>::Random(2000, 2000);
    w[7] = EigenDiskVectorMat<ComplexD>::Random(2000, 2000);
    n    = w[2];
    LOG(Message) << "w[2] correct? " 
                 << ((m == n) ? "yes" : "no" ) << std::endl;
    LOG(Message) << "hit ratio " << w.hitRatio() << std::endl;

    Grid_finalize();
    
    return EXIT_SUCCESS;
}
