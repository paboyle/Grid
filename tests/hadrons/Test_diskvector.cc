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

    SerializableDiskVector<Object, TestReader, TestWriter> v("diskvector_test", 1000, 2);

    obj.e = Enum::red;
    random(rng, obj.scm);
    v[32] = obj;
    random(rng, obj.scm);
    v[2] = obj;
    v2w  = obj;
    random(rng, obj.scm);
    v[6] = obj;
    random(rng, obj.scm);
    v[12] = obj;
    random(rng, obj.scm);
    v[13] = obj;
    v13w  = obj;

    v2r = v[2];
    LOG(Message) << "v[2] correct? " 
                 << ((v2r == v2w) ? "yes" : "no" ) << std::endl;
    v13r = v[13];
    LOG(Message) << "v[13] correct? " 
                 << ((v13r == v13w) ? "yes" : "no" ) << std::endl;

    Grid_finalize();
    
    return EXIT_SUCCESS;
}
