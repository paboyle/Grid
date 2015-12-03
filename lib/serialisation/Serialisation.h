#ifndef GRID_SERIALISATION_READER_H
#define GRID_SERIALISATION_READER_H

#include <serialisation/MacroMagic.h>
#include <serialisation/BaseIO.h>
#include <stdint.h>

//////////////////////////////////////////
// Todo:
//////////////////////////////////////////
#include <serialisation/BinaryIO.h>
#include <serialisation/TextIO.h>
//#include <serialisation/JsonIO.h>
//#include <serialisation/YamlIO.h>
#include <serialisation/XmlIO.h>

//////////////////////////////////////////
// Select the default serialiser use ifdef's
//////////////////////////////////////////
namespace Grid {
  typedef XmlReader DefaultReader;
  typedef XmlWriter DefaultWriter;
}
#endif
