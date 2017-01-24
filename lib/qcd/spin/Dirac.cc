#include <Grid.h>

namespace Grid {
namespace QCD {

#include "GammaMulTable.h"

const std::array<const char *, Gamma::nGamma> Gamma::name = {{
  "-Gamma5      ",
  "Gamma5       ",
  "-GammaT      ",
  "GammaT       ",
  "-GammaTGamma5",
  "GammaTGamma5 ",
  "-GammaX      ",
  "GammaX       ",
  "-GammaXGamma5",
  "GammaXGamma5 ",
  "-GammaY      ",
  "GammaY       ",
  "-GammaYGamma5",
  "GammaYGamma5 ",
  "-GammaZ      ",
  "GammaZ       ",
  "-GammaZGamma5",
  "GammaZGamma5 ",
  "-Identity    ",
  "Identity     ",
  "-SigmaXT     ",
  "SigmaXT      ",
  "-SigmaXY     ",
  "SigmaXY      ",
  "-SigmaXZ     ",
  "SigmaXZ      ",
  "-SigmaYT     ",
  "SigmaYT      ",
  "-SigmaYZ     ",
  "SigmaYZ      ",
  "-SigmaZT     ",
  "SigmaZT      "}};

}}
