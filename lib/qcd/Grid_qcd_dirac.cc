#include <Grid.h>

namespace Grid {

  namespace QCD {

    Gamma::GammaMatrix  Gamma::GammaMatrices [] = {
      Gamma::Identity,
      Gamma::GammaX,
      Gamma::GammaY,
      Gamma::GammaZ,
      Gamma::GammaT,
      Gamma::Gamma5,
      Gamma::MinusIdentity,
      Gamma::MinusGammaX,
      Gamma::MinusGammaY,
      Gamma::MinusGammaZ,
      Gamma::MinusGammaT,
      Gamma::MinusGamma5
    };
    const char *Gamma::GammaMatrixNames[] = { 
      "Identity ",
      "GammaX   ",
      "GammaY   ",
      "GammaZ   ",
      "GammaT   ",
      "Gamma5   ",
      "-Identity",
      "-GammaX  ",
      "-GammaY  ",
      "-GammaZ  ",
      "-GammaT  ",
      "-Gamma5  ",
      "         "
    };


  }
}
