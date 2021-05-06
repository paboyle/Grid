#include <Grid/Grid.h>

NAMESPACE_BEGIN(Grid);

const std::array<const GparityFlavour, 3> GparityFlavour::sigma_mu = {{
    GparityFlavour(GparityFlavour::Algebra::SigmaX),
    GparityFlavour(GparityFlavour::Algebra::SigmaY),
    GparityFlavour(GparityFlavour::Algebra::SigmaZ)
    }};

const std::array<const GparityFlavour, 6> GparityFlavour::sigma_all = {{
  GparityFlavour(GparityFlavour::Algebra::Identity),
  GparityFlavour(GparityFlavour::Algebra::SigmaX),
  GparityFlavour(GparityFlavour::Algebra::SigmaY),
  GparityFlavour(GparityFlavour::Algebra::SigmaZ),
  GparityFlavour(GparityFlavour::Algebra::ProjPlus),
  GparityFlavour(GparityFlavour::Algebra::ProjMinus)
}};

const std::array<const char *, GparityFlavour::nSigma> GparityFlavour::name = {{
    "SigmaX",
    "MinusSigmaX",
    "SigmaY",
    "MinusSigmaY",
    "SigmaZ",
    "MinusSigmaZ",
    "Identity",
    "MinusIdentity",
    "ProjPlus",
    "MinusProjPlus",
    "ProjMinus",
    "MinusProjMinus"}};

NAMESPACE_END(Grid);
