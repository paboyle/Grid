#ifndef Hadrons_MScalarSUN_Utils_hpp_
#define Hadrons_MScalarSUN_Utils_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>

BEGIN_HADRONS_NAMESPACE

BEGIN_MODULE_NAMESPACE(MScalarSUN)

GRID_SERIALIZABLE_ENUM(DiffType, undef, forward, 1, backward, 2, central, 3);

template <typename Field>
inline void dmu(Field &out, const Field &in, const unsigned int mu, const DiffType type)
{
    auto & env = Environment::getInstance();

    if (mu >= env.getNd())
    {
        HADRON_ERROR(Range, "Derivative direction out of range");
    }
    switch(type)
    {
        case DiffType::backward:
            out = in - Cshift(in, mu, -1);
            break;
        case DiffType::forward:
            out = Cshift(in, mu, 1) - in;
            break;
        case DiffType::central:
            out = 0.5*(Cshift(in, mu, 1) - Cshift(in, mu, -1));
            break;
        default:
            HADRON_ERROR(Argument, "Derivative type invalid");
            break;
    }
}

template <typename Field>
inline void dmuAcc(Field &out, const Field &in, const unsigned int mu, const DiffType type)
{
    auto & env = Environment::getInstance();

    if (mu >= env.getNd())
    {
        HADRON_ERROR(Range, "Derivative direction out of range");
    }
    switch(type)
    {
        case DiffType::backward:
            out += in - Cshift(in, mu, -1);
            break;
        case DiffType::forward:
            out += Cshift(in, mu, 1) - in;
            break;
        case DiffType::central:
            out += 0.5*(Cshift(in, mu, 1) - Cshift(in, mu, -1));
            break;
        default:
            HADRON_ERROR(Argument, "Derivative type invalid");
            break;
    }
}

inline std::string varName(const std::string name, const unsigned int mu)
{
    return name + "_" + std::to_string(mu);
}

inline std::string varName(const std::string name, const unsigned int mu, 
                           const unsigned int nu)
{
    return name + "_" + std::to_string(mu) + "_" + std::to_string(nu);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MScalarSUN_Utils_hpp_
