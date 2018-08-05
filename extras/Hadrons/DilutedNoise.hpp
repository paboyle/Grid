#ifndef Hadrons_DilutedNoise_hpp_
#define Hadrons_DilutedNoise_hpp_

#include <Grid/Hadrons/Global.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                   Abstract container for diluted noise                     *
 ******************************************************************************/
template <typename FImpl>
class DilutedNoise
{
public:
    typedef typename FImpl::FermionField FermionField;
public:
    // constructor/destructor
    DilutedNoise(GridCartesian *g);
    DilutedNoise(GridCartesian *g, const unsigned int nNoise);
    virtual ~DilutedNoise(void) = default;
    // access
    std::vector<FermionField> &       getNoise(void);
    const std::vector<FermionField> & getNoise(void) const;
    const FermionField &              operator[](const unsigned int i) const;
    FermionField &                    operator[](const unsigned int i);
    void                              resize(const unsigned int nNoise);
    unsigned int                      size(void) const;
    GridCartesian                     *getGrid(void) const;
    // generate noise (pure virtual)
    virtual void generateNoise(GridParallelRNG &rng) = 0;
private:
    std::vector<FermionField> noise_;
    GridCartesian             *grid_;
    unsigned int              nNoise_;
};

template <typename FImpl>
class TimeDilutedSpinColorDiagonalNoise: public DilutedNoise<FImpl>
{
public:
    typedef typename FImpl::FermionField FermionField;
public:
    // constructor/destructor
    TimeDilutedSpinColorDiagonalNoise(GridCartesian *g);
    virtual ~TimeDilutedSpinColorDiagonalNoise(void) = default;
    // generate noise
    virtual void generateNoise(GridParallelRNG &rng);
private:
    unsigned int nt_;
};

/******************************************************************************
 *                    DilutedNoise template implementation                    *
 ******************************************************************************/
template <typename FImpl>
DilutedNoise<FImpl>::DilutedNoise(GridCartesian *g)
: grid_(g)
{}

template <typename FImpl>
DilutedNoise<FImpl>::DilutedNoise(GridCartesian *g,
                                  const unsigned int nNoise)
: DilutedNoise(g)
{
    resize(nNoise);
}

template <typename FImpl>
std::vector<typename DilutedNoise<FImpl>::FermionField> & DilutedNoise<FImpl>::
getNoise(void)
{
    return noise_;
}

template <typename FImpl>
const std::vector<typename DilutedNoise<FImpl>::FermionField> & DilutedNoise<FImpl>::
getNoise(void) const
{
    return noise_;
}

template <typename FImpl>
const typename DilutedNoise<FImpl>::FermionField & 
DilutedNoise<FImpl>::operator[](const unsigned int i) const
{
    return noise_[i];
}

template <typename FImpl>
typename DilutedNoise<FImpl>::FermionField & 
DilutedNoise<FImpl>::operator[](const unsigned int i)
{
    return noise_[i];
}

template <typename FImpl>
void DilutedNoise<FImpl>::resize(const unsigned int nNoise)
{
    nNoise_ = nNoise;
    noise_.resize(nNoise, grid_);
}

template <typename FImpl>
unsigned int DilutedNoise<FImpl>::size(void) const
{  
    return noise_.size();
}

template <typename FImpl>
GridCartesian * DilutedNoise<FImpl>::getGrid(void) const
{
    return grid_;
}

/******************************************************************************
 *        TimeDilutedSpinColorDiagonalNoise template implementation           *
 ******************************************************************************/
template <typename FImpl>
TimeDilutedSpinColorDiagonalNoise<FImpl>::
TimeDilutedSpinColorDiagonalNoise(GridCartesian *g)
: DilutedNoise<FImpl>(g)
{
    nt_ = this->getGrid()->GlobalDimensions().back();
    this->resize(nt_*Ns*FImpl::Dimension);
}

template <typename FImpl>
void TimeDilutedSpinColorDiagonalNoise<FImpl>::generateNoise(GridParallelRNG &rng)
{
    typedef decltype(peekColour((*this)[0], 0)) SpinField;

    auto                       &noise = *this;
    auto                       g      = this->getGrid();
    auto                       nd     = g->GlobalDimensions().size();
    auto                       nc     = FImpl::Dimension;
    Complex                    shift(1., 1.);
    Lattice<iScalar<vInteger>> tLat(g);
    LatticeComplex             eta(g), etaCut(g);
    SpinField                  etas(g);
    unsigned int               i = 0;

    LatticeCoordinate(tLat, nd - 1);
    bernoulli(rng, eta);
    eta = (2.*eta - shift)*(1./::sqrt(2.));
    std::cout << eta << std::endl;
    for (unsigned int t = 0; t < nt_; ++t)
    {
        etaCut = where((tLat == t), eta, 0.*eta);
        for (unsigned int s = 0; s < Ns; ++s)
        {
            etas = zero;
            pokeSpin(etas, etaCut, s);
            for (unsigned int c = 0; c < nc; ++c)
            {
                noise[i] = zero;
                pokeColour(noise[i], etas, c);
                std::cout << noise[i] << std::endl;
                i++;
            }
        }
    }
}

END_HADRONS_NAMESPACE

#endif // Hadrons_DilutedNoise_hpp_