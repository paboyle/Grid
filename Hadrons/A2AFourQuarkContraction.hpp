#ifndef A2A_Four_Quark_Contraction_hpp_
#define A2A_Four_Quark_Contraction_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Environment.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *           Class compute four quark contractions using disk vectors         *
 ******************************************************************************/
template <typename FImpl>
class A2AUtilsDV
{
  public:
    typedef typename FImpl::ComplexField ComplexField;
    typedef typename FImpl::FermionField FermionField;
    typedef typename FImpl::PropagatorField PropagatorField;

    typedef typename FImpl::SiteSpinor vobj;
    typedef typename vobj::scalar_object sobj;
    typedef typename vobj::scalar_type scalar_type;
    typedef typename vobj::vector_type vector_type;

    typedef iSpinMatrix<vector_type> SpinMatrix_v;
    typedef iSpinMatrix<scalar_type> SpinMatrix_s;
    typedef iSinglet<vector_type> Scalar_v;
    typedef iSinglet<scalar_type> Scalar_s;

    typedef iSpinColourMatrix<vector_type> SpinColourMatrix_v;

    template <typename TensorType>
    static void ContractWWVVDiskVector(std::vector<PropagatorField> &WWVV,
                                        const TensorType &WW_sd,
                                        const FermionField *vs,
                                        const FermionField *vd);
};
template <class FImpl>
template <typename TensorType>
void A2AUtilsDV<FImpl>::ContractWWVVDiskVector(std::vector<PropagatorField> &WWVV,
                                                const TensorType &WW_sd,
                                                const FermionField *vs,
                                                const FermionField *vd)
{
    GridBase *grid = vs[0]._grid;

    int nd = grid->_ndimension;
    int Nsimd = grid->Nsimd();
    int N_t = WW_sd.dimensions()[0];
    int N_s = WW_sd.dimensions()[1];
    int N_d = WW_sd.dimensions()[2];

    int d_unroll = 32; // Empirical optimisation

    for (int t = 0; t < N_t; t++)
    {
        std::cout << "Contraction t = " << t << std::endl;
        EigenDiskVector<ComplexD>::Matrix buf = WW_sd[t];
        parallel_for(int ss = 0; ss < grid->oSites(); ss++)
        {
            for (int d_o = 0; d_o < N_d; d_o += d_unroll)
            {
                for (int s = 0; s < N_s; s++)
                {
                    auto tmp1 = vs[s]._odata[ss];
                    vobj tmp2 = tmp1 * 0.0; // would rather write zero
                    vobj tmp3 = tmp1 * 0.0; // would rather write zero

                    for (int d = d_o; d < MIN(d_o + d_unroll, N_d); d++)
                    {
                        Scalar_v coeff = buf(s, d);
                        tmp3 = conjugate(vd[d]._odata[ss]);
                        mac(&tmp2, &coeff, &tmp3);
                    }

                    //////////////////////////
                    // Fast outer product of tmp1 with a sum of terms suppressed by d_unroll
                    //////////////////////////
                    for (int s1 = 0; s1 < Ns; s1++)
                    {
                        for (int s2 = 0; s2 < Ns; s2++)
                        {
                            WWVV[t]._odata[ss]()(s1, s2)(0, 0) += tmp1()(s1)(0) * tmp2()(s2)(0);
                            WWVV[t]._odata[ss]()(s1, s2)(0, 1) += tmp1()(s1)(0) * tmp2()(s2)(1);
                            WWVV[t]._odata[ss]()(s1, s2)(0, 2) += tmp1()(s1)(0) * tmp2()(s2)(2);
                            WWVV[t]._odata[ss]()(s1, s2)(1, 0) += tmp1()(s1)(1) * tmp2()(s2)(0);
                            WWVV[t]._odata[ss]()(s1, s2)(1, 1) += tmp1()(s1)(1) * tmp2()(s2)(1);
                            WWVV[t]._odata[ss]()(s1, s2)(1, 2) += tmp1()(s1)(1) * tmp2()(s2)(2);
                            WWVV[t]._odata[ss]()(s1, s2)(2, 0) += tmp1()(s1)(2) * tmp2()(s2)(0);
                            WWVV[t]._odata[ss]()(s1, s2)(2, 1) += tmp1()(s1)(2) * tmp2()(s2)(1);
                            WWVV[t]._odata[ss]()(s1, s2)(2, 2) += tmp1()(s1)(2) * tmp2()(s2)(2);
                        }
                    }
                }
            }
        }
    }
}

END_HADRONS_NAMESPACE

#endif // A2A_Four_Quark_Contraction_hpp_