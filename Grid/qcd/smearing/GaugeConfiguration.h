/*!
  @file GaugeConfiguration.h

  @brief Declares the GaugeConfiguration class
*/
#pragma once

NAMESPACE_BEGIN(Grid);


//trivial class for no smearing
template< class Impl >
class NoSmearing : public ConfigurationBase<typename Impl::Field>
{
public:
  INHERIT_FIELD_TYPES(Impl);

  Field* ThinLinks;

  NoSmearing(): ThinLinks(NULL) {}

  virtual void set_Field(Field& U) { ThinLinks = &U; }

  virtual void smeared_force(Field&) {}

  virtual Field& get_SmearedU() { return *ThinLinks; }

  virtual Field &get_U(bool smeared = false)
  {
    return *ThinLinks;
  }
};

/*!
  @brief Smeared configuration container

  It will behave like a configuration from the point of view of
  the HMC update and integrators.
  An "advanced configuration" object that can provide not only the
  data to store the gauge configuration but also operations to manipulate
  it, like smearing.

  It stores a list of smeared configurations.
*/
template <class Gimpl>
class SmearedConfiguration : public ConfigurationBase<typename Gimpl::Field>
{
public:
  INHERIT_GIMPL_TYPES(Gimpl);

protected:
  const unsigned int smearingLevels;
  Smear_Stout<Gimpl> *StoutSmearing;
  std::vector<GaugeField> SmearedSet;
public:
  GaugeField*  ThinLinks; /* Pointer to the thin links configuration */ // move to base???
protected:
  
  // Member functions
  //====================================================================

  // Overridden in masked version
  virtual void fill_smearedSet(GaugeField &U)
  {
    ThinLinks = &U;  // attach the smearing routine to the field U

    // check the pointer is not null
    if (ThinLinks == NULL)
      std::cout << GridLogError
                << "[SmearedConfiguration] Error in ThinLinks pointer\n";

    if (smearingLevels > 0)
    {
      std::cout << GridLogDebug
                << "[SmearedConfiguration] Filling SmearedSet\n";
      GaugeField previous_u(ThinLinks->Grid());

      previous_u = *ThinLinks;
      for (int smearLvl = 0; smearLvl < smearingLevels; ++smearLvl)
      {
        StoutSmearing->smear(SmearedSet[smearLvl], previous_u);
        previous_u = SmearedSet[smearLvl];

        // For debug purposes
        RealD impl_plaq = WilsonLoops<Gimpl>::avgPlaquette(previous_u);
        std::cout << GridLogDebug
                  << "[SmearedConfiguration] Plaq: " << impl_plaq << std::endl;
      }
    }
  }

  //overridden in masked verson
  virtual GaugeField AnalyticSmearedForce(const GaugeField& SigmaKPrime,
					  const GaugeField& GaugeK) const 
  {
    GridBase* grid = GaugeK.Grid();
    GaugeField C(grid), SigmaK(grid), iLambda(grid);
    GaugeLinkField iLambda_mu(grid);
    GaugeLinkField iQ(grid), e_iQ(grid);
    GaugeLinkField SigmaKPrime_mu(grid);
    GaugeLinkField GaugeKmu(grid), Cmu(grid);

    StoutSmearing->BaseSmear(C, GaugeK);
    SigmaK = Zero();
    iLambda = Zero();

    for (int mu = 0; mu < Nd; mu++)
    {
      Cmu = peekLorentz(C, mu);
      GaugeKmu = peekLorentz(GaugeK, mu);
      SigmaKPrime_mu = peekLorentz(SigmaKPrime, mu);
      iQ = Ta(Cmu * adj(GaugeKmu));
      set_iLambda(iLambda_mu, e_iQ, iQ, SigmaKPrime_mu, GaugeKmu);
      pokeLorentz(SigmaK, SigmaKPrime_mu * e_iQ + adj(Cmu) * iLambda_mu, mu);
      pokeLorentz(iLambda, iLambda_mu, mu);
    }
    StoutSmearing->derivative(SigmaK, iLambda,
                             GaugeK);  // derivative of SmearBase
    return SigmaK;
  }

  /*! @brief Returns smeared configuration at level 'Level' */
public:
  const GaugeField &get_smeared_conf(int Level) const
  {
    return SmearedSet[Level];
  }

  //====================================================================
  void set_iLambda(GaugeLinkField& iLambda, GaugeLinkField& e_iQ,
                   const GaugeLinkField& iQ, const GaugeLinkField& Sigmap,
                   const GaugeLinkField& GaugeK) const 
  {
    GridBase* grid = iQ.Grid();
    GaugeLinkField iQ2(grid), iQ3(grid), B1(grid), B2(grid), USigmap(grid);
    GaugeLinkField unity(grid);
    unity = 1.0;

    LatticeComplex u(grid), w(grid);
    LatticeComplex f0(grid), f1(grid), f2(grid);
    LatticeComplex xi0(grid), xi1(grid), tmp(grid);
    LatticeComplex u2(grid), w2(grid), cosw(grid);
    LatticeComplex emiu(grid), e2iu(grid), qt(grid), fden(grid);
    LatticeComplex r01(grid), r11(grid), r21(grid), r02(grid), r12(grid);
    LatticeComplex r22(grid), tr1(grid), tr2(grid);
    LatticeComplex b10(grid), b11(grid), b12(grid), b20(grid), b21(grid),
      b22(grid);
    LatticeComplex LatticeUnitComplex(grid);

    LatticeUnitComplex = 1.0;

    // Exponential
    iQ2 = iQ * iQ;
    iQ3 = iQ * iQ2;
    StoutSmearing->set_uw(u, w, iQ2, iQ3);
    StoutSmearing->set_fj(f0, f1, f2, u, w);
    e_iQ = f0 * unity + timesMinusI(f1) * iQ - f2 * iQ2;

    // Getting B1, B2, Gamma and Lambda
    // simplify this part, reduntant calculations in set_fj
    xi0 = StoutSmearing->func_xi0(w);
    xi1 = StoutSmearing->func_xi1(w);
    u2 = u * u;
    w2 = w * w;
    cosw = cos(w);

    emiu = cos(u) - timesI(sin(u));
    e2iu = cos(2.0 * u) + timesI(sin(2.0 * u));

    r01 = (2.0 * u + timesI(2.0 * (u2 - w2))) * e2iu +
      emiu * ((16.0 * u * cosw + 2.0 * u * (3.0 * u2 + w2) * xi0) +
	      timesI(-8.0 * u2 * cosw + 2.0 * (9.0 * u2 + w2) * xi0));

    r11 = (2.0 * LatticeUnitComplex + timesI(4.0 * u)) * e2iu +
      emiu * ((-2.0 * cosw + (3.0 * u2 - w2) * xi0) +
	      timesI((2.0 * u * cosw + 6.0 * u * xi0)));

    r21 =
      2.0 * timesI(e2iu) + emiu * (-3.0 * u * xi0 + timesI(cosw - 3.0 * xi0));

    r02 = -2.0 * e2iu +
      emiu * (-8.0 * u2 * xi0 +
	      timesI(2.0 * u * (cosw + xi0 + 3.0 * u2 * xi1)));

    r12 = emiu * (2.0 * u * xi0 + timesI(-cosw - xi0 + 3.0 * u2 * xi1));

    r22 = emiu * (xi0 - timesI(3.0 * u * xi1));

    fden = LatticeUnitComplex / (2.0 * (9.0 * u2 - w2) * (9.0 * u2 - w2));

    b10 = 2.0 * u * r01 + (3.0 * u2 - w2) * r02 - (30.0 * u2 + 2.0 * w2) * f0;
    b11 = 2.0 * u * r11 + (3.0 * u2 - w2) * r12 - (30.0 * u2 + 2.0 * w2) * f1;
    b12 = 2.0 * u * r21 + (3.0 * u2 - w2) * r22 - (30.0 * u2 + 2.0 * w2) * f2;

    b20 = r01 - (3.0 * u) * r02 - (24.0 * u) * f0;
    b21 = r11 - (3.0 * u) * r12 - (24.0 * u) * f1;
    b22 = r21 - (3.0 * u) * r22 - (24.0 * u) * f2;

    b10 *= fden;
    b11 *= fden;
    b12 *= fden;
    b20 *= fden;
    b21 *= fden;
    b22 *= fden;

    B1 = b10 * unity + timesMinusI(b11) * iQ - b12 * iQ2;
    B2 = b20 * unity + timesMinusI(b21) * iQ - b22 * iQ2;
    USigmap = GaugeK * Sigmap;

    tr1 = trace(USigmap * B1);
    tr2 = trace(USigmap * B2);

    GaugeLinkField QUS = iQ * USigmap;
    GaugeLinkField USQ = USigmap * iQ;

    GaugeLinkField iGamma = tr1 * iQ - timesI(tr2) * iQ2 +
      timesI(f1) * USigmap + f2 * QUS + f2 * USQ;

    iLambda = Ta(iGamma);
  }

  //====================================================================
public:

  /* Standard constructor */
  SmearedConfiguration(GridCartesian* UGrid, unsigned int Nsmear,
                       Smear_Stout<Gimpl>& Stout)
      : smearingLevels(Nsmear), StoutSmearing(&Stout), ThinLinks(NULL)
  {
    for (unsigned int i = 0; i < smearingLevels; ++i)
      SmearedSet.push_back(*(new GaugeField(UGrid)));
  }

  /*! For just thin links */
  SmearedConfiguration()
    : smearingLevels(0), StoutSmearing(nullptr), SmearedSet(), ThinLinks(NULL) {}

  // attach the smeared routines to the thin links U and fill the smeared set
  virtual void set_Field(GaugeField &U)
  {
    double start = usecond();
    fill_smearedSet(U);
    double end = usecond();
    double time = (end - start)/ 1e3;
    std::cout << GridLogMessage << "Smearing in " << time << " ms" << std::endl;  
  }

  //====================================================================
  virtual void smeared_force(GaugeField &SigmaTilde) 
  {
    if (smearingLevels > 0)
    {
      double start = usecond();
      GaugeField force = SigmaTilde; // actually = U*SigmaTilde
      GaugeLinkField tmp_mu(SigmaTilde.Grid());

      for (int mu = 0; mu < Nd; mu++)
      {
        // to get just SigmaTilde
        tmp_mu = adj(peekLorentz(SmearedSet[smearingLevels - 1], mu)) * peekLorentz(force, mu);
        pokeLorentz(force, tmp_mu, mu);
      }

      for (int ismr = smearingLevels - 1; ismr > 0; --ismr)
        force = AnalyticSmearedForce(force, get_smeared_conf(ismr - 1));

      force = AnalyticSmearedForce(force, *ThinLinks);

      for (int mu = 0; mu < Nd; mu++)
      {
        tmp_mu = peekLorentz(*ThinLinks, mu) * peekLorentz(force, mu);
        pokeLorentz(SigmaTilde, tmp_mu, mu);
      }
      double end = usecond();
      double time = (end - start)/ 1e3;
      std::cout << GridLogMessage << " GaugeConfiguration: Smeared Force chain rule took " << time << " ms" << std::endl;
    }  // if smearingLevels = 0 do nothing
    SigmaTilde=Gimpl::projectForce(SigmaTilde); // Ta
      
  }
  //====================================================================

  virtual GaugeField& get_SmearedU() { return SmearedSet[smearingLevels - 1]; }

  virtual GaugeField &get_U(bool smeared = false)
  {
    // get the config, thin links by default
    if (smeared)
    {
      if (smearingLevels)
      {
        RealD impl_plaq =
	  WilsonLoops<Gimpl>::avgPlaquette(SmearedSet[smearingLevels - 1]);
        std::cout << GridLogDebug << "getting Usmr Plaq: " << impl_plaq
                  << std::endl;
        return get_SmearedU();
      }
      else
      {
        RealD impl_plaq = WilsonLoops<Gimpl>::avgPlaquette(*ThinLinks);
        std::cout << GridLogDebug << "getting Thin Plaq: " << impl_plaq
                  << std::endl;
        return *ThinLinks;
      }
    }
    else
    {
      RealD impl_plaq = WilsonLoops<Gimpl>::avgPlaquette(*ThinLinks);
      std::cout << GridLogDebug << "getting Thin Plaq: " << impl_plaq
                << std::endl;
      return *ThinLinks;
    }
  }
};

NAMESPACE_END(Grid);

