#ifndef AHADIC_Tools_KT_Selector_H
#define AHADIC_Tools_KT_Selector_H

#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Message.H"


namespace AHADIC {
  class Hadronisation_Reweighting;

  class KT_Selector {
  private:
    std::vector<double> m_sigma, m_ktmax;
    size_t m_n_kt_variations {0};
    Hadronisation_Reweighting * p_reweighting;
    double SelectKT(const double ktmax);
    double Gaussian(const double, const double);
    double Erf(const double, const double);
  public:
    KT_Selector(Hadronisation_Reweighting * reweighting);
    ~KT_Selector();
    void Init();

    double operator()(const double & ktmax, const bool vary_ptmax);
  };
}

#endif
