#ifndef AHADIC_Tools_KT_Selector_H
#define AHADIC_Tools_KT_Selector_H

#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Message.H"


namespace AHADIC {
  class Ahadic_Reweighting;

  class KT_Selector {
  private:
    std::vector<double> m_sigma;
    Ahadic_Reweighting * p_reweighting;
    double WeightFunction(const double & kt);
    double SelectKT(const double ktmax);
    double Gaussian(const double, const double);
    double Erf(const double, const double);
    void KTAccepted(const double kt);
    void KTRejected(const double kt);
  public:
    KT_Selector(Ahadic_Reweighting * reweighting);
    ~KT_Selector();
    void Init();

    double operator()(const double & ktmax);
  };
}

#endif
