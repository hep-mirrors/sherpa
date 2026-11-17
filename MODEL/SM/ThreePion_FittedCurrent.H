#ifndef MODEL_SM_ThreePion_FittedCurrent_H
#define MODEL_SM_ThreePion_FittedCurrent_H

#include "METOOLS/FormFactors/Propagator.H"
#include <memory>
#include <string>

namespace METOOLS {

  /*!
    THREE_PION_FORM_FACTOR: Current, which gamma* -> pi+ pi- pi0 current the
    hep-ph/0512180 eq. (10), with
      BW_V    = M^2/(M^2 - q^2 - i M Gamma), constant width, at the masses of
                the fit: omega(782) and phi(1020) at the particle-table masses
                of the fitting branch, omega' = omega(1375) and omega'' =
                omega(1631) at hep-ph/0512180's own masses and widths;
      H_k     = the rho_k line shape summed over the three pion pairs, all plus:
                the neutral rho at s_+-, the charged one at s_+0 and s_-0;
      rho(770)  Gounaris-Sakurai with the imaginary part M Gamma(s)
                (GounarisSakuraiM of the fitting branch), M = 0.77,
                Gamma = 0.1507, pi+ mass 0.13957;
      rho(1450), rho(1700)
                M^2/(M^2 - s - i sqrt(s) Gamma(s)) with the running width of
                the line-shape registry (METOOLS/FormFactors/Rho_Decays.C, the
                same code as on the fitting branch).
  */
  class Fitted_ThreePion_Current {
  public:
    //! the three products of hep-ph/0512180 eq. (10), named by their rho
    enum class Tower { rho770 = 0, rho1450 = 1, rho1700 = 2 };

    //! THREE_PION_FORM_FACTOR: Current: Fitted (true) or LowEnergy (false)
    static bool Selected();
    static Fitted_ThreePion_Current & Instance();
    ~Fitted_ThreePion_Current();

    Propagator_Base * NewIsoscalar(const Tower & k,const Complex & scale) const;
    Complex RhoLineShape(const Tower & k,const bool & charged,const double & s);
    Complex operator()(const double & q2,const double & spm,
                       const double & sp0,const double & sm0);

    std::string Describe() const;

    Fitted_ThreePion_Current(const Fitted_ThreePion_Current &) = delete;
    Fitted_ThreePion_Current & operator=(const Fitted_ThreePion_Current &) = delete;
  private:
    Fitted_ThreePion_Current();
    struct Impl;
    std::unique_ptr<Impl> p_impl;
  };

}

#endif
