#ifndef BEAM_Main_Collider_Weight_H
#define BEAM_Main_Collider_Weight_H

#include "BEAM/Main/Beam_Rejection.H"
#include "BEAM/Main/Weight_Base.H"

namespace BEAM {
  class Collider_Weight: public Weight_Base {
  private:

    collidermode m_mode;
    ATOOLS::Info_Key m_sprimekey, m_ykey, m_xkey;

    int m_rejection;
    Beam_Rejection_Base* p_rejector;
    // random number for azimuth between the two impact parameters
    double m_eran;

    // BEAM_SPREAD_CORRELATION: correlated Gaussian beam energy spread.
    // p_gauss[i] are the two Gaussian spectra, m_rho
    // the correlation and m_rhonorm the normalisation of the correlated
    // density on the truncated box (see CorrelationFactor).
    const class Gaussian* p_gauss[2] = {nullptr, nullptr};
    double m_rho = 0., m_rhonorm = 1.;
    double CorrelationFactor();

  public:

    static double CorrelationWeight(double d1, double d2, double rho,
                                    double rhonorm);
    // Probability of the correlated standard Gaussian in |d1| < n1, |d2| < n2.
    static double BoxProbability(double rho, double n1, double n2);

    Collider_Weight(Kinematics_Base* kinematics);
    ~Collider_Weight();

    void AssignKeys(ATOOLS::Integration_Info* const info) override;
    bool Calculate(const double& scale) override;
    double operator()() override;
    void SetERan(double eran) override { m_eran = eran; }
    double OverlapWeight();
  };
} // namespace BEAM
#endif
