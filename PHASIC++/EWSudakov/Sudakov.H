#ifndef PHASIC__EWSudakov__Sudakov_H
#define PHASIC__EWSudakov__Sudakov_H

#include "PHASIC++/EWSudakov/Comix_Interface.H"

#include "PHASIC++/Process/Process_Base.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"

namespace PHASIC {
  class Sudakov {
  public:
    Sudakov(Process_Base&);
    double EWSudakov(const ATOOLS::Vec4D_Vector& mom);
  private:
    Process_Base& m_proc;
    ATOOLS::Cluster_Amplitude *p_ampl;
    Comix_Interface m_ci;
    std::vector<Spin_Amplitudes> m_spinampls;
    std::vector<double> m_coeffs;
    double m_sw2, m_cw2, m_sw, m_cw;

    static ATOOLS::Cluster_Amplitude* CreateAmplitude(Process_Base&);
    void UpdateAmplitude(const ATOOLS::Vec4D_Vector& mom);

    /// calculate coeffs for each m_spinampls entry
    void CalculateSpinAmplitudeCoeffs();

    /// calculate double-logarithmic (DL) coefficient, \sigma_DL * L(s)
    double DoubleLogCoeff(std::vector<int> spincombination) const;

    /// calculate effective electroweak Casimir operator, cf. eq. (B.10, B.25)
    double DiagonalCew(const Flavour& flav, int pol) const;
    double NondiagonalCew(kf_code from, kf_code to) const;  // only Z/\gamma
  };
}
#endif
