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
    static ATOOLS::Cluster_Amplitude* CreateAmplitude(Process_Base&);

    /// calculate effective electroweak Casimir operator, cf. eq. (B.10)
    static double Cew(const Flavour& flav, int pol);

    /// calculate coefficients for L(s) terms for each m_spinampls entry
    double CalculateDoubleLogarithmicCoefficient() const;
  };
}
#endif
