#ifndef PHASIC__EWSudakov__Sudakov_H
#define PHASIC__EWSudakov__Sudakov_H

#include "PHASIC++/EWSudakov/Comix_Interface.H"
#include "PHASIC++/EWSudakov/EWSudakov_Amplitudes.H"

#include "PHASIC++/Process/Process_Base.H"

namespace PHASIC {
  class Sudakov {
  public:
    Sudakov(Process_Base&);
    double EWSudakov(const ATOOLS::Vec4D_Vector& mom);
  private:

    /// calculate coeffs for each m_spinampls entry
    void CalculateSpinAmplitudeCoeffs();

    /// calculate double-logarithmic (DL) coefficient, \sigma_DL * L(s)
    Complex DoubleLogCoeff(const Spin_Amplitudes&, size_t spinidx);

    /// calculate effective electroweak Casimir operator, cf. eq. (B.10, B.25)
    double DiagonalCew(const Flavour& flav, int pol) const;
    double NondiagonalCew() const;  // only for neutral gauge bosons, eq. (B.25)

    // check coefficients against the literature
    bool CheckCoeffs();
    using HelicityCoeffMap = std::map<std::vector<int>, double>;
    const HelicityCoeffMap& ReferenceCoeffs();

    Process_Base& m_proc;
    EWSudakov_Amplitudes m_ampls;
    const Comix_Interface m_comixinterface;

    const double m_sw2, m_cw2, m_sw, m_cw;
    const size_t m_check;

    // intermediate results
    mutable std::vector<Complex> m_coeffs;
    mutable std::map<size_t, std::vector<Spin_Amplitudes>> m_SU2rotatedspinampls;
    mutable std::vector<Spin_Amplitudes> m_spinampls;
  };
}

#endif
