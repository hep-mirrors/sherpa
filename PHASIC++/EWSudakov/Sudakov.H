#ifndef PHASIC__EWSudakov__Sudakov_H
#define PHASIC__EWSudakov__Sudakov_H

#include "PHASIC++/EWSudakov/Comix_Interface.H"
#include "PHASIC++/EWSudakov/EWSudakov_Amplitudes.H"

#include "PHASIC++/Process/Process_Base.H"

namespace PHASIC {
  class Sudakov {
  public:
    Sudakov(Process_Base*);
    double EWSudakov(const ATOOLS::Vec4D_Vector& mom);
  private:

    /// calculate coeffs for each m_spinampls entry
    void CalculateSpinAmplitudeCoeffs();

    /// calculate LSC L(s) coefficient
    Complex LsCoeff(Complex amplvalue,
                    std::vector<int> spincombination,
                    size_t spinidx);

    /// calculate Z contribution to LSC l(s) coefficient
    Complex lsZCoeff(Complex amplvalue,
                     std::vector<int> spincombination,
                     size_t spinidx);

    /// calculate SSC l(s) coefficients
    std::map<std::pair<size_t, size_t>, Complex> lsLogROverSCoeffs(
        Complex amplvalue,
        std::vector<int> spincombination,
        size_t spinidx);

    /// calculate effective electroweak Casimir operator, cf. eq. (B.10, B.25)
    double DiagonalCew(const Flavour&, int pol) const;

    /// the same, but non-diagonal neutral gauge bosons terms, eq. (B.25)
    double NondiagonalCew() const;

    /// calculate the squared coupling to Z
    double IZ2(const Flavour&, int pol) const;

    /// calculate the coupling to Z
    double IZ(const Flavour&, int pol) const;

    /// calculate the coupling to W plus/minus
    double Ipm(const Flavour&, int pol, const std::string& sign) const;

    /// combine m_coeff with logs to give delta
    Complex deltaEW(const double s);

    Process_Base* p_proc;
    EWSudakov_Amplitudes m_ampls;
    const Comix_Interface m_comixinterface;

    const double m_sw2, m_cw2, m_sw, m_cw,
      m_mw2, m_mz2;
    const size_t m_check;
    Complex m_born;

    // intermediate results
    mutable std::map<std::string, std::vector<Complex>> m_coeffs;
    mutable std::map<size_t, std::vector<Spin_Amplitudes>> m_SU2rotatedspinampls;
    mutable std::vector<Spin_Amplitudes> m_spinampls;
  };
}

#endif
