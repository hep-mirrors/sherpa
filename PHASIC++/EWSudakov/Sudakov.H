#ifndef PHASIC__EWSudakov__Sudakov_H
#define PHASIC__EWSudakov__Sudakov_H

#include "PHASIC++/EWSudakov/Comix_Interface.H"
#include "PHASIC++/EWSudakov/EWSudakov_Amplitudes.H"

#include "PHASIC++/Process/Process_Base.H"

namespace PHASIC {

  using Coeff_Map_Key = std::pair<EWSudakov_Log_Type, Leg_Index_Set>;
  using Coeff_Map_Value = std::vector<Complex>;
  using Coeff_Map = std::map<Coeff_Map_Key, Coeff_Map_Value>;

  std::ostream& operator<<(std::ostream& os, const Coeff_Map_Key&);

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
    Complex lsLogROverSCoeffs(Complex amplvalue,
                              std::vector<int> spincombination,
                              size_t spinidx,
                              const Two_Leg_Indizes& indizes);

    /// calculate effective electroweak Casimir operator, cf. eq. (B.10, B.25)
    double DiagonalCew(const Flavour&, int pol) const;

    /// the same, but non-diagonal neutral gauge bosons terms, eq. (B.25)
    double NondiagonalCew() const;

    /// calculate the squared coupling to Z
    double IZ2(const Flavour&, int pol) const;

    /// calculate the coupling to Z
    double IZ(const Flavour&, int pol) const;

    /// calculate the coupling to W plus/minus
    double Ipm(const Flavour&, int pol, bool isplus) const;

    /// combine m_coeff with logs to give delta
    Complex deltaEW(const double s);

    Process_Base* p_proc;
    const std::set<EWSudakov_Log_Type> m_activecoeffs;
    EWSudakov_Amplitudes m_ampls;
    const Comix_Interface m_comixinterface;

    const double m_sw2, m_cw2, m_sw, m_cw,
      m_mw2, m_mz2;
    const size_t m_check;
    Complex m_born;

    // intermediate results
    mutable Coeff_Map m_coeffs;
    mutable std::map<size_t, std::vector<Spin_Amplitudes>> m_lsczspinampls;
    mutable std::map<Two_Leg_Indizes, std::vector<Spin_Amplitudes>> m_sscwspinampls;
    mutable std::vector<Spin_Amplitudes> m_spinampls;
  };
}

#endif
