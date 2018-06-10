#ifndef PHASIC__EWSudakov__Sudakov_H
#define PHASIC__EWSudakov__Sudakov_H

#include "PHASIC++/EWSudakov/Comix_Interface.H"
#include "PHASIC++/EWSudakov/EWGroupConstants.H"
#include "PHASIC++/EWSudakov/EWSudakov_Amplitudes.H"

#include "PHASIC++/Process/Process_Base.H"

namespace PHASIC {

  using Coeff_Map_Key = std::pair<EWSudakov_Log_Type, Leg_Index_Set>;
  using Coeff_Value = std::pair<Complex, Complex>;
  using Coeff_Map_Value = std::vector<Coeff_Value>;
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
    Coeff_Value LsCoeff(Complex amplvalue,
                        std::vector<int> spincombination);

    /// calculate Z contribution to LSC l(s) coefficient
    Coeff_Value lsZCoeff(Complex amplvalue,
                         std::vector<int> spincombination);

    /// calculate SSC l(s) coefficients
    Coeff_Value lsLogROverSCoeffs(Complex amplvalue,
                                  std::vector<int> spincombination,
                                  const Two_Leg_Indizes& indizes);

    /// combine m_coeff with logs to give delta
    Complex deltaEW(const double s);

    Process_Base* p_proc;
    const std::set<EWSudakov_Log_Type> m_activecoeffs;
    EWSudakov_Amplitudes m_ampls;
    const Comix_Interface m_comixinterface;
    const EWGroupConstants m_ewgroupconsts;

    const double m_mw2, m_mz2;
    const size_t m_check;
    Complex m_born;

    // intermediate results
    mutable Coeff_Map m_coeffs;
    mutable std::map<size_t, std::vector<Spin_Amplitudes>> m_lsczspinampls;
    mutable std::map<Leg_Set, std::vector<Spin_Amplitudes>> m_sscwspinampls;
    mutable std::vector<Spin_Amplitudes> m_spinampls;
  };
}

#endif
