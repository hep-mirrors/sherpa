#ifndef PHASIC__EWSudakov__Sudakov_H
#define PHASIC__EWSudakov__Sudakov_H

#include "ATOOLS/Math/Histogram.H"
#include "PHASIC++/EWSudakov/HE_Comix_Interface.H"
#include "PHASIC++/EWSudakov/EWGroupConstants.H"
#include "PHASIC++/EWSudakov/EWSudakov_Amplitudes.H"
#include "MODEL/Main/Running_AlphaQED.H"

#include "PHASIC++/Process/Process_Base.H"

namespace PHASIC {

  using Coeff_Map_Key = std::pair<EWSudakov_Log_Type, Leg_Index_Set>;
  using Coeff_Value = Complex;
  using Coeff_Map_Value = std::vector<Coeff_Value>;
  using Coeff_Map = std::map<Coeff_Map_Key, Coeff_Map_Value>;

  std::ostream& operator<<(std::ostream& os, const Coeff_Map_Key&);

  class Sudakov {
  public:
    Sudakov(Process_Base*);
    ~Sudakov();
    double KFactor(const ATOOLS::Vec4D_Vector& mom);
  private:

    bool IsInHighEnergyLimit();
    void ClearSpinAmplitudes();
    void FillBaseSpinAmplitudes();

    /// calculate coeffs for each m_spinampls entry
    void CalculateSpinAmplitudeCoeffs();

    /// helper method to update values needed for replacing Z_L/W_L -> chi/phi
    void UpdateGolstoneSpincombinationAndMEPrefactor();

    /// calculate LSC L(s) coefficient
    Coeff_Value LsCoeff();

    /// calculate Z contribution to LSC l(s) coefficient
    Coeff_Value lsZCoeff();

    /// calculate SSC l(s) coefficients
    Coeff_Value lsLogROverSCoeffs(const Two_Leg_Indizes&);

    /// calculate l_C(s) coefficients
    Coeff_Value lsCCoeff();

    /// calculate l_Yuk(s) coefficients
    Coeff_Value lsYukCoeff();

    /// calculate l_PR(s) coefficients
    Coeff_Value lsPRCoeff();

    double KFactor();

    /// calculate transformed amplitude value
    Complex TransformedAmplitudeValue(const Leg_Kfcode_Map&,
                                      const std::vector<int>&,
                                      const Comix_Interface* interface=nullptr);

    Process_Base* p_proc;
    // const std::set<EWSudakov_Log_Type> m_activecoeffs;
    std::set<EWSudakov_Log_Type> m_activecoeffs;
    EWSudakov_Amplitudes m_ampls;

    /// the default COMIX interface used for all logs except the PR logs
    const Comix_Interface m_comixinterface;

    /// the COMIX interface for using dynamic EW parameters used for the PR logs
    HE_Comix_Interface m_comixinterface_he;

    const EWGroupConstants m_ewgroupconsts;

    // const double m_mw2;
    double m_threshold;
    //    MODEL::Running_AlphaQED m_runaqed;
    size_t m_check;
    Complex m_born;

    // intermediate results
    mutable Coeff_Map m_coeffs;
    mutable std::vector<Spin_Amplitudes> m_spinampls;
    mutable std::map<Leg_Kfcode_Map, std::vector<Spin_Amplitudes>> m_transformedspinampls;
    mutable Complex m_current_me_value;
    mutable std::vector<int> m_current_spincombination;
    mutable std::vector<int> m_current_goldstone_spincombination;
    mutable Complex m_current_goldstone_me_prefactor;

    // self-monitoring
    static Histogram m_kfachisto;
    static size_t m_numonshellwarning;
  };
}

#endif
