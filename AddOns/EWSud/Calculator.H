#ifndef EWSudakov__Calculator_H
#define EWSudakov__Calculator_H

#include "ATOOLS/Math/Histogram.H"
#include "AddOns/EWSud/HE_Comix_Interface.H"
#include "AddOns/EWSud/EWGroupConstants.H"
#include "AddOns/EWSud/Amplitudes.H"
#include "MODEL/Main/Running_AlphaQED.H"

#include "PHASIC++/Process/Process_Base.H"

namespace EWSud {

  using Coeff_Map_Key = std::pair<EWSudakov_Log_Type, Leg_Index_Set>;
  using Coeff_Value = Complex;
  using Coeff_Map_Value = std::vector<Coeff_Value>;
  using Coeff_Map = std::map<Coeff_Map_Key, Coeff_Map_Value>;
  using HighEnergySchemes = EWSudakovHighEnergySchemes;

  std::ostream& operator<<(std::ostream& os, const Coeff_Map_Key&);

  class Calculator {

  public:

    Calculator(PHASIC::Process_Base*);
    ~Calculator();
    EWSudakov_Log_Corrections_Map
    CorrectionsMap(const ATOOLS::Vec4D_Vector& mom);
    const std::set<EWSudakov_Log_Type>& ActiveLogTypes() const
    { return m_activelogtypes; }
    const EWGroupConstants& GetEWGroupConstants() const {
      return m_ewgroupconsts;
    }
    Comix_Interface& GetComixInterface() {
      return m_comixinterface;
    }

  private:
    HighEnergySchemes m_helimitscheme;
    bool IsInHighEnergyLimit();
    void ClearSpinAmplitudes();
    void FillBaseSpinAmplitudes();

    /// calculate coeffs for each m_spinampls entry
    void CalculateSpinAmplitudeCoeffs();

    /// calculate the conversion factor for the GBET (deprecated, see impl.)
    Complex GBETConversionFactor(Leg_Kfcode_Map legs = {});

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

    EWSudakov_Log_Corrections_Map CorrectionsMap();
    /// Set High Energy limit scheme from string
    void SetHighEnergyScheme(const std::string&);
    /// calculate transformed amplitude value
    Complex TransformedAmplitudeValue(const Leg_Kfcode_Map&,
                                      const std::vector<int>&,
                                      const Comix_Interface* interface=nullptr);

    PHASIC::Process_Base* p_proc;
    std::set<EWSudakov_Log_Type> m_activelogtypes;
    bool m_c_coeff_ignores_vector_bosons;
    Amplitudes m_ampls;

    /// the default COMIX interface used for all logs except the PR logs
    Comix_Interface m_comixinterface;

    /// the COMIX interface for using dynamic EW parameters used for the PR logs
    HE_Comix_Interface m_comixinterface_he;

    const EWGroupConstants m_ewgroupconsts;

    double m_threshold;
    bool m_checkinvariantratios;
    bool m_checkcoeff, m_checkkfac;
    std::string m_checklogfile;
    bool m_includesubleading;

    // intermediate results
    mutable Coeff_Map m_coeffs;
    mutable std::vector<METOOLS::Spin_Amplitudes> m_spinampls;
    mutable Complex m_current_me_value;
    mutable std::vector<int> m_current_spincombination;
    mutable Complex m_current_goldstone_me_prefactor;
    mutable double m_ewscale2;

    // self-monitoring
    static ATOOLS::Histogram m_kfachisto;
    static size_t m_numonshellwarning;
  };
}

#endif
