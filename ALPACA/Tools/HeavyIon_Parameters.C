#include "ALPACA/Tools/HeavyIon_Parameters.H" 

#include "ATOOLS/Org/Run_Parameter.H" 
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Org/Message.H"

using namespace ALPACA;

HeavyIon_Parameters ALPACA::HIPars;

HeavyIon_Parameters::HeavyIon_Parameters() {}

HeavyIon_Parameters::~HeavyIon_Parameters() {}

void HeavyIon_Parameters::Init() {
    RegisterDefaults();
    ReadParameters();
}

void HeavyIon_Parameters::RegisterDefaults() const {
    ATOOLS::Settings& s = ATOOLS::Settings::GetMainSettings();

    s["Lambda"].SetDefault(1.);
    s["Xsec_Form"].SetDefault("BlackDisc");
    s["Tau_Max"].SetDefault(4.);
    s["Fixed_Sigma"].SetDefault("False");
    s["Fixed_Sigma_Value"].SetDefault(1.);
    s["Fixed_Gamma"].SetDefault("False");
    s["Fixed_Gamma_Value"].SetDefault(1.);
    s["Only_Gluons"].SetDefault("False");
    s["Splitting_Merging"].SetDefault("True");
    s["Show_Event_Information"].SetDefault("True");
    s["Show_Progress_Bar"].SetDefault("True");
    s["Tau_Restart"].SetDefault(5000.);
    s["tsample_Min"].SetDefault(0.);
    s["tsample_Max"].SetDefault(4.);
    s["Include_Bose_Factors"].SetDefault("True");
    s["p_Min"].SetDefault(0.05);
    s["N_Include"].SetDefault(5);
    s["Elastic_Scattering"].SetDefault("True");
    s["Do_Kinematics"].SetDefault("True");
    s["Tau_Max_Scaling_Limit"].SetDefault(1.);
    s["Alpha_S"].SetDefault(0.2652);
    s["Gaussian_kT2"].SetDefault("False");
    s["kT2_Reg"].SetDefault(0.0001);
    s["f_r"].SetDefault(2.);
    s["f_Delta_p"].SetDefault(0.3);
    s["f_Shell"].SetDefault("True");
    s["f_N_Max"].SetDefault(5);
    s["Formation_Time"].SetDefault(1.);
    s["Timekeeper"].SetDefault(1);
    s["OE_Mult_Scatter"].SetDefault(1.);
    s["OE_Mult_Merge"].SetDefault(1.);
    s["OE_Mult_Split"].SetDefault(1.);
    s["Test_Double"].SetDefault(0.);
    s["Test_Bool"].SetDefault("False");
    s["m2_Min_Scale"].SetDefault(1.);
    
}

void HeavyIon_Parameters::ReadParameters() {
    ATOOLS::Settings& s = ATOOLS::Settings::GetMainSettings();

    std::string xsecform = s["Xsec_Form"].Get<std::string>();
    if (xsecform==std::string("BlackDisc")){
      m_xsecform = xsec_form::BlackDisc;
    } else{
      m_xsecform = xsec_form::Gauss;
    }

    m_taumax = s["Tau_Max"].Get<double>();
    m_fixed_sigma = s["Fixed_Sigma"].Get<bool>();
    m_fixed_sigma_val = s["Fixed_Sigma_Value"].Get<double>();
    m_fixed_gamma = s["Fixed_Gamma"].Get<bool>();
    m_fixed_gamma_val = s["Fixed_Gamma_Value"].Get<double>();
    m_only_gluons = s["Only_Gluons"].Get<bool>();
    m_splitting_merging = s["Splitting_Merging"].Get<bool>();
    m_show_event_information = s["Show_Event_Information"].Get<bool>();
    m_show_progress_bar = s["Show_Progress_Bar"].Get<bool>();
    m_tau_restart = s["Tau_Restart"].Get<double>();
    m_tsample_min = s["tsample_Min"].Get<double>();
    m_tsample_max = s["tsample_Max"].Get<double>();
    m_include_bose_factors = s["Include_Bose_Factors"].Get<bool>();
    m_p_min = s["p_Min"].Get<double>();
    m_N_include = s["N_Include"].Get<int>();
    m_elastic_scattering = s["Elastic_Scattering"].Get<bool>();
    m_tau_max_scaling_limit = s["Tau_Max_Scaling_Limit"].Get<double>();
    m_alpha_s = s["Alpha_S"].Get<double>();
    m_gaussian_kT2 = s["Gaussian_kT2"].Get<bool>();
    m_kT2_reg = s["kT2_Reg"].Get<double>();
    m_f_r = s["f_r"].Get<double>();
    m_f_delta_p = s["f_Delta_p"].Get<double>();
    m_f_shell = s["f_Shell"].Get<bool>();
    m_f_N_max = s["f_N_Max"].Get<int>();
    m_formation_time = s["Formation_Time"].Get<double>();
    m_timekeeper = s["Timekeeper"].Get<int>();
    m_lambda = s["Lambda"].Get<double>();
    m_m2_min_scale = s["m2_Min_Scale"].Get<double>();

    m_OE_mult_scatter = s["OE_Mult_Scatter"].Get<double>();
    m_OE_mult_merge = s["OE_Mult_Merge"].Get<double>();
    m_OE_mult_split = s["OE_Mult_Split"].Get<double>();

    m_test_double = s["Test_Double"].Get<double>();
    m_test_bool = s["Test_Bool"].Get<bool>();
}
  
  double HeavyIon_Parameters::operator()(std::string keyword) {
  if (m_params.find(keyword)==m_params.end()) {
    msg_Error()<<"Error in "<<METHOD<<"("<<keyword<<") : "<<std::endl
	       <<"   Keyword not found, will exit the run."<<std::endl;
    abort();
  }
  return m_params[keyword];
}

std::ostream & ALPACA::operator<<(std::ostream & s, const xsec_form::code & xsecf) {
  if (xsecf==xsec_form::Gauss)          s<<"Gaussian";
  else if (xsecf==xsec_form::BlackDisc) s<<"black disc";
  else                                  s<<"unknown";
  return s;
}
