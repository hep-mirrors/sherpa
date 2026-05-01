#include "AMISIC++/Tools/MI_Parameters.H"
#include "AMISIC++/Tools/Lookup_Tables.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Scoped_Settings.H"

using namespace AMISIC;
using namespace ATOOLS;
using namespace std;


const MI_Parameters * AMISIC::mipars = nullptr;

MI_Parameters::MI_Parameters() :
  m_pt02ref(0.), m_ptmin2ref(0.), m_Eref(0.), m_Sref(0.), m_Ecms(0.),
  m_Scms(0.), m_eta(0.)
{
  auto s = Settings::GetMainSettings()["AMISIC"];
  m_evttype =
    s["EVENT_TYPE"].SetDefault(evt_type::Perturbative).Get<evt_type::code>();
  m_parameters_vector[string("pt_0(ref)")]
    = s["PT_0(ref)"].SetDefault({2.05}).GetVector<double>();
  m_parameters[string("pt_0(ref)")]
    = m_parameters_vector[string("pt_0(ref)")][0];
  m_parameters[string("pt_0(IR)")]
    = s["PT_0(IR)"].SetDefault(0.5).Get<double>();
  m_parameters_vector[string("pt_min(ref)")]
    = s["PT_Min(ref)"].SetDefault({2.5}).GetVector<double>();
  m_parameters[string("pt_min(ref)")]
    = m_parameters_vector[string("pt_min(ref)")][0];
  m_parameters[string("pt_min(IR)")]
    = s["PT_Min(IR)"].SetDefault(1.00).Get<double>();
  m_parameters[string("Ecms(ref)")]
    = s["E(ref)"].SetDefault(7000.).Get<double>();
  m_parameters_vector[string("eta")]
    = s["Eta"].SetDefault({0.08}).GetVector<double>();
  m_parameters[string("eta")]
    = m_parameters_vector[string("eta")][0];
  m_pt02ref   = sqr(m_parameters[string("pt_0(ref)")]);
  m_pt02IR    = sqr(m_parameters[string("pt_0(IR)")]);
  for (size_t i=0; i<m_parameters_vector[string("pt_0(ref)")].size(); ++i)
  {
    m_pt02ref_variations.push_back(sqr(m_parameters_vector[string("pt_0(ref)")][i]));
  }
  for (size_t i=0; i<m_parameters_vector[string("pt_min(ref)")].size(); ++i)
  {
    m_ptmin2ref_variations.push_back(sqr(m_parameters_vector[string("pt_min(ref)")][i]));
  }
  m_pt02ref   = m_pt02ref_variations[0];
  m_ptmin2ref = m_ptmin2ref_variations[0];
  m_ptmin2IR  = sqr(m_parameters[string("pt_min(IR)")]);
  m_Sref      = sqr(m_Eref = m_parameters[string("Ecms(ref)")]);
  m_Scms      = sqr(m_Ecms = rpa->gen.Ecms());
  m_eta_variations = m_parameters_vector[string("eta")];
  m_eta       = m_parameters[string("eta")];
  std::vector<double> ptmin_variations;
  for (size_t ivar=0; ivar<std::max({m_parameters_vector[string("pt_min(ref)")].size(),
                               m_parameters_vector[string("eta")].size()}); ++ivar)
  {
    ptmin_variations.push_back(sqrt(CalculatePTmin2(m_Scms, ivar)));
  }
  double pt_min = ptmin_variations[0];
  std::vector<double> pt0_variations;
  for (size_t ivar=0; ivar<std::max({m_parameters_vector[string("pt_0(ref)")].size(),
                               m_parameters_vector[string("eta")].size()}); ++ivar)
  {
    pt0_variations.push_back(sqrt(CalculatePT02(m_Scms, ivar)));
  }
  double pt_0 = pt0_variations[0];
  m_parameters_vector[string("pt_min")]
    = s["PT_Min"].SetDefault(ptmin_variations).GetVector<double>();
  m_parameters[string("pt_min")]
    = m_parameters_vector[string("pt_min")][0];
  m_parameters_vector[string("pt_0")]
    = s["PT_0"].SetDefault(pt0_variations).GetVector<double>();
  m_parameters[string("pt_0")]
    = m_parameters_vector[string("pt_0")][0];
  m_scaleRscheme = s["MU_R_SCHEME"].SetDefault("PT").Get<scale_scheme>();
  m_scaleFscheme = s["MU_F_SCHEME"].SetDefault("PT").Get<scale_scheme>();
  m_parameters[string("RenScale_Factor")]
    = s["MU_R_FACTOR"].SetDefault(0.5).Get<double>();
  m_parameters[string("FacScale_Factor")]
    = s["MU_F_FACTOR"].SetDefault(1.0).Get<double>();
  m_parameters[string("E_min")] = s["E_Min"].SetDefault(0.25).Get<double>();
  m_parameters_vector[string("SigmaND_Norm")]
    = s["SIGMA_ND_NORM"].SetDefault({1.0}).GetVector<double>();
  m_parameters[string("SigmaND_Norm")]
    = m_parameters_vector[string("SigmaND_Norm")][0];
  m_parameters[string("PomeronIntercept")]
    = s["PomeronIntercept"].SetDefault(0.0808).Get<double>();
  m_parameters[string("PomeronSlope")]
    = s["PomeronSlope"].SetDefault(0.25).Get<double>();
  m_parameters[string("TriplePomeronCoupling")]
    = s["TriplePomeronCoupling"].SetDefault(0.318).Get<double>();
  m_parameters[string("ReggeonIntercept")]
    = s["ReggeonIntercept"].SetDefault(-0.4525).Get<double>();
  m_parameters[string("Diffractive_cres")]
    = s["Diffractive_cres"].SetDefault(2.).Get<double>();
  m_parameters[string("Diffractive_Mres")]
    = s["Diffractive_Mres"].SetDefault(2.).Get<double>();
  m_parameters[string("Diffractive_s1")]
    = s["Diffractive_s1"].SetDefault(400.).Get<double>();
  m_parameters[string("ElasticSlope_c0")]
    = s["ElasticSlope_c0"].SetDefault(2.24).Get<double>();
  m_parameters[string("ElasticSlope_c1")]
    = s["ElasticSlope_c1"].SetDefault(2.1).Get<double>();
  m_parameters[string("f_omega")]
    = s["f_omega"].SetDefault(0.15).Get<double>();
  m_parameters[string("phi_omega")]
    = s["phi_omega"].SetDefault(-0.25).Get<double>();
  m_parameters[string("f_nr")]
    = s["f_nr"].SetDefault(0.1675).Get<double>();
  m_parameters[string("Lambda_nr")]
    = s["Lambda_nr"].SetDefault(0.15).Get<double>();
  m_parameters[string("delta_nr")]
    = s["delta_nr"].SetDefault(0.75).Get<double>();
  m_parameters[string("max_reweight_factor")]
    = s["MAX_REWEIGHT_FACTOR"].SetDefault(-1.).Get<double>();

  m_flags[string("nPT_bins")]
    = s["nPT_bins"].SetDefault(200).Get<size_t>();
  m_flags[string("nMC_points")]
    = s["nMC_points"].SetDefault(2000).Get<size_t>();
  m_flags[string("nS_bins")]
    = s["nS_bins"].SetDefault(40).Get<size_t>();
  m_flags[string("nB_bins")]
    = s["nB_bins"].SetDefault(20).Get<size_t>();
  m_flags[string("nMaxScatters")]
    = s["N_MaxScatters"].SetDefault(10000).Get<size_t>();
  m_flags[string("weight_output")]
    = s["WEIGHT_OUTPUT"].SetDefault(0).Get<size_t>();
  m_flags[string("nB_samples")]
    = s["nB_samples"].SetDefault(0).Get<size_t>();


  size_t twopions = s["TwoPionInterference"].SetDefault(0).Get<size_t>();
  switch (twopions) {
  case 4:  m_twopions = two_pions::cont_only;      break;
  case 3:  m_twopions = two_pions::rho_omega_cont; break;
  case 2:  m_twopions = two_pions::rho_omega;      break;
  case 1:  m_twopions = two_pions::rho_only;       break;
  case 0:
  default: m_twopions = two_pions::none;           break;
  }
  m_triggerflavs = s["TRIGGER"].SetDefault(vector<int>()).GetVector<int>();
  if (m_triggerflavs.size()>2)
    THROW(critical_error,"Too many trigger flavours specified.");
}

double MI_Parameters::CalculatePT02(const double & s) const {
  return Max(m_pt02IR, m_pt02ref * pow((s<0 ? m_Scms : s)/m_Sref,2*m_eta));
}

double MI_Parameters::CalculatePT02(const double & s, size_t ivar) const {
  double pt02ref = m_pt02ref;
  double eta = m_eta;
  if (ivar < m_pt02ref_variations.size()) pt02ref = m_pt02ref_variations[ivar];
  if (ivar < m_eta_variations.size()) eta = m_eta_variations[ivar];
  return Max(m_pt02IR, pt02ref * pow((s<0 ? m_Scms : s)/m_Sref,2*eta));
}

double MI_Parameters::CalculatePTmin2(const double & s) const {
  return Max(m_ptmin2IR, m_ptmin2ref * pow((s<0 ? m_Scms : s)/m_Sref,2*m_eta));
}

double MI_Parameters::CalculatePTmin2(const double & s, size_t ivar) const {
  double ptmin2ref = m_ptmin2ref;
  double eta = m_eta;
  if (ivar < m_ptmin2ref_variations.size()) ptmin2ref = m_ptmin2ref_variations[ivar];
  if (ivar < m_eta_variations.size()) eta = m_eta_variations[ivar];
  return Max(m_ptmin2IR, ptmin2ref * pow((s<0 ? m_Scms : s)/m_Sref,2*eta));
}


double MI_Parameters::operator()(const string& keyword) const
{
  map<string,double>::const_iterator piter = m_parameters.find(keyword);
  if (piter!=m_parameters.end()) return piter->second;
  THROW(fatal_error,"Keyword not found in MI_Parameters.");
}

int MI_Parameters::operator[](const string& keyword) const
{
  map<string,int>::const_iterator piter = m_flags.find(keyword);
  if (piter!=m_flags.end()) return piter->second;
  THROW(fatal_error,"Keyword not found in MI_Parameters.");
}

std::vector<double> MI_Parameters::GetVariationVector(const std::string& keyword) const
{
  auto piter = m_parameters_vector.find(keyword);
  if (piter != m_parameters_vector.end()) return piter->second;
  THROW(fatal_error,"Keyword not found in MI_Parameters vector map.");
}

std::ostream& AMISIC::operator<<(std::ostream& s, const evt_type::code& f)
{
  switch (f) {
  case evt_type::code::Perturbative:  return s << "Perturbative";
  case evt_type::code::Elastic:       return s << "Elastic";
  case evt_type::code::DiffractiveA:  return s << "DiffractiveA";
  case evt_type::code::DiffractiveB:  return s << "DiffractiveB";
  case evt_type::code::DiffractiveAB: return s << "DiffractiveAB";
  case evt_type::code::QuasiElastic:  return s << "QuasiElastic";
  case evt_type::NonPerturbative:     return s << "Non-Perturbative";
  case evt_type::AllMinimumBias:      return s << "MinimumBias";
  }
  return s;
}

std::istream& AMISIC::operator>>(std::istream& s, evt_type::code& f)
{
  std::string tag;
  s >> tag;
  if (tag == "Perturbative")
    f = evt_type::code::Perturbative;
  else if (tag == "Elastic")
    f = evt_type::code::Elastic;
  else if (tag == "DiffractiveA")
    f = evt_type::code::DiffractiveA;
  else if (tag == "DiffractiveB")
    f = evt_type::code::DiffractiveB;
  else if (tag == "DiffractiveAB")
    f = evt_type::code::DiffractiveAB;
  else if (tag == "QuasiElastic")
    f = evt_type::code::QuasiElastic;
  else
    THROW(fatal_error, "Unknown overlap form \"" + tag + "\"");
  return s;
}

std::ostream& AMISIC::operator<<(std::ostream& os, const scale_scheme& sc)
{
  switch (sc) {
    case scale_scheme::PT: return os << "PT";
    case scale_scheme::PT_with_Raps: return os << "PT modified with rapidities";
  }
  return os;
}

std::istream& AMISIC::operator>>(std::istream& is, scale_scheme& sc)
{
  std::string tag;
  is >> tag;
  if (tag == "PT") sc = scale_scheme::PT;
  else if (tag == "PT_with_Raps")
    sc = scale_scheme::PT_with_Raps;
  else
    THROW(fatal_error, "Unknown scale scheme \"" + tag + "\"");
  return is;
}
