#include "AMISIC++/Tools/MI_Parameters.H"
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
  m_parameters[string("pt_0(ref)")]
    = s["PT_0(ref)"].SetDefault(2.05).Get<double>();
  m_parameters[string("pt_0(IR)")]
    = s["PT_0(IR)"].SetDefault(0.5).Get<double>();
  m_parameters[string("pt_min(ref)")]
    = s["PT_Min(ref)"].SetDefault(2.25).Get<double>();
  m_parameters[string("Ecms(ref)")]
    = s["E(ref)"].SetDefault(7000.).Get<double>();
  m_parameters[string("eta")]
    = s["Eta"].SetDefault(0.16).Get<double>();
  m_pt02ref   = sqr(m_parameters[string("pt_0(ref)")]);
  m_pt02IR    = sqr(m_parameters[string("pt_0(IR)")]);
  m_ptmin2ref = sqr(m_parameters[string("pt_min(ref)")]);
  m_Sref      = sqr(m_Eref = m_parameters[string("Ecms(ref)")]);
  m_Scms      = sqr(m_Ecms = rpa->gen.Ecms());
  m_eta       = m_parameters[string("eta")];
  double pt_0 = sqrt(CalculatePT02(m_Scms));
  m_parameters[string("pt_min")]
    = s["PT_Min"].SetDefault(m_parameters[string("pt_min(ref)")]).Get<double>();
  m_parameters[string("pt_0")]
    = s["PT_0"].SetDefault(pt_0).Get<double>();
  m_scalescheme = s["MU_R_SCHEME"].SetDefault("PT").Get<scale_scheme>();
  m_parameters[string("RenScale_Factor")]
    = s["MU_R_FACTOR"].SetDefault(0.5).Get<double>();
  m_parameters[string("FacScale_Factor")]
    = s["MU_F_FACTOR"].SetDefault(1.0).Get<double>();
  m_parameters[string("SigmaND_Norm")]
    = s["SIGMA_ND_NORM"].SetDefault(1.02).Get<double>();
  m_parameters[string("nPT_bins")]
    = s["nPT_bins"].SetDefault(200).Get<size_t>();
  m_parameters[string("nMC_points")]
    = s["nMC_points"].SetDefault(1000).Get<size_t>();
  m_parameters[string("nS_bins")]
    = s["nS_bins"].SetDefault(40).Get<size_t>();
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
    = s["Diffractive_s1"].SetDefault(20.).Get<double>();
  m_parameters[string("ElasticSlope_c0")]
    = s["ElasticSlope_c0"].SetDefault(2.24).Get<double>();
  m_parameters[string("ElasticSlope_c1")]
    = s["ElasticSlope_c1"].SetDefault(2.1).Get<double>();
}

double MI_Parameters::CalculatePT02(const double & s) const {
  return Max(m_pt02IR, m_pt02ref * pow((s<0 ? m_Scms : s)/m_Sref,m_eta));
}


double MI_Parameters::operator()(const string& keyword) const
{
  map<string,double>::const_iterator piter = m_parameters.find(keyword);
  if (piter!=m_parameters.end()) return piter->second;
  THROW(fatal_error,"Keyword not found in MI_Parameters.");
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

std::ostream& AMISIC::operator<<(std::ostream& s, const overlap_form& f)
{
  switch (f) {
    case overlap_form::Single_Gaussian: return s << "Single_Gaussian";
    case overlap_form::Double_Gaussian: return s << "Double_Gaussian";
    case overlap_form::unknown: return s << "Unknown";
  }
  return s;
}

std::istream& AMISIC::operator>>(std::istream& s, overlap_form& f)
{
  std::string tag;
  s >> tag;
  if (tag == "Single_Gaussian") f = overlap_form::Single_Gaussian;
  else if (tag == "Double_Gaussian")
    f = overlap_form::Double_Gaussian;
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
