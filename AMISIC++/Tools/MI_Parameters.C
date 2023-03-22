#include "AMISIC++/Tools/MI_Parameters.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Scoped_Settings.H"

using namespace AMISIC;
using namespace ATOOLS;
using namespace std;


const MI_Parameters * AMISIC::mipars = NULL;

MI_Parameters::MI_Parameters() :
  m_pt02ref(0.), m_ptmin2ref(0.),
  m_Sref(0.), m_Eref(0.), m_Scms(0.), m_Ecms(0.), m_eta(0.)
{
  auto s = Settings::GetMainSettings()["AMISIC"];
  m_parameters[string("pt_0(ref)")]
    = s["PT_0(ref)"].SetDefault(2.2).Get<double>();
  m_parameters[string("pt_0(IR)")]
    = s["PT_0(IR)"].SetDefault(0.5).Get<double>();
  m_parameters[string("pt_min(ref)")]
    = s["PT_Min(ref)"].SetDefault(2.5).Get<double>();
  m_parameters[string("Ecms(ref)")]
    = s["E(ref)"].SetDefault(7000.).Get<double>();
  m_parameters[string("eta")]
    = s["Eta"].SetDefault(0.08).Get<double>();
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
  m_scalescheme = s["MU_R_SCHEME"].SetDefault("PT").Get<scale_scheme::code>();
  m_parameters[string("RenScale_Factor")]
    = s["MU_R_FACTOR"].SetDefault(0.5).Get<double>();
  m_parameters[string("FacScale_Factor")]
    = s["MU_F_FACTOR"].SetDefault(1.0).Get<double>();
  m_parameters[string("SigmaND_Norm")]
    = s["SIGMA_ND_NORM"].SetDefault(1.0).Get<double>();
  m_parameters[string("Matter_Fraction1")]
    = s["MATTER_FRACTION1"].SetDefault(0.5).Get<double>();
  m_parameters[string("Matter_Radius1")]
    = s["MATTER_RADIUS1"].SetDefault(1.0).Get<double>();
  m_parameters[string("Matter_Radius2")]
    = s["MATTER_RADIUS2"].SetDefault(2.0).Get<double>();
  m_overlapform = s["MATTER_FORM"]
	  .SetDefault(overlap_form::code::Single_Gaussian)
	  .Get<overlap_form::code>();
  m_parameters[string("nPT_bins")]
    = s["nPT_bins"].SetDefault(200).Get<size_t>();
  m_parameters[string("nMC_points")]
    = s["nMC_points"].SetDefault(1000).Get<size_t>();
  m_parameters[string("nS_bins")]
    = s["nS_bins"].SetDefault(100).Get<size_t>();
  m_parameters[string("PomeronIntercept")]
    = s["PomeronIntercept"].SetDefault(0.0808).Get<double>();
  m_parameters[string("PomeronSlope")]
    = s["PomeronSlope"].SetDefault(0.25).Get<double>();
  m_parameters[string("TriplePomeronCoupling")]
    = s["TriplePomeronCoupling"].SetDefault(0.318).Get<double>();
  m_parameters[string("ReggeonIntercept")]
    = s["ReggeonIntercept"].SetDefault(-0.4525).Get<double>();
}

double MI_Parameters::CalculatePT02(const double & s) const {
  return Max(m_pt02IR, m_pt02ref * pow((s<0 ? m_Scms : s)/m_Sref,m_eta));
}

double MI_Parameters::operator()(const string& keyword) const
{
  map<string,double>::const_iterator piter = m_parameters.find(keyword);
  if (piter!=m_parameters.end()) return piter->second;
  THROW(fatal_error,"Keyword "+keyword+" not found.");
  return 0.;
}

std::ostream& AMISIC::operator<<(std::ostream& s, const overlap_form::code& f)
{
  switch (f) {
    case overlap_form::code::Single_Gaussian: return s << "Single_Gaussian";
    case overlap_form::code::Double_Gaussian: return s << "Double_Gaussian";
  }
  return s;
}

std::istream& AMISIC::operator>>(std::istream& s, overlap_form::code& f)
{
  std::string tag;
  s >> tag;
  if (tag == "Single_Gaussian")
    f = overlap_form::code::Single_Gaussian;
  else if (tag == "Double_Gaussian")
    f = overlap_form::code::Double_Gaussian;
  else
    THROW(fatal_error, "Unknown overlap form \"" + tag + "\"");
  return s;
}

std::ostream& AMISIC::operator<<(std::ostream& os, const scale_scheme::code& sc)
{
  switch (sc) {
    case scale_scheme::code::PT:           return os << "PT";
    case scale_scheme::code::PT_with_Raps: return os << "PT modified with rapidities";
  }
  return os;
}

std::istream& AMISIC::operator>>(std::istream& is, scale_scheme::code& sc)
{
  std::string tag;
  is >> tag;
  if (tag == "PT")
    sc = scale_scheme::code::PT;
  else if (tag == "PT_with_Raps")
    sc = scale_scheme::code::PT_with_Raps;
  else
    THROW(fatal_error, "Unknown scale scheme \"" + tag + "\"");
  return is;
}
