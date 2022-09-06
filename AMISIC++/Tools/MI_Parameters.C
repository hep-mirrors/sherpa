#include "AMISIC++/Tools/MI_Parameters.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Org/Exception.H"

using namespace AMISIC;
using namespace ATOOLS;
using namespace std;


MI_Parameters * AMISIC::mipars = NULL;

MI_Parameters::MI_Parameters() {}

bool MI_Parameters::Init()
{
  ReadParameters();
  return true;
}


void MI_Parameters::ReadParameters()
{
  auto s = Settings::GetMainSettings()["AMISIC"];
  m_parameters[string("pt_0(ref)")]
    = s["PT_0(ref)"].SetDefault(1.).Get<double>();
  m_parameters[string("pt_min(ref)")]
    = s["PT_Min(ref)"].SetDefault(1.5).Get<double>();
  m_parameters[string("eta")]
    = s["Eta"].SetDefault(0.16).Get<double>();
  m_parameters[string("Ecms(ref)")]
    = s["E(ref)"].SetDefault(7000.).Get<double>();
  double pt_0   = CalculatePT(m_parameters[string("pt_0(ref)")]);
  double pt_min = CalculatePT(m_parameters[string("pt_min(ref)")]);
  m_parameters[string("pt_min")]
    = s["PT_Min"].SetDefault(pt_min).Get<double>();
  m_parameters[string("pt_0")]
    = s["PT_0"].SetDefault(pt_0).Get<double>();
  m_scalescheme = s["MU_R_SCHEME"].SetDefault("PT").Get<scale_scheme::code>();
  m_parameters[string("RenScale_Factor")]
    = s["MU_R_FACTOR"].SetDefault(0.5).Get<double>();
  m_parameters[string("FacScale_Factor")]
    = s["MU_F_FACTOR"].SetDefault(1.0).Get<double>();
  m_parameters[string("SigmaND_Norm")]
    = s["SIGMA_ND_NORM"].SetDefault(0.4).Get<double>();
  m_parameters[string("Matter_Fraction1")]
    = s["MATTER_FRACTION1"].SetDefault(0.5).Get<double>();
  m_parameters[string("Matter_Radius1")]
    = s["MATTER_RADIUS1"].SetDefault(0.4).Get<double>();
  m_parameters[string("Matter_Radius2")]
    = s["MATTER_RADIUS2"].SetDefault(1.0).Get<double>();
  m_overlapform = s["MATTER_FORM"]
	  .SetDefault(overlap_form::Single_Gaussian)
	  .Get<overlap_form::code>();
  m_parameters[string("nPT_bins")]
    = s["nPT_bins"].SetDefault(200).Get<int>();
  m_parameters[string("nMC_points")]
    = s["nMC_points"].SetDefault(1000).Get<int>();
  m_parameters[string("nS_bins")]
    = s["nS_bins"].SetDefault(100).Get<int>();
}

double MI_Parameters::CalculatePT(const double & pt) {
  return pt * pow(rpa->gen.Ecms()/m_parameters[string("Ecms(ref)")],
		   m_parameters[string("eta")]);
}


double MI_Parameters::operator()(string keyword) 
{
  map<string,double>::iterator piter = m_parameters.find(keyword);
  if (piter!=m_parameters.end()) return piter->second;
  msg_Error()<<"Error in MI_Parameters("<<keyword<<") "
	     <<"in "<<m_parameters.size()<<".\n"
	     <<"   Keyword not found. Return 0 and hope for the best.\n";
  exit(1);
  return 0.;
}

std::ostream& AMISIC::operator<<(std::ostream& s, const overlap_form::code& f)
{
  switch (f) {
    case overlap_form::Single_Gaussian: return s << "Single_Gaussian";
    case overlap_form::Double_Gaussian: return s << "Double_Gaussian";
  }
  return s;
}

std::istream& AMISIC::operator>>(std::istream& s, overlap_form::code& f)
{
  std::string tag;
  s >> tag;
  if (tag == "Single_Gaussian")
    f = overlap_form::Single_Gaussian;
  else if (tag == "Double_Gaussian")
    f = overlap_form::Double_Gaussian;
  else
    THROW(fatal_error, "Unknown overlap form \"" + tag + "\"");
  return s;
}

std::ostream& AMISIC::operator<<(std::ostream& os, const scale_scheme::code& sc)
{
  switch (sc) {
    case scale_scheme::PT:           return os << "PT";
    case scale_scheme::PT_with_Raps: return os << "PT modified with rapidities";
  }
  return os;
}

std::istream& AMISIC::operator>>(std::istream& is, scale_scheme::code& sc)
{
  std::string tag;
  is >> tag;
  if (tag == "PT")
    sc = scale_scheme::PT;
  else if (tag == "PT_with_Raps")
    sc = scale_scheme::PT_with_Raps;
  else
    THROW(fatal_error, "Unknown scale scheme \"" + tag + "\"");
  return is;
}
