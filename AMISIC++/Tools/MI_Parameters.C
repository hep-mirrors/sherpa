#include "AMISIC++/Tools/MI_Parameters.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"

using namespace AMISIC;
using namespace ATOOLS;
using namespace std;


MI_Parameters * AMISIC::mipars = NULL;

MI_Parameters::MI_Parameters() {}

bool MI_Parameters::Init(Default_Reader *const defaultreader)
{
  ReadParameters(defaultreader);
  return true;
}


void MI_Parameters::ReadParameters(Default_Reader *const defaultreader)
{
  m_parameters[string("pt_0(ref)")] =
    defaultreader->GetValue<double>("AMISIC::PT_0(ref)",1.);
  m_parameters[string("pt_min(ref)")]    =
    defaultreader->GetValue<double>("AMISIC::PT_Min(ref)",1.5);
  m_parameters[string("eta")]       =
    defaultreader->GetValue<double>("AMISIC::Eta",0.16);
  m_parameters[string("Ecms(ref)")] =
    defaultreader->GetValue<double>("AMISIC::E(ref)",7000.);
  double pt_0   = CalculatePT(m_parameters[string("pt_0(ref)")]);
  double pt_min = CalculatePT(m_parameters[string("pt_min(ref)")]);
  m_parameters[string("pt_min")]    =
    defaultreader->GetValue<double>("AMISIC::PT_Min",pt_min);
  m_parameters[string("pt_0")]      =
    defaultreader->GetValue<double>("AMISIC::PT_0",pt_0);
  string scheme = defaultreader->GetValue<string>("AMISIC::MU_R_SCHEME","PT");
  m_scalescheme = ToType<scale_scheme::code>(scheme);
  m_parameters[string("RenScale_Factor")] = 
    defaultreader->GetValue<double>("AMISIC::MU_R_FACTOR",0.5);
  m_parameters[string("FacScale_Factor")] = 
    defaultreader->GetValue<double>("AMISIC::MU_F_FACTOR",1.0);
  m_parameters[string("Matter_Fraction1")] =
    defaultreader->GetValue<double>("AMISIC::MATTER_FRACTION1",0.5);
  m_parameters[string("Matter_Radius1")] =
    defaultreader->GetValue<double>("AMISIC::MATTER_RADIUS1",0.4);
  m_parameters[string("Matter_Radius2")] =
    defaultreader->GetValue<double>("AMISIC::MATTER_RADIUS2",1.0);
  string form = defaultreader->GetValue<string>("AMISIC::MATTER_FORM",
						"Double_Gaussian");
  m_overlapform = ToType<overlap_form::code>(form);
  m_parameters[string("nPT_bins")]    =
    defaultreader->GetValue<int>("AMISIC::nPT_bins",200);
  m_parameters[string("nMC_points")]    =
    defaultreader->GetValue<int>("AMISIC::nMC_points",1000);
}

double MI_Parameters::CalculatePT(const double & pt) {
  return pt * (pow(rpa->gen.Ecms()/m_parameters[string("Ecms(ref)")],
		   m_parameters[string("eta")]));
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
