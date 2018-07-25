#include "REMNANTS/Tools/Remnant_Parameters.H"
#include "ATOOLS/Org/Default_Reader.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"

using namespace REMNANTS;
using namespace ATOOLS;
using namespace std;


Remnant_Parameters * REMNANTS::repars = NULL;

Remnant_Parameters::Remnant_Parameters() :
  m_defform("gauss_limited"),
  m_defmean(0.0), m_defsigma(1.0), m_refE(7000.), m_scaleexpo(0.55),
  m_defQ2(0.77), m_defktmax(5.), m_defeta(5.)
{}

bool Remnant_Parameters::Init(Default_Reader *const defaultreader)
{
  msg_Out()<<"In Remnant_Parameters::Init"<<endl;
  ReadParameters(defaultreader);
  return true;
}

void Remnant_Parameters::ReadParameters(Default_Reader *const defaultreader)
{
  // Setting the mean and width of the Gaussian - we could discuss
  // more functional forms of the distribution here.  Default for leptons
  // is zero transverse momentum, only hadrons have a kT distribution of
  // remnants - again, this may have to change or become more
  // sophisticated for photons with hadronic structure.
  for (size_t beam=0;beam<2;beam++) {
    string beamtag   = to_string(beam);
    string formtag   = "K_PERP_FORM_"+to_string(beam);
    string meantag   = "K_PERP_MEAN_"+to_string(beam);
    string sigmatag  = "K_PERP_SIGMA_"+to_string(beam);
    string Q2tag     = "K_PERP_Q2_"+to_string(beam);
    string refEtag   = "K_PERP_REFE_"+to_string(beam);
    string expotag   = "K_PERP_SCALE_EXPO_"+to_string(beam);
    string ktmaxtag  = "K_PERP_MAX_"+to_string(beam);
    string ktexpotag = "K_PERP_CUT_EXPO_"+to_string(beam);
    double refE      = defaultreader->GetValue<double>(refEtag,m_refE);
    double expo      = defaultreader->GetValue<double>(expotag,m_scaleexpo);
    double mean      = defaultreader->GetValue<double>(meantag,m_defmean);
    double sigma     = defaultreader->GetValue<double>(sigmatag,m_defsigma);
    double Q2        = defaultreader->GetValue<double>(Q2tag,m_defQ2);
    m_form[beam]     = SelectForm(defaultreader->GetValue<string>(formtag,m_defform));
    m_parameters["kt_mean_"+beamtag]  = mean;
    m_parameters["kt_sigma_"+beamtag] = sigma * pow((rpa->gen.Ecms()/refE),expo);
    m_parameters["Q2_"+beamtag]       = Q2 * pow((rpa->gen.Ecms()/refE),expo);
    m_parameters["kt_max_"+beamtag]   = defaultreader->GetValue<double>(ktmaxtag,m_defktmax);
    m_parameters["eta_"+beamtag]      = defaultreader->GetValue<double>(ktexpotag,m_defeta);
  }
  // this is for the distribution of the beam remnant partons/spectators in the
  // transverse plane - relevant for colour reconnections.
  m_parameters[string("Matter_Fraction1")] =
    defaultreader->GetValue<double>("AMISIC::MATTER_FRACTION1",0.8);
  m_parameters[string("Matter_Radius1")] =
    defaultreader->GetValue<double>("AMISIC::MATTER_RADIUS1",1.0);
  m_parameters[string("Matter_Radius2")] =
    defaultreader->GetValue<double>("AMISIC::MATTER_RADIUS2",1.0);
  string form = defaultreader->GetValue<string>("AMISIC::MATTER_FORM",
						"Single_Gaussian");
  if (form=="Double_Gaussian") m_overlapform = overlap_form::Double_Gaussian;
  if (form=="Single_Gaussian") m_overlapform = overlap_form::Single_Gaussian;
}

prim_kperp_form::code Remnant_Parameters::SelectForm(const std::string & form) {
  prim_kperp_form::code pkf = prim_kperp_form::undefined;
  if      (form=="gauss")          pkf = prim_kperp_form::gauss;
  else if (form=="gauss_limited")  pkf = prim_kperp_form::gauss_limited;
  else if (form=="dipole")         pkf = prim_kperp_form::dipole;
  else if (form=="dipole_limited") pkf = prim_kperp_form::dipole_limited;
  return pkf;
}

double Remnant_Parameters::operator()(string keyword) 
{
  map<string,double>::iterator piter = m_parameters.find(keyword);
  if (piter!=m_parameters.end()) return piter->second;
  msg_Error()<<"Error in Remnant_Parameters("<<keyword<<") "
	     <<"in "<<m_parameters.size()<<".\n"
	     <<"   Keyword not found. Return 0 and hope for the best.\n";
  exit(1);
  return 0.;
}

