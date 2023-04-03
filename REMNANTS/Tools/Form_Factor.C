#include "REMNANTS/Tools/Form_Factor.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Org/Exception.H"

using namespace REMNANTS;
using namespace ATOOLS;

std::ostream& REMNANTS::operator<<(std::ostream& s, const REMNANTS::matter_form::code& f)
{
  switch (f) {
    case matter_form::code::Single_Gaussian: return s << "Single_Gaussian";
    case matter_form::code::Double_Gaussian: return s << "Double_Gaussian";
    case matter_form::code::unknown: return s << "Unknown";
    }
  return s;
}

std::istream& REMNANTS::operator>>(std::istream& s, REMNANTS::matter_form::code& f)
{
  std::string tag;
  s >> tag;
  if (tag == "Single_Gaussian")
    f = matter_form::code::Single_Gaussian;
  else if (tag == "Double_Gaussian")
    f = matter_form::code::Double_Gaussian;
  else
    THROW(fatal_error, "Unknown matter form \"" + tag + "\"");
  return s;
}


Form_Factor::Form_Factor(const Flavour & flav) :
  m_flav(flav), m_form(matter_form::code::Single_Gaussian),
  m_fraction1(1.), m_radius1(1.), m_radius2(0.)
{
  SetDefaults(m_flav);
  Initialise();
}

void Form_Factor::Initialise() {
  Scoped_Settings data, alldata = Settings::GetMainSettings()["REMNANTS"];
  std::vector<std::string> pids = alldata.GetKeys(), keys;
  for (std::vector<std::string>::iterator pid=pids.begin();pid!=pids.end();pid++) {
    kf_code kf = ToType<kf_code>((*pid));
    if (kf!=m_flav.Kfcode()) continue;
    keys       = alldata[(*pid)].GetKeys();
    /////////////////////////////////////////////////////////////////////////////
    // Fix the overall form of the matter distribution first - this is a bit
    // unelegant here.
    /////////////////////////////////////////////////////////////////////////////
    bool found = false;
    for (std::vector<std::string>::iterator key=keys.begin();key!=keys.end();key++) {
      if ((*key)=="FORM") found = true;
    }
    if (found) {
      data = alldata[(*pid)]["FORM"];
      std::string form = data.SetDefault("Single_Gaussian").Get<std::string>();
      if (form=="Double_Gaussian")      m_form = matter_form::code::Double_Gaussian;
      else if (form=="Single_Gaussian") m_form = matter_form::code::Single_Gaussian;
      else {
	msg_Error()<<"Error in "<<METHOD<<"(pid = "<<(*pid)<<"): "
		   <<"matter form "<<form<<" not recognised, will assume Single_Gaussian.\n";
	m_form = matter_form::code::Single_Gaussian;
      }
    }
    /////////////////////////////////////////////////////////////////////////////
    // Now set the parameters.
    // Make sure that they are consistent with the overall form
    /////////////////////////////////////////////////////////////////////////////
    for (std::vector<std::string>::iterator key=keys.begin();key!=keys.end();key++) {
      data = alldata[(*pid)][(*key)];
      if ((*key)=="FORM") continue;
      else if ((*key)=="RADIUS1")   m_radius1   = data.SetDefault(m_radius1).Get<double>();
      else if ((*key)=="RADIUS2")   m_radius2   = ((m_form==matter_form::code::Double_Gaussian) ?
						   data.SetDefault(m_radius2).Get<double>() : 0.);
      else if ((*key)=="FRACTION1") m_fraction1 = ((m_form==matter_form::code::Double_Gaussian) ?
						   data.SetDefault(m_fraction1).Get<double>() : 1.);
    }
  }
}

void Form_Factor::SetDefaults(const Flavour & flav) {
  if (flav.IsNucleon()) {
    /////////////////////////////////////////////////////////////////////////////
    // we use the mean charge radius of the proton here - may need to revise
    /////////////////////////////////////////////////////////////////////////////
    m_form      = matter_form::code::Single_Gaussian;
    m_fraction1 = 1.;
    m_radius1   = 0.86;
    m_radius2   = 0.00;
  }
  else if (flav.IsPhoton()) {
    /////////////////////////////////////////////////////////////////////////////
    // we use the mean charge radius of the rho here - may need to revise
    /////////////////////////////////////////////////////////////////////////////
    m_form      = matter_form::code::Single_Gaussian;
    m_fraction1 = 1.;
    m_radius1   = 0.75;
    m_radius2   = 0.00;
  }
}


Vec4D Form_Factor::operator()() {
  ///////////////////////////////////////////////////////////////////////////////
  // Generate a position distributed according to the form-factor 
  ///////////////////////////////////////////////////////////////////////////////
  double radius = (m_form==matter_form::code::Double_Gaussian &&
		   ran->Get()<=m_fraction1) ? m_radius1 : m_radius2;
  double x1 = ran->GetGaussian(), x2 = ran->GetGaussian();
  return Vec4D(0.,radius*x1,radius*x2,0.);
}
