#include "MCATNLO/Showers/SF_Lorentz.H"

#include "ATOOLS/Math/MathTools.H"
#include "MCATNLO/Showers/Splitting_Function_Base.H"
#include "MODEL/Main/Single_Vertex.H"
#include "ATOOLS/Org/Scoped_Settings.H"

using namespace MCATNLO;
using namespace ATOOLS;

double SF_Lorentz::s_kappa=2.0/3.0;

SF_Lorentz::SF_Lorentz(const SF_Key &key):
  p_ms(key.p_ms), p_cf(key.p_cf), m_col(0), m_pdfmin{ key.m_pdfmin }
{
  m_flavs[0]=key.p_v->in[0].Bar();
  if (key.m_mode==0) {
    m_flavs[1]=key.p_v->in[1];
    m_flavs[2]=key.p_v->in[2];
  }
  else {
    m_flavs[1]=key.p_v->in[2];
    m_flavs[2]=key.p_v->in[1];
  }
  Settings& s = Settings::GetMainSettings();
  std::string dipole_string = s["DIPOLES"]["CASE"].Get<std::string>();
  if      (dipole_string == "CS")    m_dipole_case = EXTAMP::DipoleCase::CS;
  else if (dipole_string == "IDa")   m_dipole_case = EXTAMP::DipoleCase::IDa;  // ee > bbWW with mapping a
  else if (dipole_string == "IDb")   m_dipole_case = EXTAMP::DipoleCase::IDb;  // ee > bbWW with mapping b
  else if (dipole_string == "IDin")  m_dipole_case = EXTAMP::DipoleCase::IDin; // pp > bbWW
  else if (dipole_string == "RES")   m_dipole_case = EXTAMP::DipoleCase::RES;  // ee > guu
  else if (dipole_string == "ID")    m_dipole_case = EXTAMP::DipoleCase::ID;   // ee > guu
  else                               m_dipole_case = EXTAMP::DipoleCase::CS;
}

SF_Lorentz::~SF_Lorentz() {}

double SF_Lorentz::Lambda
(const double &a,const double &b,const double &c) const
{
  return a*a+b*b+c*c-2.*(a*b+a*c+b*c);
}

double SF_Lorentz::AsymmetryFactor
(const double z,const double y,const double Q2)
{
  return 1.0;
}

double SF_Lorentz::Scale(const double z,const double y,
			 const double scale,const double Q2) const
{
  return scale;
}

double SF_Lorentz::JFF(const double &y,const double &mui2,
		       const double &muj2,const double &muk2,
		       const double &muij2)
{
  m_lastJ = (1.-y)*sqr(1.0-mui2-muj2-muk2)/sqrt(Lambda(1.0,muij2,muk2));
  return m_lastJ;
}

double SF_Lorentz::JFI(const double &y,const double &eta,
		       const double &scale,const Cluster_Amplitude *const sub)
{
  if (sub) {
    m_lastJ = 1.0;
  } else {
    const double fresh = p_sf->GetXPDF(scale,eta/(1.0-y),m_flspec,m_beam);
    const double old = p_sf->GetXPDF(scale,eta,m_flspec,m_beam);
    if (fresh < 0.0 || old < 0.0 || !PDFValueAllowedAsDenominator(old, eta))
      m_lastJ = 0.;
    else
      m_lastJ = (1.0 - y) * fresh / old;
  }
  return m_lastJ;
}

double SF_Lorentz::JIF(const double &z,const double &y,const double &eta,
		       const double &scale,const Cluster_Amplitude *const sub)
{
  if (sub) {
    m_lastJ = 1.0 / z;
  } else {
    const double fresh = p_sf->GetXPDF(scale,eta/z,m_flavs[0],m_beam);
    const double old = p_sf->GetXPDF(scale,eta,m_flavs[1],m_beam);
    if (fresh < 0.0 || old < 0.0 || !PDFValueAllowedAsDenominator(old, eta))
      m_lastJ = 0.;
    else
      m_lastJ = fresh / old;
  }
  return m_lastJ;
}

double SF_Lorentz::JII(const double &z,const double &y,const double &eta,
		       const double &scale,const Cluster_Amplitude *const sub)
{
  if (sub) {
    m_lastJ = 1.0/z;
  } else {
    const double fresh = p_sf->GetXPDF(scale,eta/z,m_flavs[0],m_beam);
    const double old = p_sf->GetXPDF(scale,eta,m_flavs[1],m_beam);
    if (fresh < 0.0 || old < 0.0 || !PDFValueAllowedAsDenominator(old, eta))
      m_lastJ = 0.;
    else
      m_lastJ = fresh / old;
  }
  return m_lastJ;
}

bool SF_Lorentz::PDFValueAllowedAsDenominator(const double& val,
                                              const double& eta)
{
  const double dynamic_pdf_threshold{
    m_pdfmin.first * log(1.0 - eta) / log(1.0 - m_pdfmin.second) };
  return (std::abs(val) > dynamic_pdf_threshold);
}
