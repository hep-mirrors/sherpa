#include "Spectrum_EW.H"
#include "MathTools.H"
#include "Run_Parameter.H"

using namespace AMEGIC;
using namespace AORGTOOLS;
using namespace APHYTOOLS;

void Spectrum_EW::FillYukawas()
{
  if (!rpa.me.UsingModelMass()) return;

  Data_Read dr(rpa.GetPath()+std::string("/")+rpa.me.ModelFile());

  Flavour flav;
  flav = Flavour(kf::e); flav.SetYuk(dr.GetValue<double>("m_e-"));
  flav = Flavour(kf::mu); flav.SetYuk(dr.GetValue<double>("m_mu"));
  flav = Flavour(kf::tau); flav.SetYuk(dr.GetValue<double>("m_tau"));

  flav = Flavour(kf::d); flav.SetYuk(dr.GetValue<double>("m_down"));
  flav = Flavour(kf::u); flav.SetYuk(dr.GetValue<double>("m_up"));
  flav = Flavour(kf::s); flav.SetYuk(dr.GetValue<double>("m_strange"));
  flav = Flavour(kf::c); flav.SetYuk(dr.GetValue<double>("m_charm"));
  flav = Flavour(kf::b); flav.SetYuk(dr.GetValue<double>("m_bottom"));
  flav = Flavour(kf::t); flav.SetYuk(dr.GetValue<double>("m_top"));
  flav = Flavour(kf::h); flav.SetYuk(dr.GetValue<double>("m_H_SM"));

  double v         = dr.GetValue<double>("v");
  double Aqed      = dr.GetValue<double>("alpha_QED(MZ)");
  double sin2      = dr.GetValue<double>("sinTW^2");

  double mass_Z     = sqrt(M_PI*Aqed)*v/(sqrt(sin2)*sqrt(1.-sin2));
  double mass_W     = sqrt(M_PI*Aqed)*v/sqrt(sin2);

  flav = Flavour(kf::Z);flav.SetYuk(mass_Z);
  flav = Flavour(kf::W);flav.SetYuk(mass_W);
}











