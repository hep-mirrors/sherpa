#include "BEAM/Spectra/EPA.H"
#include "BEAM/Spectra/EPA_FF.H"
#include "BEAM/Spectra/EPA_Spectra_Plotter.H"

#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Settings.H"

#include <fstream>
#include <string>

using namespace BEAM;
using namespace ATOOLS;
using namespace std;

EPA::EPA(const Flavour& beam, const double energy,
	 const double pol, const int dir) :
  Beam_Base(beamspectrum::EPA, beam, energy, pol, dir),
  m_type(EPA_ff_type::point), p_ff(nullptr),
  m_mass(m_beam.Mass(true)), m_charge(m_beam.Charge()),
  m_gamma(m_energy/m_mass), m_plotting(0)
{
  Initialise();
  m_Nbunches   = 2;
  m_bunches.resize(m_Nbunches);
  m_bunches[0] = Flavour(kf_photon);
  m_bunches[1] = m_beam;
  m_vecouts.resize(m_Nbunches);
  m_vecouts[0] = Vec4D(m_energy, 0., 0., m_dir * m_energy);
  m_vecouts[1] = Vec4D(0.,0.,0.,0.);
  m_on         = true;
}

bool EPA::CalculateWeight(double x, double q2) {
  m_x = x; m_q2 = q2;
  m_weight = (x >= m_xmin && x <= m_xmax)
                     ? ATOOLS::Max(0., m_pref / x * p_ff->N(x))
                     : 0.;
  if (IsNan(m_weight)) msg_Out()<<"Boink! "<<METHOD<<"(x = "<<x<<") yields NaN.\n";
  return true;
}

void EPA::FixPosition() {
  double R = p_ff->SelectB(m_x), phi = 2.*M_PI*ran->Get();
  m_position = R * Vec4D(0., cos(phi), sin(phi), 0.);
}

void EPA::SetOutMomentum(const ATOOLS::Vec4D &out, const size_t & i) {
  if (i==0) {
    m_vecouts[0] = out;
    m_vecouts[1] = m_lab-out;
    m_x          = out[0]/m_lab[0];
    m_q2         = dabs(out.Abs2());
  }
}

void EPA::Initialise() {
  Settings &s = Settings::GetMainSettings();
  RegisterDefaults();
  m_aqed      = s["EPA"]["AlphaQED"].Get<double>();
  m_pref      = sqr(m_charge)*m_aqed/M_PI;
  m_approx    = s["EPA"]["Approximation"].Get<size_t>();
  m_analytic  = s["EPA"]["AnalyticFF"].Get<bool>();
  m_plotting  = s["EPA"]["PlotSpectra"].Get<bool>();
  m_q2max     = ExtractParameter(s,"Q2Max");
  m_q2min     = ExtractParameter(s,"Q2Min");
  m_theta_max = ExtractParameter(s,"ThetaMax");
  m_pt2max    = sqr(m_energy*m_theta_max);
  m_pt2min    = ExtractParameter(s,"PT2Min");
  m_xmin      = ExtractParameter(s,"xMin");
  m_xmax      = ExtractParameter(s,"xMax");
  m_bmin      = ExtractParameter(s,"bMin");
  m_bmax      = ExtractParameter(s,"bMax");
  m_WSd       = ExtractParameter(s,"WoodSaxon_d");
  m_mu        = ExtractParameter(s,"MagneticMu");
  m_Lambda2   = ExtractParameter(s,"Lambda2");
  m_Q02       = ExtractParameter(s,"Q02");
  m_nxbins    = s["EPA"]["xBins"].Get<int>();
  m_nbbins    = s["EPA"]["bBins"].Get<int>();

  InitFormFactor(s);
  //WriteDebugFiles(s);
  if (m_plotting>0) {
    EPA_Spectra_Plotter plotter(this,string("Spectra"));
    plotter(m_plotting);
    Tests();
  }
}

void EPA::RegisterDefaults() const {
  const auto &s = Settings::GetMainSettings()["EPA"];
  s["Q2Max"].SetDefault(3.0);
  s["Q2Min"].SetDefault(-1.);
  s["xMax"].SetDefault(1.);
  s["xMin"].SetDefault(1.e-6);
  s["xBins"].SetDefault(12);
  s["bMin"].SetDefault(0.);
  s["bMax"].SetDefault(1.e12);
  s["bBins"].SetDefault(10);
  s["PT2Min"].SetDefault(0.0);
  s["Form_Factor"].SetDefault(
   size_t(m_beam.IsIon()       ? EPA_ff_type::WoodSaxon
          : m_beam.IsNucleon() ? EPA_ff_type::dipole
          : m_beam.IsMeson()   ? EPA_ff_type::dipole
                               : EPA_ff_type::point));
  s["MagneticMu"].SetDefault(m_beam.IsNucleon() ? 2.79 : 0.);
  s["Lambda2"].SetDefault(0.71);
  s["Q02"].SetDefault(2.*rpa->hBar_c()/m_beam.Radius());
  s["WoodSaxon_d"].SetDefault(0.5);
  s["AlphaQED"].SetDefault(0.0072992701);
  s["ThetaMax"].SetDefault(0.3);
  s["Approximation"].SetDefault(1);
  s["AnalyticFF"].SetDefault(true);
  s["PlotSpectra"].SetDefault(0);
  s["Debug"].SetDefault(false);
  s["Debug_Files"].SetDefault("EPA_debugOutput");
}

double EPA::ExtractParameter(Settings &s,const std::string & tag) {
  std::vector<double> parms = s["EPA"][tag].GetVector<double>();
  if (parms.size()!=1 && parms.size()!=2)
    THROW(fatal_error,
	  "Specify either one or two values for 'EPA:"+tag+"'.  Will exit.");
  double parm = (m_dir > 0) ? parms.front() : parms.back();
  if (tag=="PTMin" && parm>1.0) {
    /* pt2min > 1 - according to approximation of
       'qmi' calculation in CalculateWeight */
    THROW(critical_error, "Too big p_T cut 'EPA:"+tag+"'.  Will exit.");
  }
  return parm;
}

void EPA::InitFormFactor(Settings &s) {
  std::vector<int> formfactors = s["EPA"]["Form_Factor"].GetVector<int>();
  if (formfactors.size()!=1 && formfactors.size()!=2)
    THROW(fatal_error,
          "Specify either one or two values for `EPA:Form_Factor'.");
  int formfactor = (m_dir > 0) ? formfactors.front() : formfactors.back();
  switch (formfactor) {
  case  0:
    m_type = EPA_ff_type::point;
    p_ff   = new EPA_Point(m_beam);
    break;
  case  1:
    m_type = EPA_ff_type::Gauss;
    p_ff   = new EPA_Gauss(m_beam);
    p_ff->SetParam("Q02",m_Q02);
    p_ff->SetParam("Mu2",sqr(m_mu));
    break;
  case  2:
    m_type = EPA_ff_type::dipole;
    p_ff   = new EPA_Dipole(m_beam);
    p_ff->SetParam("Lambda2",m_Lambda2);
    p_ff->SetParam("Mu2",sqr(m_mu));
    break;
  case 13:
    m_type = EPA_ff_type::WoodSaxon;
    p_ff   = new EPA_WoodSaxon(m_beam);
    p_ff->SetParam("WSd",m_WSd);
    break;
  default:
    THROW(fatal_error,
          "unspecified EPA form factor: "+ToString(formfactor));
  }
  p_ff->SetSwitch("approximation",m_approx);
  p_ff->SetSwitch("analytic",m_analytic);
  p_ff->SetQ2Range(m_q2min,m_q2max);
  p_ff->SetPT2Range(m_pt2min,m_pt2max);
  p_ff->FillTables(m_nxbins,m_nbbins);
}
