#include "PDF/Main/PDF_Base.H"

#define COMPILE__Getter_Function
#define OBJECT_TYPE PDF::PDF_Base
#define PARAMETER_TYPE PDF::PDF_Arguments
#include "ATOOLS/Org/Getter_Function.C"

#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Scoped_Settings.H"

using namespace PDF;
using namespace ATOOLS;

namespace PDF {
  PDF_Defaults *pdfdefs(nullptr);
}

PDF_Defaults::PDF_Defaults()
{
  #ifdef USING__LHAPDF
    m_deflib[kf_p_plus]  = "LHAPDFSherpa";
    m_defset[kf_p_plus]  = "PDF4LHC21_40_pdfas";
    m_allowed_flavours.insert(kf_p_plus);
    m_deflib[kf_n]       = "LHAPDFSherpa";
    m_defset[kf_n]       = "PDF4LHC21_40_pdfas";
    m_allowed_flavours.insert(kf_n);
    m_deflib[kf_Lambda]  = "LHAPDFSherpa";
    m_defset[kf_Lambda]  = "PDF4LHC21_40_pdfas";
    m_allowed_flavours.insert(kf_Lambda);
    m_deflib[kf_pi_plus] = "LHAPDFSherpa";
    m_defset[kf_pi_plus] = "JAM21PionPDFnlo";
    m_allowed_flavours.insert(kf_pi_plus);
    m_deflib[kf_K_plus]  = "LHAPDFSherpa";
    m_defset[kf_K_plus]  = "JAM21PionPDFnlo";
    m_allowed_flavours.insert(kf_K_plus);
    m_deflib[kf_reggeon] = "LHAPDFSherpa";
    m_defset[kf_reggeon] = "GRVPI0";
    m_allowed_flavours.insert(kf_reggeon);
  #else
    m_deflib[kf_p_plus]  = "NNPDFSherpa";
    m_defset[kf_p_plus]  = "NNPDF31_nnlo_as_0118_mc";
    m_allowed_flavours.insert(kf_p_plus);
    m_deflib[kf_n]       = "NNPDFSherpa";
    m_defset[kf_n]       = "NNPDF31_nnlo_as_0118_mc";
    m_allowed_flavours.insert(kf_n);
  #endif
  m_deflib[kf_e]       = "PDFESherpa";
  m_defset[kf_e]       = "PDFe";
  m_allowed_flavours.insert(kf_e);
  m_deflib[kf_mu]      = "PDFESherpa";
  m_defset[kf_mu]      = "PDFe";
  m_allowed_flavours.insert(kf_mu);
  m_deflib[kf_photon]  = "SASGSherpa";
  m_defset[kf_photon]  = "SAS1M";
  m_allowed_flavours.insert(kf_photon);
  m_deflib[kf_pomeron] = "H1Sherpa";
  m_defset[kf_pomeron] = "FitB";
  m_allowed_flavours.insert(kf_pomeron);
}

std::ostream &PDF::operator<<(std::ostream &ostr,const PDF::PDF_AS_Info &asi)
{
  return ostr<<"\\alpha_s of order="<<asi.m_order
             <<" with \\alpha_s(\\mu="<<sqrt(asi.m_mz2)<<")="<<asi.m_asmz;
}

void PDF_Id_Maps::InitMaps() {
  // Nucleons (octets only) derived from the proton PDF:
  // - neutrons are isospin-invariant (up to QED effects)
  // - for all others we assume SU_F(3) symmetry.  This is a manifestly
  //   bad approximation but better than nothing
  m_maps[Flavour(kf_p_plus)]         = PDF_Id_List();
  m_maps[Flavour(kf_n)]              = PDF_Id_List();
  m_maps[Flavour(kf_n)][1]           = 2;  // d -> u
  m_maps[Flavour(kf_n)][2]           = 1;  // u -> d
  m_maps[Flavour(kf_Sigma_minus)]    = PDF_Id_List();
  m_maps[Flavour(kf_Sigma_minus)][1] = 2;  // d -> u
  m_maps[Flavour(kf_Sigma_minus)][2] = 3;  // u -> s
  m_maps[Flavour(kf_Sigma_minus)][3] = 1;  // s -> d
  m_maps[Flavour(kf_Sigma)]          = PDF_Id_List();
  m_maps[Flavour(kf_Sigma)][2]       = 1;  // u -> d
  m_maps[Flavour(kf_Sigma)][3]       = 1;  // s -> d
  m_maps[Flavour(kf_Lambda)]         = PDF_Id_List();
  m_maps[Flavour(kf_Lambda)][2]      = 1;  // u -> d
  m_maps[Flavour(kf_Lambda)][3]      = 1;  // s -> d
  m_maps[Flavour(kf_Sigma_plus)]     = PDF_Id_List();
  m_maps[Flavour(kf_Sigma_plus)][1]  = 3;  // d -> s
  m_maps[Flavour(kf_Sigma_plus)][3]  = 1;  // s -> d
  m_maps[Flavour(kf_Xi)]             = PDF_Id_List();
  m_maps[Flavour(kf_Xi)][1]          = 3;  // d -> s
  m_maps[Flavour(kf_Xi)][2]          = 1;  // u -> d
  m_maps[Flavour(kf_Xi)][3]          = 2;  // s -> u
  m_maps[Flavour(kf_Xi_minus)]       = PDF_Id_List();
  m_maps[Flavour(kf_Xi_minus)][3]    = 2;  // s -> u
  m_maps[Flavour(kf_Xi_minus)][2]    = 3;  // u -> s
  // Pion and kaon PDFs (only K^{+-} at the moment):
  // - again we assume SU_F(3) symmetry
  // - TODO: for K_L we need to average d and s PDFs or make
  //         a decision based on wave function
  m_maps[Flavour(kf_pi_plus).Bar()]   = PDF_Id_List();
  m_maps[Flavour(kf_K_plus).Bar()]    = PDF_Id_List();
  m_maps[Flavour(kf_K_plus).Bar()][1] = 3;  // d -> s
  m_maps[Flavour(kf_K_plus).Bar()][3] = 1;  // s -> d
  // Lepton and photon PDFs
  m_maps[Flavour(kf_e)]             = PDF_Id_List();
  m_maps[Flavour(kf_mu)]            = PDF_Id_List();
  m_maps[Flavour(kf_photon).Bar()]  = PDF_Id_List();
  // Pomerons and reggeons:
  m_maps[Flavour(kf_pomeron).Bar()] = PDF_Id_List();
  m_maps[Flavour(kf_reggeon).Bar()] = PDF_Id_List();
}

PDF_Base::PDF_Base() :
  m_type("none"), m_member(0), m_lhef_number(-1), m_nf(-1), m_exponent(1.),
  m_rescale(1.),m_xmin(1.), m_xmax(0.), m_q2min(1.e12), m_q2max(0.),
  m_rescX(false)
{
  RegisterDefaults();
  Settings& s   = Settings::GetMainSettings();
  m_lhef_number = s["LHEF_PDF_NUMBER"].Get<int>();
  m_force_4f    = s["PDF_FORCE_4F"].Get<int>();
}

PDF_Base::~PDF_Base() = default;

void PDF_Base::RegisterDefaults()
{
  Settings& s = Settings::GetMainSettings();
  s["LHEF_PDF_NUMBER"].SetDefault(-1);
  s["PDF_FORCE_4F"].SetDefault(0);

  Scoped_Settings lhapdfsettings{ s["LHAPDF"] };
  lhapdfsettings["NUMBER_OF_FLAVOURS"].SetDefault(5);
  lhapdfsettings["GRID_PATH"].SetDefault("");
  lhapdfsettings.DeclareVectorSettingsWithEmptyDefault({ "DISALLOW_FLAVOUR" });
}

double PDF_Base::AlphaSPDF(const double &q2)
{
  THROW(not_implemented, "PDF doesn't implement alpha_s running.")
}

void PDF_Base::Calculate(double x,double Q2)
{
  if(Q2<m_q2min) {
    static double lasterr(-1.0);
    if (Q2!=lasterr)
    msg_Error()<<METHOD<<"(): Q-range violation Q = "<<sqrt(Q2)
	       <<" < "<<sqrt(m_q2min)<<". Set Q -> "
	       <<sqrt(m_q2min)<<"."<<std::endl;
    lasterr=Q2;
    Q2=1.000001*m_q2min;
  }
  if(Q2>m_q2max) {
    static double lasterr(-1.0);
    if (Q2!=lasterr)
    msg_Error()<<METHOD<<"(): Q-range violation Q = "<<sqrt(Q2)
	       <<" > "<<sqrt(m_q2max)<<". Set Q -> "
	       <<sqrt(m_q2max)<<"."<<std::endl;
    lasterr=Q2;
    Q2=0.999999*m_q2max;
  }
  double xR = x*(m_rescX?m_rescale:1.);
  if(xR<m_xmin*m_rescale) {
    static double lasterr(-1.0);
    if (xR!=lasterr)
    msg_Error()<<METHOD<<"(): x = "<<xR<<" ("<<m_rescale
	       <<") < "<<m_xmin<<". Set x -> "
	       <<m_xmin<<"."<<std::endl;
    lasterr=xR;
    xR=1.000001*m_xmin*m_rescale;
  }
  if(xR>m_xmax*m_rescale) {
    static double lasterr(-1.0);
    if (xR!=lasterr)
    msg_Error()<<METHOD<<"(): x = "<<x<<" ("<<m_rescale
	       <<") > "<<m_xmax<<". Set x -> "
	       <<m_xmax<<"."<<std::endl;
    lasterr=xR;
    xR=0.999999*m_xmax*m_rescale;
  }
  return CalculateSpec(xR,Q2);
}

void PDF_Base::ShowSyntax(const size_t i)
{
  if (!msg_LevelIsInfo() || i==0) return;
  msg_Out()<<METHOD<<"(): {\n\n"
	   <<"   // available PDF sets ...\n"
	   <<"   // specified by PDF_SET: <set for both beams>\n"
	   <<"   // or PDF_SET: [<set for beam_1>, <set for beam_2>])\n"
	   <<"   // Default can be used as a placeholder to let Sherpa choose\n\n";
  PDF_Getter_Function::PrintGetterInfo(msg->Out(),25);
  msg_Out()<<"\n}"<<std::endl;
}


bool PDF_Base::Contains(const Flavour &flav) const
{
  if (m_force_4f && (flav.Kfcode()==kf_b || flav.Kfcode()==kf_t)) return false;

  if (m_partons.find(flav) != m_partons.end()) return true;
  for (size_t i(0); i < flav.Size(); ++i)
    if (m_partons.find(flav[i]) == m_partons.end()) return false;
  return true;
}

void PDF_Base::FillHisto(const Flavour & flav) {
  msg_Out()<<"Parton map for "<<flav<<": ";
  for (size_t i=1;i<4;i++) msg_Out()<<Flavour(i)<<" --> "<<Flavour(m_kfmap[i])<<" ";
  msg_Out()<<"\n";
  double Q2 = 10., x=1.e-4, nsteps=10.;
  for (size_t i=0;i<int(nsteps);i++) {
    Calculate(x,Q2);
    msg_Out()<<"PDF("<<flav<<", "
	     <<"x = "<<std::setprecision(6)<<std::setw(8)<<x<<", "
	     <<"Q2 = "<<std::setprecision(6)<<std::setw(8)<<Q2<<"): "
	     <<"u = "<<std::setprecision(6)<<std::setw(8)<<GetXPDF(Flavour(kf_u))<<", "
	     <<"ubar = "<<std::setprecision(6)<<std::setw(8)<<GetXPDF(Flavour(kf_u).Bar())<<", "
	     <<"d = "<<std::setprecision(6)<<std::setw(8)<<GetXPDF(Flavour(kf_d))<<", "
	     <<"dbar = "<<std::setprecision(6)<<std::setw(8)<<GetXPDF(Flavour(kf_d).Bar())<<", "
	     <<"s = "<<std::setprecision(6)<<std::setw(8)<<GetXPDF(Flavour(kf_s))<<", "
	     <<"sbar = "<<std::setprecision(6)<<std::setw(8)<<GetXPDF(Flavour(kf_s).Bar())
	     <<"\n";
    x *= exp((log(0.99) - log(1.e-4))/nsteps);
  }
}
