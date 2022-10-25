#include "PDF/Main/PDF_Base.H"

#define COMPILE__Getter_Function
#define OBJECT_TYPE PDF::PDF_Base
#define PARAMETER_TYPE PDF::PDF_Arguments
#include "ATOOLS/Org/Getter_Function.C"

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
  m_deflib[kf_p_plus] = "LHAPDFSherpa";
  m_defset[kf_p_plus] = "PDF4LHC21_40_pdfas";

  m_deflib[kf_e] = "PDFESherpa";
  m_defset[kf_e] = "PDFe";

  m_deflib[kf_photon] = "GRVSherpa";
  m_defset[kf_photon] = "GRV";
}

std::ostream &PDF::operator<<(std::ostream &ostr,const PDF::PDF_AS_Info &asi)
{
  return ostr<<"\\alpha_s of order="<<asi.m_order
             <<" with \\alpha_s(\\mu="<<sqrt(asi.m_mz2)<<")="<<asi.m_asmz;
}

PDF_Base::PDF_Base():
  m_type("none"), m_member(0), m_lhef_number(-1), m_exponent(1.),
  m_rescale(1.), m_rescX(false),
  m_xmin(1.), m_xmax(0.), m_q2min(1.e12), m_q2max(0.),
  m_nf(-1)
{
  RegisterDefaults();
  Settings& s = Settings::GetMainSettings();
  m_lhef_number = s["LHEF_PDF_NUMBER"].Get<int>();
}

PDF_Base::~PDF_Base() = default;

void PDF_Base::RegisterDefaults()
{
  Settings& s = Settings::GetMainSettings();
  s["LHEF_PDF_NUMBER"].SetDefault(-1);
  s["INCLUDE_PHOTON_IN_PHOTON_PDF"].SetDefault(false);

  Scoped_Settings lhapdfsettings{ s["LHAPDF"] };
  lhapdfsettings["NUMBER_OF_FLAVOURS"].SetDefault(5);
  lhapdfsettings["GRID_PATH"].SetDefault("");
  lhapdfsettings.DeclareVectorSettingsWithEmptyDefault({ "DISALLOW_FLAVOUR" });
}

double PDF_Base::GetDefaultAlpha()
{
 return -1.0;
}

double PDF_Base::GetDefaultScale()
{
 return -1.0;
}

int PDF_Base::GetFlavourScheme()
{
 return 0;
}

double PDF_Base::AlphaSPDF(const double &q2)
{
  THROW(not_implemented, "PDF doesn't implement alpha_s running.");
  return -1.0;
}

PDF_Base *PDF_Base::GetBasicPDF() 
{
  return this;
}

void PDF_Base::SetBounds()
{
  m_rq2min=m_q2min;
  m_rq2max=m_q2max;
}

void PDF_Base::SetAlphaSInfo()
{
}

void PDF_Base::SetPDFMember()
{
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
