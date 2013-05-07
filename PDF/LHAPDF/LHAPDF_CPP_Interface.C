#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/Library_Loader.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Math/Random.H"
#include <cstring>
#include <cstring>

#include "LHAPDF/LHAPDF.h"
#include "PDF/Main/PDF_Base.H"
#include "ATOOLS/Phys/Flavour.H"

namespace PDF {
  class LHAPDF_CPP_Interface : public PDF_Base {
  private:
    LHAPDF::PDF * p_pdf;
    std::string   m_set;
    int           m_smember;
    int           m_anti;
    std::map<int, double> m_f;
    std::map<int, double> m_fv;
  public:
    LHAPDF_CPP_Interface(const ATOOLS::Flavour,std::string,int);
    ~LHAPDF_CPP_Interface();
    PDF_Base * GetCopy();

    void   CalculateSpec(double,double);
    double AlphaSPDF(const double &);
    double GetXPDF(const ATOOLS::Flavour);

    void SetPDFMember();

  };
}

using namespace PDF;
using namespace ATOOLS;
using namespace LHAPDF;


LHAPDF_CPP_Interface::LHAPDF_CPP_Interface(const ATOOLS::Flavour _bunch,
                                           const std::string _set,
                                           const int _member) :
  m_set(_set), m_anti(1)
{
  m_smember=_member;
  m_type="LHA["+m_set+"]";

  m_bunch = _bunch; 
  if (m_bunch==Flavour(kf_p_plus).Bar()) m_anti=-1;
  static std::set<std::string> s_init;
  if (s_init.find(m_set)==s_init.end()) {
    if (m_smember!=0) msg_Info()<<METHOD<<"(): Init member "<<m_smember<<"."<<std::endl;
    m_member=abs(m_smember);
    p_pdf = LHAPDF::mkPDF(m_set,m_smember);
    // TODO: get alphaS info
    m_asinfo.m_order=p_pdf->info().metadata<int>("AlphaS_OrderQCD");
    int nf(p_pdf->info().metadata<int>("EvolutionNf"));
    if (nf<0) m_asinfo.m_flavs.resize(5);
    else      m_asinfo.m_flavs.resize(nf);
    // for now assume thresholds are equal to masses
    for (size_t i(0);i<m_asinfo.m_flavs.size();++i) {
      m_asinfo.m_flavs[i]=PDF_Flavour((kf_code)i+1);
      if      (i==0)
        m_asinfo.m_flavs[i].m_mass=m_asinfo.m_flavs[i].m_thres
                                  =p_pdf->info().metadata<double>("MDown");
      else if (i==1)
        m_asinfo.m_flavs[i].m_mass=m_asinfo.m_flavs[i].m_thres
                                  =p_pdf->info().metadata<double>("MUp");
      else if (i==2)
        m_asinfo.m_flavs[i].m_mass=m_asinfo.m_flavs[i].m_thres
                                  =p_pdf->info().metadata<double>("MStrange");
      else if (i==3)
        m_asinfo.m_flavs[i].m_mass=m_asinfo.m_flavs[i].m_thres
                                  =p_pdf->info().metadata<double>("MCharm");
      else if (i==4)
        m_asinfo.m_flavs[i].m_mass=m_asinfo.m_flavs[i].m_thres
                                  =p_pdf->info().metadata<double>("MBottom");
      else if (i==5)
        m_asinfo.m_flavs[i].m_mass=m_asinfo.m_flavs[i].m_thres
                                  =p_pdf->info().metadata<double>("MTop");
    }
    m_asinfo.m_asmz=p_pdf->info().metadata<double>("AlphaS_MZ");
  }

  // get x,Q2 ranges from PDF
  m_xmin=p_pdf->info().metadata<double>("XMin");
  m_xmax=p_pdf->info().metadata<double>("XMax");
  m_q2min=p_pdf->info().metadata<double>("Q2Min");
  m_q2max=p_pdf->info().metadata<double>("Q2Max");
  
  for (int i=1;i<6;i++) {
    m_partons.insert(Flavour((kf_code)(i)));
    m_partons.insert(Flavour((kf_code)(i)).Bar());
  }
  m_partons.insert(Flavour(kf_gluon));
  m_partons.insert(Flavour(kf_jet));
  m_partons.insert(Flavour(kf_quark));
  m_partons.insert(Flavour(kf_quark).Bar());

  m_lhef_number = p_pdf->lhapdfID();
}

LHAPDF_CPP_Interface::~LHAPDF_CPP_Interface()
{
  if (p_pdf) { delete p_pdf; p_pdf=NULL; }
}


PDF_Base * LHAPDF_CPP_Interface::GetCopy() 
{
  return new LHAPDF_CPP_Interface(m_bunch,m_set,m_smember);
}

double LHAPDF_CPP_Interface::AlphaSPDF(const double &scale2) {
  return p_pdf->alphaS(scale2);
}

void LHAPDF_CPP_Interface::SetPDFMember()
{
  if (m_smember<0) {
    THROW(not_implemented,"Not implemented yet.")
    double rn=ran->Get();
    m_member=1+Min((int)(rn*abs(m_smember)),-m_smember-1);
    //p_pdf->initPDF(m_member);
  }
}

void LHAPDF_CPP_Interface::CalculateSpec(double x,double Q2) {
  x/=m_rescale;
  m_fv=p_pdf->xfxQ2(x,Q2);
}

double LHAPDF_CPP_Interface::GetXPDF(const ATOOLS::Flavour infl) {
  int kfc = m_anti*int(infl);
  return m_rescale*m_fv[kfc];
}

DECLARE_PDF_GETTER(LHAPDF_Getter);

PDF_Base *LHAPDF_Getter::operator()
  (const Parameter_Type &args) const
{
  if (!args.m_bunch.IsHadron()) return NULL;
  int mode=args.p_read->GetValue<int>("PDF_SET_VERSION",0);
  int ibeam=args.m_ibeam;
  mode=args.p_read->GetValue<int>("PDF_SET_VERSION_"+ToString(ibeam+1),mode);
  return new LHAPDF_CPP_Interface(args.m_bunch,m_key,mode);
}

void LHAPDF_Getter::PrintInfo
(std::ostream &str,const size_t width) const
{
  str<<"LHAPDF interface";
}

std::vector<LHAPDF_Getter*> p_get_lhapdf;

extern "C" void InitPDFLib()
{
  p_get_lhapdf.resize(1);
  p_get_lhapdf[0] = new LHAPDF_Getter("CT10nlo");
//  Data_Reader read(" ",";","!","=");
//  read.AddComment("#");
//  read.AddWordSeparator("\t");
//  std::string path;
//  if (read.ReadFromFile(path,"LHAPDF_DATA_PATH"))
//    path=std::string(LHAPDF_PATH)+"share/LHAPDF/"; 
//  std::vector<std::string> files=LHAPDF_ScanDir(LHAPDF::pdfsetsPath());
//  p_get_lhapdf.resize(files.size());
//  for (size_t i(0);i<files.size();++i) p_get_lhapdf[i] = new LHAPDF_Getter(files[i]);
}

extern "C" void ExitPDFLib()
{
  for (size_t i(0);i<p_get_lhapdf.size();++i) delete p_get_lhapdf[i];
}
