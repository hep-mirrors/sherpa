#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Math/Random.H"
#include "PDF/Main/PDF_Base.H"
#include "ATOOLS/Phys/Flavour.H"

#include "LHAPDF/LHAPDF.h"

namespace PDF {
  class LHAPDF_CPP_Interface : public PDF_Base {
  private:
    LHAPDF::PDF * p_pdf;
    int           m_smember;
    int           m_anti;
    std::map<int, double> m_xfx;
    std::map<int, bool>   m_calculated;
    double        m_x,m_Q2;
  public:
    LHAPDF_CPP_Interface(const ATOOLS::Flavour,std::string,int);
    ~LHAPDF_CPP_Interface();
    PDF_Base * GetCopy();

    void   CalculateSpec(const double&,const double&);
    double AlphaSPDF(const double &);
    double GetXPDF(const ATOOLS::Flavour&);
    double GetXPDF(const kf_code&, bool);

    void SetPDFMember();

  };
}

using namespace PDF;
using namespace ATOOLS;


LHAPDF_CPP_Interface::LHAPDF_CPP_Interface(const ATOOLS::Flavour _bunch,
                                           const std::string _set,
                                           const int _member) :
  p_pdf(NULL), m_anti(1)
{
  m_set=_set;
  m_smember=_member;
  m_type="LHA["+m_set+"]";

  m_bunch = _bunch;
  if (m_bunch==Flavour(kf_p_plus).Bar()) m_anti=-1;
  static std::set<std::string> s_init;
  if (s_init.find(m_set)==s_init.end()) {
    m_member=abs(m_smember);
    p_pdf = LHAPDF::mkPDF(m_set,m_smember);
    // TODO: get alphaS info
    m_asinfo.m_order=p_pdf->info().get_entry_as<int>("AlphaS_OrderQCD");
    int nf(p_pdf->info().get_entry_as<int>("NumFlavors"));
    if (nf<0) {
      Data_Reader read(" ",";","#","=");
      int nf(read.GetValue<int>("LHAPDF_NUMBER_OF_FLAVOURS",5));
      msg_Info()<<METHOD<<"(): No nf info. Set nf = "<<nf<<"\n";
      m_asinfo.m_flavs.resize(nf);
    }
    else      m_asinfo.m_flavs.resize(nf);
    // for now assume thresholds are equal to masses, as does LHAPDF-6.0.0
    for (size_t i(0);i<m_asinfo.m_flavs.size();++i) {
      m_asinfo.m_flavs[i]=PDF_Flavour((kf_code)i+1);
      if      (i==0)
        m_asinfo.m_flavs[i].m_mass=m_asinfo.m_flavs[i].m_thres
            =p_pdf->info().get_entry_as<double>("MDown");
      else if (i==1)
        m_asinfo.m_flavs[i].m_mass=m_asinfo.m_flavs[i].m_thres
            =p_pdf->info().get_entry_as<double>("MUp");
      else if (i==2)
        m_asinfo.m_flavs[i].m_mass=m_asinfo.m_flavs[i].m_thres
            =p_pdf->info().get_entry_as<double>("MStrange");
      else if (i==3)
        m_asinfo.m_flavs[i].m_mass=m_asinfo.m_flavs[i].m_thres
            =p_pdf->info().get_entry_as<double>("MCharm");
      else if (i==4)
        m_asinfo.m_flavs[i].m_mass=m_asinfo.m_flavs[i].m_thres
            =p_pdf->info().get_entry_as<double>("MBottom");
      else if (i==5)
        m_asinfo.m_flavs[i].m_mass=m_asinfo.m_flavs[i].m_thres
            =p_pdf->info().get_entry_as<double>("MTop");
    }
    m_asinfo.m_asmz=p_pdf->info().get_entry_as<double>("AlphaS_MZ");
    m_asinfo.m_mz2=sqr(p_pdf->info().get_entry_as<double>("MZ"));
  }

  // get x,Q2 ranges from PDF
  m_xmin=p_pdf->xMin();
  m_xmax=p_pdf->xMax();
  m_q2min=p_pdf->q2Min();
  m_q2max=p_pdf->q2Max();

  // initialise all book-keep arrays etc.
  std::vector<int> kfcs;
  kfcs.push_back(kf_d);
  kfcs.push_back(-kf_d);
  kfcs.push_back(kf_u);
  kfcs.push_back(-kf_u);
  kfcs.push_back(kf_s);
  kfcs.push_back(-kf_s);
  kfcs.push_back(kf_c);
  kfcs.push_back(-kf_c);
  kfcs.push_back(kf_b);
  kfcs.push_back(-kf_b);
  kfcs.push_back(kf_t);
  kfcs.push_back(-kf_t);
  kfcs.push_back(kf_gluon);
  kfcs.push_back(kf_photon);
  for (int i=0;i<kfcs.size();i++) if (p_pdf->hasFlavor(kfcs[i])) {
    m_partons.insert(Flavour(abs(kfcs[i]),kfcs[i]<0));
    m_xfx[kfcs[i]]=0.;
    m_calculated[kfcs[i]]=false;
  }
  if (p_pdf->hasFlavor( kf_d))    m_partons.insert(Flavour(kf_quark));
  if (p_pdf->hasFlavor(-kf_d))    m_partons.insert(Flavour(kf_quark).Bar());
  if (p_pdf->hasFlavor(kf_gluon)) m_partons.insert(Flavour(kf_jet));

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
  return p_pdf->alphasQ2(scale2);
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

void LHAPDF_CPP_Interface::CalculateSpec(const double& x,const double& Q2) {
  for (std::map<int,bool>::iterator it=m_calculated.begin();
       it!=m_calculated.end();++it) it->second=false;
  m_x=x;
  m_Q2=Q2;
}

double LHAPDF_CPP_Interface::GetXPDF(const ATOOLS::Flavour& infl) {
  int kfc = m_anti*int(infl);
  if (int(infl)==kf_gluon || int(infl)==kf_photon)
    kfc = int(infl);
  if (!m_calculated[kfc]) {
    m_xfx[kfc]=p_pdf->xfxQ2(kfc,m_x,m_Q2);
    m_calculated[kfc]=true;
  }
  return m_rescale*m_xfx[kfc];
}

double LHAPDF_CPP_Interface::GetXPDF(const kf_code& kf, bool anti) {
  int kfc = m_anti*(anti?-kf:kf);
  if (kf==kf_gluon || kf==kf_photon)
    kfc = kf;
  if (!m_calculated[kfc]) {
    m_xfx[kfc]=p_pdf->xfxQ2(kfc,m_x,m_Q2);
    m_calculated[kfc]=true;
  }
  return m_rescale*m_xfx[kfc];
}

DECLARE_PDF_GETTER(LHAPDF_Getter);

PDF_Base *LHAPDF_Getter::operator()
  (const Parameter_Type &args) const
{
  if (!args.m_bunch.IsHadron()) return NULL;
  return new LHAPDF_CPP_Interface(args.m_bunch,args.m_set,args.m_member);
}

void LHAPDF_Getter::PrintInfo
(std::ostream &str,const size_t width) const
{
  str<<"LHAPDF interface";
}

std::vector<LHAPDF_Getter*> p_get_lhapdf;

extern "C" void InitPDFLib()
{
  Data_Reader read(" ",";","!","=");
  read.AddComment("#");
  read.AddWordSeparator("\t");
  std::string path;
  if (read.ReadFromFile(path,"LHAPDF_GRID_PATH")) LHAPDF::setPaths(path);
  const std::vector<std::string>& sets(LHAPDF::availablePDFSets());
  msg_Debugging()<<METHOD<<"(): LHAPDF paths: "<<LHAPDF::paths()<<std::endl;
  msg_Debugging()<<METHOD<<"(): LHAPDF sets: "<<sets<<std::endl;
  for (size_t i(0);i<sets.size();++i)
    p_get_lhapdf.push_back(new LHAPDF_Getter(sets[i]));
}

extern "C" void ExitPDFLib()
{
  for (size_t i(0);i<p_get_lhapdf.size();++i) delete p_get_lhapdf[i];
}
