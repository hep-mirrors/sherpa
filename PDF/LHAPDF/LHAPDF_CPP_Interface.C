#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Math/Random.H"
#include "PDF/Main/PDF_Base.H"
#include "ATOOLS/Phys/Flavour.H"

#include "LHAPDF/LHAPDF.h"

namespace PDF {
  class LHAPDF_CPP_Interface : public PDF_Base {
  private:
    LHAPDF::PDF * p_pdf;
    int           m_smember;
    std::map<int, double> m_xfx;
    std::map<int, bool>   m_calculated;
    double        m_x,m_Q2;
    std::vector<int> m_disallowedflavour;
  public:
    LHAPDF_CPP_Interface(const ATOOLS::Flavour,std::string,int);
    ~LHAPDF_CPP_Interface();
    PDF_Base * GetCopy();

    void   CalculateSpec(const double&,const double&);
    double AlphaSPDF(const double &);
    double GetXPDF(const ATOOLS::Flavour&);
    double GetXPDF(const kf_code&, bool);
    double GetDefaultAlpha();
    double GetDefaultScale();
    int GetFlavourScheme();
    
    void SetAlphaSInfo();
    void SetPDFMember();

  };
}

using namespace PDF;
using namespace ATOOLS;


LHAPDF_CPP_Interface::LHAPDF_CPP_Interface(const ATOOLS::Flavour _bunch,
                                           const std::string _set,
                                           const int _member) :
  p_pdf(NULL)
{
  m_set=_set;
  m_smember=_member;
  m_type="LHA["+m_set+"]";

  Scoped_Settings s{ Settings::GetMainSettings()["LHAPDF"] };

  m_bunch = _bunch;
  static std::set<std::string> s_init;
  if (s_init.find(m_set)==s_init.end()) {
    m_member=abs(m_smember);
    int lhapdfverb(LHAPDF::verbosity());
    LHAPDF::setVerbosity(msg_LevelIsDebugging()?lhapdfverb:0);
    p_pdf = LHAPDF::mkPDF(m_set,m_smember);
    LHAPDF::setVerbosity(lhapdfverb);
    SetAlphaSInfo();
  }

  auto q2lim = s["USE_Q2LIMIT"].SetDefault(1).Get<int>();
  // get x,Q2 ranges from PDF
  m_xmin=p_pdf->xMin();
  m_xmax=p_pdf->xMax();
  m_q2min=q2lim?p_pdf->q2Min():0.0;
  m_q2max=q2lim?p_pdf->q2Max():1.0e37;
  m_nf=m_asinfo.m_nf;

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
  kfcs.push_back(kf_e);
  kfcs.push_back(-kf_e);
  kfcs.push_back(kf_mu);
  kfcs.push_back(-kf_mu);
  kfcs.push_back(kf_tau);
  kfcs.push_back(-kf_tau);
  kfcs.push_back(kf_gluon);
  kfcs.push_back(kf_photon);
  kfcs.push_back(kf_Z);
  kfcs.push_back(kf_Wplus);
  kfcs.push_back(-kf_Wplus);
  kfcs.push_back(kf_h0);
  for (int i=0;i<kfcs.size();i++) if (p_pdf->hasFlavor(kfcs[i])) {
    m_partons.insert(Flavour(abs(kfcs[i]),kfcs[i]<0));
    m_xfx[kfcs[i]]=0.;
    m_calculated[kfcs[i]]=false;
  }
  if (p_pdf->hasFlavor(kf_d)) {
    m_partons.insert(Flavour(kf_quark));
    if (p_pdf->hasFlavor(kf_gluon)) {
      m_partons.insert(Flavour(kf_jet));
      if (p_pdf->hasFlavor(kf_photon)) {
        m_partons.insert(Flavour(kf_ewjet));
      }
    }
  }

  m_lhef_number = p_pdf->lhapdfID();

  m_disallowedflavour = s["DISALLOW_FLAVOUR"].GetVector<int>();
  if (m_disallowedflavour.size()) {
    msg_Info()<<METHOD<<"(): Set PDF for the following flavours to zero: ";
    for (size_t i(0);i<m_disallowedflavour.size();++i)
      msg_Info()<<Flavour(abs(m_disallowedflavour[i]),m_disallowedflavour[i]<0)
                <<" ";
    msg_Info()<<std::endl;
  }

  rpa->gen.AddCitation(1,"LHAPDF6 is published under \\cite{Buckley:2014ana}.");
}

double LHAPDF_CPP_Interface::GetDefaultAlpha()
{
 return m_asinfo.m_asmz;
}

int LHAPDF_CPP_Interface::GetFlavourScheme()
{
 int nflav=p_pdf->info().get_entry_as<int>("NumFlavors");
 if (p_pdf->info().get_entry_as<std::string>("FlavorScheme")=="variable") {
  if (nflav==6){
   return -1;
  }
  nflav+=10;
 }
 return nflav;
}

double LHAPDF_CPP_Interface::GetDefaultScale()
{
 return m_asinfo.m_mz2;
}

void LHAPDF_CPP_Interface::SetAlphaSInfo()
{
  if (m_asinfo.m_order>=0) return;
  // TODO: get alphaS info
  m_asinfo.m_order=p_pdf->info().get_entry_as<int>("AlphaS_OrderQCD");
  m_asinfo.m_nf=p_pdf->info().get_entry_as<int>("NumFlavors",-1);
  if (m_asinfo.m_nf<0) {
    Scoped_Settings s{ Settings::GetMainSettings()["LHAPDF"] };
    const int nf(s["NUMBER_OF_FLAVOURS"].Get<int>());
    msg_Info()<<METHOD<<"(): No nf info. Set nf = "<<nf<<"\n";
    m_asinfo.m_nf=nf;
  }
  m_asinfo.m_flavs.resize(m_asinfo.m_nf);
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

LHAPDF_CPP_Interface::~LHAPDF_CPP_Interface()
{
  if (p_pdf) { delete p_pdf; p_pdf=NULL; }
}


PDF_Base * LHAPDF_CPP_Interface::GetCopy() 
{
  return new LHAPDF_CPP_Interface(m_bunch,m_set,m_smember);
}

double LHAPDF_CPP_Interface::AlphaSPDF(const double &scale2) {
  if (IsBad(scale2) || scale2<0.0) {
    msg_Error()<<METHOD<<"(): t = "<<scale2<<". Returning zero."<<std::endl;
    return 0.0;
  }
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
  m_x=x/m_rescale;
  m_Q2=Q2;
}

double LHAPDF_CPP_Interface::GetXPDF(const ATOOLS::Flavour& infl) {
  return GetXPDF(infl.Kfcode(), infl.IsAnti());
}

double LHAPDF_CPP_Interface::GetXPDF(const kf_code& kf, bool anti) {
  if (IsBad(m_x) || IsBad(m_Q2) || m_Q2<0.0) {
    msg_Error()<<METHOD<<"(): Encountered bad (x,Q2)=("<<m_x<<","<<m_Q2<<"), "
                       <<"returning zero."<<std::endl;
    return 0.;
  }
  int kfc = (m_bunch.IsAnti()?-1:1)*(anti?-kf:kf);
  if (kf==kf_gluon || kf==kf_photon)
    kfc = kf;
  for (size_t i(0);i<m_disallowedflavour.size();++i) {
    if (kfc==m_disallowedflavour[i]) {
      m_xfx[kfc]=0.;
      m_calculated[kfc]=true;
      break;
    }
  }
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
  if (!args.m_bunch.IsHadron() && !args.m_bunch.IsPhoton()) return NULL;
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
  Scoped_Settings s{ Settings::GetMainSettings()["LHAPDF"] };
  if (s["GRID_PATH"].IsSetExplicitly())
    LHAPDF::setPaths(s["GRID_PATH"].Get<std::string>());
  const std::vector<std::string>& sets(LHAPDF::availablePDFSets());
  msg_Debugging()<<METHOD<<"(): LHAPDF paths: "<<LHAPDF::paths()<<std::endl;
  msg_Debugging()<<METHOD<<"(): LHAPDF sets: "<<sets<<std::endl;
  for (size_t i(0);i<sets.size();++i)
    p_get_lhapdf.push_back(new LHAPDF_Getter(sets[i]));
}

extern "C" void ExitPDFLib()
{
  // prevent LHAPDF citation info from appearing
  LHAPDF::setVerbosity(0);
  for (size_t i(0);i<p_get_lhapdf.size();++i) delete p_get_lhapdf[i];
}
