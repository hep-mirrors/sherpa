#include "PDF/LHAPDF/LHAPDF_Fortran_Interface.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/Library_Loader.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Math/Random.H"
#include <cstring>
#include <dirent.h>
#include <cstring>

#ifdef DARWIN
#if __MAC_OS_X_VERSION_MIN_REQUIRED < __MAC_10_8
#define DIRENT_TYPE dirent
#else
#define DIRENT_TYPE const dirent
#endif // Lion or earlier
#else //DARWIN
#define DIRENT_TYPE const dirent
#endif //DARWIN

#ifndef _D_EXACT_NAMLEN
#define _D_EXACT_NAMLEN(ENTRY) ENTRY->d_namlen
#endif

using namespace PDF;
using namespace ATOOLS;

#include "LHAPDF/LHAPDF.h"

LHAPDF_Fortran_Interface::LHAPDF_Fortran_Interface(const ATOOLS::Flavour _bunch,
						   const std::string _set,const int _member) :
  m_anti(1)
{
  m_set=set;
  m_smember=_member;
  m_type="LHA["+m_set+"]";

  m_bunch = _bunch; 
  if (m_bunch==Flavour(kf_p_plus).Bar()) m_anti=-1;
  static std::set<std::string> s_init;
  if (s_init.find(m_set)==s_init.end()) {
    if (m_smember!=0) msg_Info()<<METHOD<<"(): Init member "<<m_smember<<"."<<std::endl;
    m_member=abs(m_smember);
    LHAPDF::initPDFSet(m_set);
    LHAPDF::initPDF(m_member);
    m_asinfo.m_order=LHAPDF::getOrderAlphaS();
    if (LHAPDF::getNf()<0) {
      Data_Reader read(" ",";","#","=");
      int nf(read.GetValue<int>("LHAPDF_NUMBER_OF_FLAVOURS",5));
      msg_Info()<<METHOD<<"(): No nf info. Set nf = "<<nf<<"\n";
      m_asinfo.m_flavs.resize(nf);
    }
    else m_asinfo.m_flavs.resize(LHAPDF::getNf());
    for (size_t i(0);i<m_asinfo.m_flavs.size();++i) {
      m_asinfo.m_flavs[i]=PDF_Flavour((kf_code)i+1);
      m_asinfo.m_flavs[i].m_mass=LHAPDF::getQMass(i+1);
      m_asinfo.m_flavs[i].m_thres=LHAPDF::getThreshold(i+1);
    }
    m_asinfo.m_asmz=AlphaSPDF(sqr(Flavour(kf_Z).Mass()));
  }

  m_xmin=LHAPDF::getXmin(m_member);
  m_xmax=LHAPDF::getXmax(m_member);
  m_q2min=LHAPDF::getQ2min(m_member);
  m_q2max=LHAPDF::getQ2max(m_member);
  
  for (int i=1;i<=m_asinfo.m_flavs.size();i++) {
    m_partons.insert(Flavour((kf_code)(i)));
    m_partons.insert(Flavour((kf_code)(i)).Bar());
  }
  m_partons.insert(Flavour(kf_gluon));
  m_partons.insert(Flavour(kf_jet));
  m_partons.insert(Flavour(kf_quark));
  m_partons.insert(Flavour(kf_quark).Bar());
  if (LHAPDF::hasPhoton()) m_partons.insert(Flavour(kf_photon));

//  m_lhef_number = LHAPDF::getPDFSetInfo(m_set,m_member).id;
}

PDF_Base * LHAPDF_Fortran_Interface::GetCopy() 
{
  return new LHAPDF_Fortran_Interface(m_bunch,m_set,m_smember);
}


double LHAPDF_Fortran_Interface::AlphaSPDF(const double &scale2) {
  double scale = sqrt(scale2);
  double as    = LHAPDF::alphasPDF(scale);
  return as;
}

void LHAPDF_Fortran_Interface::SetPDFMember()
{
  if (m_smember<0) {
    double rn=ran->Get();
    m_member=1+Min((int)(rn*abs(m_smember)),-m_smember-1);
    LHAPDF::initPDF(m_member);
  }
}

void LHAPDF_Fortran_Interface::CalculateSpec(double x,double Q2) {
  double Q = sqrt(Q2);
  if (LHAPDF::hasPhoton()) m_fv=LHAPDF::xfxphoton(x,Q);
  else                     m_fv=LHAPDF::xfx(x,Q);
}

double LHAPDF_Fortran_Interface::GetXPDF(const ATOOLS::Flavour infl) {
  int kfc;
  if (int(infl) == kf_gluon) {
    kfc=0;
  } else if (LHAPDF::hasPhoton() && int(infl) == kf_photon) {
    kfc=7;
  } else {
    kfc=m_anti*int(infl);
    if (kfc<-6 || kfc>6) {
      msg_Out()<<"WARNING in "<<METHOD<<"("<<infl<<") not supported by this PDF!"<<std::endl;
      return 0.;
    }
  }
  return m_rescale*m_fv[6+kfc];
}

DECLARE_PDF_GETTER(LHAPDF_Getter);

PDF_Base *LHAPDF_Getter::operator()
  (const Parameter_Type &args) const
{
  if (!args.m_bunch.IsHadron()) return NULL;
  return new LHAPDF_Fortran_Interface(args.m_bunch,args.m_set,args.m_member);
}

void LHAPDF_Getter::PrintInfo
(std::ostream &str,const size_t width) const
{
  str<<"LHAPDF interface";
}

int LHAPDF_DummyInclude(DIRENT_TYPE *entry)
{
  return true;
}

std::vector<std::string> LHAPDF_ScanDir(const std::string &path) 
{
  msg_Debugging()<<METHOD<<"(): Scanning directory "<<path<<" {\n";
  std::vector<std::string> res;
  struct dirent **entries;
  int n(scandir(path.c_str(),&entries,&LHAPDF_DummyInclude,alphasort));
  if (n<0) {
    msg_Error()<<METHOD<<"(): Scandir error in "<<path<<". Abort."<<std::endl;
    return res;
  }
  for (int i(0);i<n;++i) {
    bool isdir(entries[i]->d_type==DT_DIR);
#ifdef ARCH_LINUX
    struct dirent **dentries;
    int n(scandir((path+"/"+entries[i]->d_name).c_str(),
		  &dentries,&LHAPDF_DummyInclude,alphasort));
    if (n>=0) isdir=true;
#endif
    if (!isdir) {
      res.push_back(entries[i]->d_name);
      msg_Debugging()<<"  "<<i<<": "<<entries[i]->d_name<<"\n";
    }
    free(entries[i]);
  }
  free(entries);
  msg_Debugging()<<"}\n";
  return res;
}

std::vector<LHAPDF_Getter*> p_get_lhapdf;

extern "C" void InitPDFLib()
{
  Data_Reader read(" ",";","!","=");
  read.AddComment("#");
  read.AddWordSeparator("\t");
  std::string path;
  if (read.ReadFromFile(path,"LHAPDF_GRID_PATH")) LHAPDF::setPDFPath(path); 
  std::vector<std::string> files=LHAPDF_ScanDir(LHAPDF::pdfsetsPath());
  p_get_lhapdf.resize(files.size());
  for (size_t i(0);i<files.size();++i) p_get_lhapdf[i] = new LHAPDF_Getter(files[i]);
}

extern "C" void ExitPDFLib()
{
  for (size_t i(0);i<p_get_lhapdf.size();++i) delete p_get_lhapdf[i];
}
