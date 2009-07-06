#include "PDF/LHAPDF/LHAPDF_Fortran_Interface.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/Library_Loader.H"
#include <cstring>
#include <dirent.h>
#include <cstring>

#ifdef ARCH_LINUX
#define DIRENT_TYPE const dirent
#else
#define DIRENT_TYPE dirent
#endif

#ifndef _D_EXACT_NAMLEN
#define _D_EXACT_NAMLEN(ENTRY) ENTRY->d_namlen
#endif

using namespace PDF;
using namespace ATOOLS;

#ifdef LHAPDF__NATIVE__WRAPPER
#include "LHAPDF/LHAPDF.h"
#else
extern "C" {
  void   lhapdreset_();
  void   lhapdfinitset_(const char *, int len);
  void   lhapdfinitsetbyname_(const char *, int len);
  void   lhapdfinit_(int &);
  void   lhapdfevolve_(double &,double &,double *);
  double lhapdfalphas_(double &);
  void   lhapdfgetdesc_();
}
#endif

LHAPDF_Fortran_Interface::LHAPDF_Fortran_Interface(const ATOOLS::Flavour _bunch,
						   const std::string _set,const int _member) :
  m_set(_set), m_member(_member), m_anti(1)
{
  m_type="LHA["+m_set+"]";

  m_bunch = _bunch; 
  if (m_bunch==Flavour(kf_p_plus).Bar()) m_anti=-1;
  static std::set<std::string> s_init;
  if (s_init.find(m_set)==s_init.end()) {
#ifdef LHAPDF__NATIVE__WRAPPER
    LHAPDF::initPDFByName(m_set, m_member);
#else
    std::string full = m_set;
    const char * help;
    help = full.c_str();
    lhapdfinitsetbyname_(help, strlen(help));
    lhapdfinit_(m_member);
#endif
  }

#ifdef LHAPDF__NATIVE__WRAPPER
  m_xmin=LHAPDF::getXmin(m_member);
  m_xmax=LHAPDF::getXmax(m_member);
  m_q2min=LHAPDF::getQ2min(m_member);
  m_q2max=LHAPDF::getQ2max(m_member);
#else
  m_xmin=0.;
  m_xmax=1.;
  m_q2min=1.;
  m_q2max=1.e12;
#endif
  
  for (int i=1;i<6;i++) {
    m_partons.insert(Flavour((kf_code)(i)));
    m_partons.insert(Flavour((kf_code)(i)).Bar());
  }
  m_partons.insert(Flavour(kf_gluon));
  m_partons.insert(Flavour(kf_jet));
  m_partons.insert(Flavour(kf_quark));
  m_partons.insert(Flavour(kf_quark).Bar());                               
}

PDF_Base * LHAPDF_Fortran_Interface::GetCopy() 
{
  return new LHAPDF_Fortran_Interface(m_bunch,m_set,m_member);
}


double LHAPDF_Fortran_Interface::AlphaSPDF(double scale2) {
  double scale = sqrt(scale2);
#ifdef LHAPDF__NATIVE__WRAPPER
  double as    = LHAPDF::alphasPDF(scale);
#else
  double as    = lhapdfalphas_(scale);
#endif
  return as;
}

void LHAPDF_Fortran_Interface::Calculate(double x,double Q2) {
  x/=m_rescale;
  double Q = sqrt(Q2*m_fac_scale_factor);
  if(Q*Q<m_q2min) {
    msg_Error()<<METHOD<<"(): Q-range violation Q = "<<Q
	       <<" < "<<sqrt(m_q2min)<<". Set Q -> "
	       <<sqrt(m_q2min)<<"."<<std::endl;
    Q=1.000001*sqrt(m_q2min);
  }
  if(Q*Q>m_q2max) {
    msg_Error()<<METHOD<<"(): Q-range violation Q = "<<Q
	       <<" > "<<sqrt(m_q2max)<<". Set Q -> "
	       <<sqrt(m_q2max)<<"."<<std::endl;
    Q=0.999999*sqrt(m_q2max);
  }
  if(x<m_xmin) {
    msg_Error()<<METHOD<<"(): x = "<<x<<" ("<<m_rescale
	       <<") < "<<m_xmin<<". Set x -> "
	       <<m_xmin<<"."<<std::endl;
    x=1.000001*m_xmin;
  }
  if(x>m_xmax) {
    msg_Error()<<METHOD<<"(): x = "<<x<<" ("<<m_rescale
	       <<") > "<<m_xmax<<". Set x -> "
	       <<m_xmax<<"."<<std::endl;
    x=0.999999*m_xmax;
  }
#ifdef LHAPDF__NATIVE__WRAPPER
  m_fv=LHAPDF::xfx(x,Q);
#else
  lhapdfevolve_(x,Q,m_f);
#endif
}

double LHAPDF_Fortran_Interface::GetXPDF(const ATOOLS::Flavour infl) {
  int kfc = m_anti*int(infl);
  if (infl == Flavour(kf_gluon)) kfc=0;
  if (kfc<-6 || kfc>6) {
    msg_Out()<<"WARNING in LHAPDF_Fortran_Interface::GetXPDF("<<infl<<") not supported by this PDF!"<<std::endl;
    return 0.;
  }
#ifdef LHAPDF__NATIVE__WRAPPER
  return m_rescale*m_fv[6+kfc];
#else
  return m_rescale*m_f[6+kfc];
#endif
}

DECLARE_PDF_GETTER(LHAPDF_Getter);

PDF_Base *LHAPDF_Getter::operator()
  (const Parameter_Type &args) const
{
  if (!args.m_bunch.IsHadron()) return NULL;
  int mode=args.p_read->GetValue<int>("PDF_SET_VERSION",1);
  return new LHAPDF_Fortran_Interface(args.m_bunch,m_key,mode);
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
    msg_Error()<<METHOD<<"(): Scandir error. Abort."<<std::endl;
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

std::vector<LHAPDF_Getter*> p_get;

extern "C" void InitPDFLib(const std::string &path)
{
  // redirect to correct lhapdf path here
  s_loader->AddPath(std::string(LHAPDF_PATH)+"/lib");
  s_loader->LoadLibrary("LHAPDF");
#ifdef LHAPDF__NATIVE__WRAPPER
  std::vector<std::string> files=LHAPDF_ScanDir(LHAPDF::pdfsetsPath());
#else
  char *lp=getenv("LHAPATH");
  std::vector<std::string> files=LHAPDF_ScanDir
    (lp==NULL?std::string(LHAPDF_PATH)+"/share/lhapdf/PDFsets":std::string(lp));
#endif
  p_get.resize(files.size());
  for (size_t i(0);i<files.size();++i) p_get[i] = new LHAPDF_Getter(files[i]);
}

extern "C" void ExitPDFLib()
{
  for (size_t i(0);i<p_get.size();++i) delete p_get[i];
}
