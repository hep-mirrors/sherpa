#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Data_Reader.H"
#include <unistd.h> 

#include "PDF/Main/PDF_Base.H"
#include "ATOOLS/Phys/Flavour.H"

using namespace PDF;
using namespace ATOOLS;


extern "C" {
  void    setct10_(int &);
  double  ct10pdf_(int &,double &, double &);
  void    shabrt_() { abort(); }
}


namespace PDF {

  class CT10_Fortran_Interface : public PDF_Base {
  private:
    int         m_anti;
    double      m_f[11], m_x, m_Q;
    bool        m_calculated[11];

  public:

    CT10_Fortran_Interface(const ATOOLS::Flavour bunch,
                           std::string set, int member) :
      m_anti(1)
    {
      m_xmin=1.e-8;
      m_xmax=1.;
      m_q2min=1.69;
      m_q2max=1.e10;

      m_set=set;
      m_type=m_set;
      m_bunch=bunch;
      m_member=member;
      if (m_bunch==Flavour(kf_p_plus).Bar()) m_anti=-1;
      int iset = 0;
      std::string path = rpa->gen.Variable("SHERPA_SHARE_PATH")+"/CT10Grid";
  
      if (m_set==std::string("ct10")) {
        iset = 100+m_member;
        m_asinfo.m_order=1;
        m_asinfo.m_asmz=0.118;
        if (m_member==0) m_lhef_number=10800;
        else m_lhef_number=10801;
      }
      double asmz[10] = {0.116, 0.117, 0.119, 0.120, 0.113, 0.114, 0.115, 0.121, 0.122, 0.123};
      for (size_t i=0; i<10; ++i) {
        if (m_set==std::string("ct10.as"+ToString(i)) && m_member==0) {
          iset = 10+i;
          m_asinfo.m_order=1;
          m_asinfo.m_asmz=asmz[i];
        }
      }
      if (m_set==std::string("ct10.3f") && m_member==0) {
        iset = 30;
        m_asinfo.m_order=1;
        m_asinfo.m_asmz=0.118;
      }
      if (m_set==std::string("ct10.4f") && m_member==0) {
        iset = 31;
        m_asinfo.m_order=1;
        m_asinfo.m_asmz=0.118;
      }

      if (m_set==std::string("ct10w")) {
        iset = 200+m_member;
        m_asinfo.m_order=1;
        m_asinfo.m_asmz=0.118;
      }
      for (size_t i=0; i<10; ++i) {
        if (m_set==std::string("ct10was"+ToString(i)) && m_member==0) {
          iset = 20+i;
          m_asinfo.m_order=1;
          m_asinfo.m_asmz=asmz[i];
        }
      }
      if (m_set==std::string("ct10w3f") && m_member==0) {
        iset = 32;
        m_asinfo.m_order=1;
        m_asinfo.m_asmz=0.118;
      }
      if (m_set==std::string("ct10w4f") && m_member==0) {
        iset = 33;
        m_asinfo.m_order=1;
        m_asinfo.m_asmz=0.118;
      }

      if (iset==0) {
        THROW(fatal_error, "PDF set "+m_set+" ("+ToString(m_member)+") not found.");
      }
      m_asinfo.m_mz2=sqr(91.1876);

      char buffer[1024];
      char * err = getcwd(buffer,1024);
      if (err==NULL) {
        msg_Error()<<"Error in CT10_Fortran_Interface.C "<<std::endl;
      }
      int stat=chdir(path.c_str());
      msg_Tracking()<<METHOD<<"(): Init Iset "<<iset<<"."<<std::endl;
      setct10_(iset);
      if (stat==0) {
        stat=chdir(buffer);
      }
      else {
        msg_Error()<<"Error in CT10_Fortran_Interface.C "<<std::endl
                   <<"   path "<<path<<" not found "<<std::endl;
      }
  
      for (int i=1;i<6;i++) {
        m_partons.insert(Flavour((kf_code)(i)));
        m_partons.insert(Flavour((kf_code)(i)).Bar());
      }
      m_partons.insert(Flavour(kf_gluon));
      m_partons.insert(Flavour(kf_jet));
      m_partons.insert(Flavour(kf_jet));
      m_partons.insert(Flavour(kf_quark));
      m_partons.insert(Flavour(kf_quark).Bar());                               
    }


    PDF_Base * GetCopy() {
      PDF_Base *copy = new CT10_Fortran_Interface(m_bunch,m_set,m_member);
      m_copies.push_back(copy);
      return copy;
    }


    void   CalculateSpec(const double& x, const double& Q2) {
      for (size_t i=0;i<11;++i) m_calculated[i]=false;
      m_x=x;
      m_Q=sqrt(Q2);
    }


    double GetXPDF(const ATOOLS::Flavour& infl) {
      if (m_x>m_xmax) return 0.;
      if (!(m_x>=0.0 && m_x<=1.0)) {
        PRINT_INFO("PDF called with x="<<m_x);
        return 0.;
      }
      int cteqindex;
      switch (infl.Kfcode()) {
      case kf_gluon: cteqindex=0;                  break;
      case kf_d:     cteqindex=m_anti*int(infl)*2; break;
      case kf_u:     cteqindex=m_anti*int(infl)/2; break;
      default:                cteqindex=m_anti*int(infl);   break;
      }
      if (!m_calculated[5-cteqindex]) {
        m_f[5-cteqindex]=ct10pdf_(cteqindex,m_x,m_Q)*m_x; 
        m_calculated[5-cteqindex]=true;
      }
      return m_rescale*m_f[5-cteqindex];     
    }

    double GetXPDF(const kf_code& kf, bool anti) {
      if (m_x>m_xmax) return 0.;
      if (!(m_x>=0.0 && m_x<=1.0)) {
        PRINT_INFO("PDF called with x="<<m_x);
        return 0.;
      }
      int cteqindex;
      switch (kf) {
      case kf_gluon: cteqindex=0;                    break;
      case kf_d:     cteqindex=m_anti*(anti?-2:2);   break;
      case kf_u:     cteqindex=m_anti*(anti?-1:1);   break;
      default:       cteqindex=m_anti*(anti?-kf:kf); break;
      }
      if (!m_calculated[5-cteqindex]) {
        m_f[5-cteqindex]=ct10pdf_(cteqindex,m_x,m_Q)*m_x;
        m_calculated[5-cteqindex]=true;
      }
      return m_rescale*m_f[5-cteqindex];
    }
  };

}


DECLARE_PDF_GETTER(CT10_Getter);

PDF_Base *CT10_Getter::operator()
  (const Parameter_Type &args) const
{
  if (!args.m_bunch.IsHadron()) return NULL;
  return new CT10_Fortran_Interface(args.m_bunch,args.m_set,args.m_member);
}

void CT10_Getter::PrintInfo
(std::ostream &str,const size_t width) const
{
  str<<"CT10 fit, see arXiv:1007.2241 [hep-ph]";
}


CT10_Getter *p_get_ct10[13];
CT10_Getter *p_get_ct10w[13];
extern "C" void InitPDFLib()
{
  p_get_ct10[0] = new CT10_Getter("ct10");
  for (size_t i=0; i<10; ++i)
    p_get_ct10[1+i] = new CT10_Getter("ct10.as"+ToString(i));
  p_get_ct10[11] = new CT10_Getter("ct10.3f");
  p_get_ct10[12] = new CT10_Getter("ct10.4f");

  p_get_ct10w[0] = new CT10_Getter("ct10w00");
  for (size_t i=0; i<10; ++i)
    p_get_ct10w[1+i] = new CT10_Getter("ct10was"+ToString(i));
  p_get_ct10w[11] = new CT10_Getter("ct10w3f");
  p_get_ct10w[12] = new CT10_Getter("ct10w4f");
}

extern "C" void ExitPDFLib()
{
  for (int i(0);i<13;++i) {
    delete p_get_ct10[i];
    delete p_get_ct10w[i];
  }
}
