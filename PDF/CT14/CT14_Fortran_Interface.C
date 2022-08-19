#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include <unistd.h>

#include "PDF/Main/PDF_Base.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/My_MPI.H"

using namespace PDF;
using namespace ATOOLS;


extern "C" {
  void    setct14_(char *);
  double  ct14pdf_(int &,double &, double &);
  void    shabrt_() { Abort(); }
}


namespace PDF {

  class CT14_Fortran_Interface : public PDF_Base {
  private:
    double      m_f[11], m_x, m_Q;
    bool        m_calculated[11];

  public:

    CT14_Fortran_Interface(const ATOOLS::Flavour bunch,
                           std::string set, int member)
    {
      m_xmin=1.e-8;
      m_xmax=1.;
      m_q2min=1.69;
      m_q2max=1.e10;

      m_set=set;
      m_type=m_set;
      m_bunch=bunch;
      m_member=member;
      std::string cset("");
      std::string path = rpa->gen.Variable("SHERPA_SHARE_PATH")+"/CT14Grid";

      std::string num;
      if (m_member<10) num="0"+ToString(m_member);
      else             num=ToString(m_member);
      std::string asmz[13] = {"0.111", "0.112", "0.113", "0.114",
                              "0.115", "0.116", "0.117", "0.118", "0.119",
                              "0.120", "0.121", "0.122", "0.123"};
      if (m_set==std::string("ct14nn")) {
        cset = std::string("ct14nn."+num+".pds");
        m_asinfo.m_order=2;
        m_asinfo.m_nf=5;
        m_asinfo.m_asmz=0.118;
        m_lhef_number=13000+m_member;
      }
      for (size_t i=0; i<13; ++i) {
        if (m_set==std::string("ct14nn.as"+asmz[i]) && m_member==0) {
          cset = std::string("ct14nn.as"+asmz[i]+".pds");
          m_asinfo.m_order=2;
          m_asinfo.m_nf=5;
          m_asinfo.m_asmz=ToType<double>(asmz[i]);
          m_lhef_number=13060+i;
        }
      }
      if (m_set==std::string("ct14n")) {
        cset = std::string("ct14n."+num+".pds");
        m_asinfo.m_order=1;
        m_asinfo.m_nf=5;
        m_asinfo.m_asmz=0.118;
        m_lhef_number=13100+m_member;
      }
      for (size_t i=0; i<13; ++i) {
        if (m_set==std::string("ct14n.as"+asmz[i]) && m_member==0) {
          cset = std::string("ct14n.as"+asmz[i]+".pds");
          m_asinfo.m_order=1;
          m_asinfo.m_nf=5;
          m_asinfo.m_asmz=ToType<double>(asmz[i]);
          m_lhef_number=13158+i;
        }
      }
      if (m_set==std::string("ct14ll")) {
        cset = std::string("ct14ll.pds");
        m_asinfo.m_order=1;
        m_asinfo.m_nf=5;
        m_asinfo.m_asmz=0.130;
        m_lhef_number=13205;
      }
      m_asinfo.m_mz2=sqr(91.1876);
      m_nf=m_asinfo.m_nf;

      if (cset=="") {
        THROW(fatal_error,"PDF set "+m_set
                          +" ("+ToString(m_member)+") not found.");
      }

      char buffer[1024];
      char * err = getcwd(buffer,1024);
      if (err==NULL) {
        msg_Error()<<"Error in CT14_Fortran_Interface.C "<<std::endl;
      }
      int stat=chdir(path.c_str());
      msg_Tracking()<<METHOD<<"(): Init cset "<<cset<<"."<<std::endl;
      char tablefile[40];
      MakeFortranString(tablefile,cset,40);
      setct14_(tablefile);
      if (stat==0) {
        chdir(buffer);
      }
      else {
        msg_Error()<<"Error in CT14_Fortran_Interface.C "<<std::endl
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
      PDF_Base *copy = new CT14_Fortran_Interface(m_bunch,m_set,m_member);
      m_copies.push_back(copy);
      return copy;
    }


    void   CalculateSpec(const double& x, const double& Q2) {
      for (size_t i=0;i<11;++i) m_calculated[i]=false;
      m_x=x/m_rescale;
      m_Q=sqrt(Q2);
    }


    double GetXPDF(const ATOOLS::Flavour& infl) {
      if (m_x>m_xmax || m_rescale<0.) return 0.;
      if (!(m_x>=0.0 && m_x<=1.0)) {
        PRINT_INFO("PDF called with x="<<m_x);
        return 0.;
      }
      int cteqindex;
      switch (infl.Kfcode()) {
      case kf_gluon: cteqindex=0;                  break;
      case kf_d:     cteqindex=(m_bunch.IsAnti()?-1:1)*int(infl)*2; break;
      case kf_u:     cteqindex=(m_bunch.IsAnti()?-1:1)*int(infl)/2; break;
      default:                cteqindex=(m_bunch.IsAnti()?-1:1)*int(infl);   break;
      }
      if (!m_calculated[5-cteqindex]) {
        m_f[5-cteqindex]=ct14pdf_(cteqindex,m_x,m_Q)*m_x;
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
      case kf_d:     cteqindex=(m_bunch.IsAnti()?-1:1)*(anti?-2:2);   break;
      case kf_u:     cteqindex=(m_bunch.IsAnti()?-1:1)*(anti?-1:1);   break;
      default:       cteqindex=(m_bunch.IsAnti()?-1:1)*(anti?-kf:kf); break;
      }
      if (!m_calculated[5-cteqindex]) {
        m_f[5-cteqindex]=ct14pdf_(cteqindex,m_x,m_Q)*m_x;
        m_calculated[5-cteqindex]=true;
      }
      return m_rescale*m_f[5-cteqindex];
    }

    inline void MakeFortranString(char *output,std::string input,
                                  unsigned int length)
    {
      for (unsigned int i=0;i<length;++i) output[i]=(char)32;
      for (size_t j=0;j<input.length();++j) output[j]=(char)input[j];
    }

  };

}


DECLARE_PDF_GETTER(CT14_Getter);

PDF_Base *CT14_Getter::operator()
  (const Parameter_Type &args) const
{
  if (!args.m_bunch.IsHadron()) return NULL;
  return new CT14_Fortran_Interface(args.m_bunch,args.m_set,args.m_member);
}

void CT14_Getter::PrintInfo
(std::ostream &str,const size_t width) const
{
  str<<"CT14 fit, see arXiv:1506.07443 [hep-ph]";
}


CT14_Getter *p_get_ct14[30];
extern "C" void InitPDFLib()
{
  p_get_ct14[0]  = new CT14_Getter("ct14nn");
  p_get_ct14[1]  = new CT14_Getter("ct14n");
  p_get_ct14[2]  = new CT14_Getter("ct14ll");

  std::string asmz[13] = {"0.111", "0.112", "0.113", "0.114",
                          "0.115", "0.116", "0.117", "0.118", "0.119",
                          "0.120", "0.121", "0.122", "0.123"};
  for (size_t i(0);i<13;++i) {
    p_get_ct14[3+i] = new CT14_Getter("ct14nnlo.as"+asmz[i]);
    if (i==0 || i==1 || i==18 || i==19 || i==20) continue;
    p_get_ct14[16+i] = new CT14_Getter("ct14n.as"+asmz[i]);
  }
}

extern "C" void ExitPDFLib()
{
  for (int i(0);i<30;++i) delete p_get_ct14[i];
}
