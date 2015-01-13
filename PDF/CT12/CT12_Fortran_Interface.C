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
  void    setct12_(char *);
  double  ct12pdf_(int &,double &, double &);
  void    shabrt_() { abort(); }
}


namespace PDF {

  class CT12_Fortran_Interface : public PDF_Base {
  private:
    int         m_anti;
    double      m_f[11], m_x, m_Q;
    bool        m_calculated[11];

  public:

    CT12_Fortran_Interface(const ATOOLS::Flavour bunch,
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
      std::string cset("");
      std::string path = rpa->gen.Variable("SHERPA_SHARE_PATH")+"/CT12Grid";

      std::string num;
      if (m_member<10) num="0"+ToString(m_member);
      else             num=ToString(m_member);
      std::string asmz[21] = {"0.110", "0.111", "0.112", "0.113", "0.114",
                              "0.115", "0.116", "0.117", "0.118", "0.119",
                              "0.120", "0.121", "0.122", "0.123", "0.124",
                              "0.125", "0.126", "0.127", "0.128", "0.129",
                              "0.130"};
      if (m_set==std::string("ct10nn")) {
        cset = std::string("ct10nn."+num+".pds");
        m_asinfo.m_order=2;
        m_asinfo.m_asmz=0.118;
        m_lhef_number=11200+m_member;
      }
      for (size_t i=0; i<21; ++i) {
        if (m_set==std::string("ct10nn.as"+asmz[i]) && m_member==0) {
          cset = std::string("ct10nn.as"+asmz[i]+".pds");
          m_asinfo.m_order=2;
          m_asinfo.m_asmz=ToType<double>(asmz[i]);
          m_lhef_number=11260+i;
        }
      }
      if (m_set==std::string("ct10n")) {
        cset = std::string("ct10n."+num+".pds");
        m_asinfo.m_order=1;
        m_asinfo.m_asmz=0.118;
        m_lhef_number=11000+m_member;
      }
      for (size_t i=0; i<21; ++i) {
        if (i==0 || i==1 || i==18 || i==19 || i==20) continue;
        if (m_set==std::string("ct10n.as"+asmz[i]) && m_member==0) {
          cset = std::string("ct10n.as"+asmz[i]+".pds");
          m_asinfo.m_order=1;
          m_asinfo.m_asmz=ToType<double>(asmz[i]);
          m_lhef_number=11060+i;
        }
      }
      if (m_set==std::string("ct10nf3") && m_member==0) {
        cset = std::string("ct10nf3.pds");
        m_asinfo.m_order=1;
        m_asinfo.m_asmz=0.1056;
        m_lhef_number=11080;
      }
      if (m_set==std::string("ct10nf32") && m_member==0) {
        cset = std::string("ct10nf32.pds");
        m_asinfo.m_order=1;
        m_asinfo.m_asmz=0.1056;
        m_lhef_number=11081;
        THROW(not_implemented,"An interface for this PDF is not yet provided.");
      }
      if (m_set==std::string("ct10nf4") && m_member==0) {
        cset = std::string("ct10nf4.pds");
        m_asinfo.m_order=1;
        m_asinfo.m_asmz=0.1127;
        m_lhef_number=11082;
      }
      if (m_set==std::string("ct10nf42") && m_member==0) {
        cset = std::string("ct10nf42.pds");
        m_asinfo.m_order=1;
        m_asinfo.m_asmz=0.1127;
        m_lhef_number=11083;
        THROW(not_implemented,"An interface for this PDF is not yet provided.");
      }

      if (m_set==std::string("ct10wn")) {
        cset = std::string("ct10wn."+num+".pds");
        m_asinfo.m_order=1;
        m_asinfo.m_asmz=0.118;
        m_lhef_number=11100+m_member;
      }
      for (size_t i=0; i<21; ++i) {
        if (i==0 || i==1 || i==18 || i==19 || i==20) continue;
        if (m_set==std::string("ct10wn.as"+asmz[i]) && m_member==0) {
          cset = std::string("ct10wn.as"+asmz[i]+".pds");
          m_asinfo.m_order=1;
          m_asinfo.m_asmz=ToType<double>(asmz[i]);
          m_lhef_number=11160+i;
        }
      }
      if (m_set==std::string("ct10wnf3") && m_member==0) {
        cset = std::string("ct10wnf3.pds");
        m_asinfo.m_order=1;
        m_asinfo.m_asmz=0.1056;
        m_lhef_number=11180;
      }
      if (m_set==std::string("ct10wnf32") && m_member==0) {
        cset = std::string("ct10wnf32.pds");
        m_asinfo.m_order=1;
        m_asinfo.m_asmz=0.1056;
        m_lhef_number=11181;
        THROW(not_implemented,"An interface for this PDF is not yet provided.");
      }
      if (m_set==std::string("ct10wnf4") && m_member==0) {
        cset = std::string("ct10wnf4.pds");
        m_asinfo.m_order=1;
        m_asinfo.m_asmz=0.1127;
        m_lhef_number=11182;
      }
      if (m_set==std::string("ct10wnf42") && m_member==0) {
        cset = std::string("ct10wnf42.pds");
        m_asinfo.m_order=1;
        m_asinfo.m_asmz=0.1127;
        m_lhef_number=11183;
        THROW(not_implemented,"An interface for this PDF is not yet provided.");
      }

      if (cset=="") {
        THROW(fatal_error,"PDF set "+m_set
                          +" ("+ToString(m_member)+") not found.");
      }

      char buffer[1024];
      char * err = getcwd(buffer,1024);
      if (err==NULL) {
        msg_Error()<<"Error in CT12_Fortran_Interface.C "<<std::endl;
      }
      int stat=chdir(path.c_str());
      msg_Tracking()<<METHOD<<"(): Init cset "<<cset<<"."<<std::endl;
      char tablefile[40];
      MakeFortranString(tablefile,cset,40);
      setct12_(tablefile);
      if (stat==0) {
        stat=chdir(buffer);
      }
      else {
        msg_Error()<<"Error in CT12_Fortran_Interface.C "<<std::endl
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
      PDF_Base *copy = new CT12_Fortran_Interface(m_bunch,m_set,m_member);
      m_copies.push_back(copy);
      return copy;
    }


    void   CalculateSpec(double x, double Q2) {
      for (size_t i=0;i<11;++i) m_calculated[i]=false;
      m_x=x;
      m_Q=sqrt(Q2);
    }


    double GetXPDF(const ATOOLS::Flavour infl) {
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
        m_f[5-cteqindex]=ct12pdf_(cteqindex,m_x,m_Q)*m_x;
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


DECLARE_PDF_GETTER(CT12_Getter);

PDF_Base *CT12_Getter::operator()
  (const Parameter_Type &args) const
{
  if (!args.m_bunch.IsHadron()) return NULL;
  return new CT12_Fortran_Interface(args.m_bunch,args.m_set,args.m_member);
}

void CT12_Getter::PrintInfo
(std::ostream &str,const size_t width) const
{
  str<<"CT12 fit, see arXiv:1007.2241 and arXiv:1302.6246 [hep-ph]";
}


CT12_Getter *p_get_ct12[63];
extern "C" void InitPDFLib()
{
  p_get_ct12[0]  = new CT12_Getter("ct10nn");
  p_get_ct12[1]  = new CT12_Getter("ct10n");
  p_get_ct12[2]  = new CT12_Getter("ct10nf3");
  p_get_ct12[3]  = new CT12_Getter("ct10nf32");
  p_get_ct12[4]  = new CT12_Getter("ct10nf4");
  p_get_ct12[5]  = new CT12_Getter("ct10nf42");

  p_get_ct12[6]  = new CT12_Getter("ct10wn");
  p_get_ct12[7]  = new CT12_Getter("ct10wnf3");
  p_get_ct12[8]  = new CT12_Getter("ct10wnf32");
  p_get_ct12[9]  = new CT12_Getter("ct10wnf4");
  p_get_ct12[10] = new CT12_Getter("ct10wnf42");

  std::string asmz[21] = {"0.110", "0.111", "0.112", "0.113", "0.114",
                          "0.115", "0.116", "0.117", "0.118", "0.119",
                          "0.120", "0.121", "0.122", "0.123", "0.124",
                          "0.125", "0.126", "0.127", "0.128", "0.129",
                          "0.130"};
  for (size_t i(0);i<21;++i) {
    p_get_ct12[10+i] = new CT12_Getter("ct10nnlo.as"+asmz[i]);
    if (i==0 || i==1 || i==18 || i==19 || i==20) continue;
    p_get_ct12[29+i] = new CT12_Getter("ct10n.as"+asmz[i]);
    p_get_ct12[45+i] = new CT12_Getter("ct10wn.as"+asmz[i]);
  }
}

extern "C" void ExitPDFLib()
{
  for (int i(0);i<63;++i) delete p_get_ct12[i];
}
