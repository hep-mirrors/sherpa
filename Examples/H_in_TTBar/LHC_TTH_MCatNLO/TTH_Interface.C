#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "MODEL/Main/Model_Base.H"
#include "MODEL/Interaction_Models/Interaction_Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"

#define TTH_NMX 4

inline int mp(const int id,const int i)
{ return i+id*TTH_NMX; }

inline void GetMom(double *m,const int n,const ATOOLS::Vec4D &p)
{ m[mp(n,0)]=p[0]; for (int i(1);i<4;++i) m[mp(n,i)]=p[i]; }

extern "C" {

inline void MakeFortranString
(char *output,std::string input,unsigned int length)
{
  for (unsigned int i=0;i<length;++i) 
    output[i]=(char)32;
  for (size_t j=0;j<input.length();++j) 
    output[j]=(char)input[j];
}

  extern struct{
    double mtop,alphas,mb;
  } topas_;
  extern struct{
    double mh;
  } smhiggs_;
  extern struct{
    double delta,muedr;
  } dimreg_;
  extern struct {
    double yukt,yukb,yukf;
  } yukawa_;

  void m2_virt_tth_(double *p,double *res,int *mode,int *type);

}

using namespace PHASIC;
using namespace ATOOLS;

namespace TTH {

  class TTH_Interface: public PHASIC::ME_Generator_Base {
  public :

    // constructor
    TTH_Interface(): ME_Generator_Base("TTH") {}

    // member functions
    bool Initialize(const std::string &path,const std::string &file,
		    MODEL::Model_Base *const model,
		    BEAM::Beam_Spectra_Handler *const beam,
		    PDF::ISR_Handler *const isr)
    {
      msg_Info()<<"#############################################################\n"
		<<"##                                                         ##\n"
		<<"##  Top anti-top Higgs  virtual corrections  computed  by  ##\n"
		<<"##  S. Dawson, C. Jackson, L. Orr, L. Reina, D. Wackeroth  ##\n"
		<<"##  Please cite  Phys. Rev. Lett. 87 (2001) 201804         ##\n"
		<<"##               Phys. Rev. D65 (2002) 053017              ##\n"
		<<"##               Phys. Rev. D67 (2003) 071503              ##\n"
		<<"##               Phys. Rev. D68 (2003) 034022              ##\n"
		<<"##                                                         ##\n"
		<<"#############################################################\n";
      smhiggs_.mh=Flavour(kf_h0).Mass();
      topas_.mb=Flavour(kf_b).Mass();
      topas_.mtop=Flavour(kf_t).Mass();
      topas_.alphas=model->ScalarFunction
	(std::string("alpha_S"),rpa->gen.CplScale());
      double gf(1.0/sqrt(2.0)/std::abs(sqr(model->ComplexConstant("cvev"))));
      yukawa_.yukt=-model->GetInteractionModel()->
	ScalarFunction("m"+Flavour(kf_t).IDName(),sqr(smhiggs_.mh))*sqrt(gf*sqrt(2.0));
      yukawa_.yukb=-model->GetInteractionModel()->
	ScalarFunction("m"+Flavour(kf_b).IDName(),sqr(smhiggs_.mh))*sqrt(gf*sqrt(2.0));
      return true;
    }

    PHASIC::Process_Base *InitializeProcess
    (const PHASIC::Process_Info &pi, bool add) { return NULL; }
    int  PerformTests() { return 1; }
    bool NewLibraries() { return false; }

    void SetClusterDefinitions(PDF::Cluster_Definitions_Base *const defs) {}

    ATOOLS::Cluster_Amplitude *ClusterConfiguration
    (PHASIC::Process_Base *const proc,const size_t &mode,const double &kt2)
    { return NULL; }

  }; // end of class TTH_Interface
 
  class TTH_Process: public PHASIC::Virtual_ME2_Base {
  protected:
    double *p_p, *p_res;
  public:
    TTH_Process(const PHASIC::Process_Info& pi,
		const ATOOLS::Flavour_Vector& flavs):
      Virtual_ME2_Base(pi,flavs)
    {
      p_p = new double[4*(TTH_NMX+1)];
      p_res = new double[4];
      rpa->gen.AddCitation
	(1,std::string("NLO ME for tth from")
	 +"\\cite{Reina:2001sf}, \\cite{Reina:2001bc}"
	 +", \\cite{Dawson:2002tg}, \\cite{Dawson:2003zu}.");
    }

    ~TTH_Process()
    {
      delete [] p_p;
      delete [] p_res;
    }

    void Calc(const ATOOLS::Vec4D_Vector &p)
    {
      dimreg_.muedr=sqrt(m_mur2);
      topas_.alphas=(*MODEL::as)(m_mur2);
      int s(m_flavs[0].IsAnti()?2:1), m(m_flavs[0].IsQuark()?1:2);
      for (size_t n(0);n<p.size();++n) GetMom(p_p,n,p[n]);
      m2_virt_tth_(p_p,p_res,&m,&s);
      double norm((2.0*M_PI)/topas_.alphas/p_res[3]);
      m_res.Finite()=(p_res[0]+sqr(M_PI)/6.0*p_res[2])*norm;
      m_res.IR()=p_res[1]*norm;
      m_res.IR2()=p_res[2]*norm;
      m_born=p_res[3];
    }

    double Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom)
    {
      return 4.0*M_PI;
    }

  };// end of class TTH_Process

}// end of namespace TTH

using namespace TTH;

DECLARE_GETTER(TTH_Interface,"TTH",ME_Generator_Base,ME_Generator_Key);
ME_Generator_Base *ATOOLS::Getter
<ME_Generator_Base,ME_Generator_Key,TTH_Interface>::
operator()(const ME_Generator_Key &key) const
{
  return new TTH_Interface();
}
void ATOOLS::Getter<ME_Generator_Base,ME_Generator_Key,TTH_Interface>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"Interface to the TTH loop ME generator"; 
}

DECLARE_VIRTUALME2_GETTER(TTH_Process,"TTH_Process")
Virtual_ME2_Base *ATOOLS::Getter
<Virtual_ME2_Base,Process_Info,TTH_Process>::
operator()(const Process_Info &pi) const
{
  if (pi.m_loopgenerator!="TTH") return NULL;
  if (MODEL::s_model->Name()!=std::string("SM")) return NULL;
  if (pi.m_oew!=1 || pi.m_fi.m_nloewtype!=nlo_type::lo)return NULL;
  if (!(pi.m_fi.m_nloqcdtype&nlo_type::loop)) return NULL;
  Flavour_Vector fl(pi.ExtractFlavours());
  if (fl[0].Strong() && fl[1].Strong() &&
      fl[2].Kfcode()==kf_h0 &&
      fl[3].Kfcode()==kf_t && fl[4]==fl[3].Bar() &&
      Flavour(kf_t).Mass()) {
    msg_Info()<<"!";
    return new TTH_Process(pi,fl);
    }
  return NULL;
}
