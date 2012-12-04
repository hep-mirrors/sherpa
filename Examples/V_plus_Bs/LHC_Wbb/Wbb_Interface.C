#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"

#define WBB_NMX 4

inline int mp(const int id,const int i)
{ return i+id*WBB_NMX; }

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
    double mb,mt,mw;
  } masses_;
  extern struct{
    double gw;
  } couplings_;
  extern struct{
    double xmur;
  } scale_;

  void m2_virt_wbb_(double *p,double *res);

}

using namespace PHASIC;
using namespace ATOOLS;

namespace WBB {

  class Wbb_Interface: public PHASIC::ME_Generator_Base {
  public :

    // constructor
    Wbb_Interface(): ME_Generator_Base("Wbb") {}

    // member functions
    bool Initialize(const std::string &path,const std::string &file,
		    MODEL::Model_Base *const model,
		    BEAM::Beam_Spectra_Handler *const beam,
		    PDF::ISR_Handler *const isr)
    {
      msg_Info()<<METHOD<<"(): {\n"
		<<"#################################################\n"
		<<"##                                             ##\n"
		<<"##  Wbb~  virtual  corrections   computed  by  ##\n"
		<<"##  F. Febres Cordero, L. Reina, D. Wackeroth  ##\n"
		<<"##  Please cite  Phys. Rev. D74 (2006) 034007  ##\n"
		<<"##                                             ##\n"
		<<"#################################################\n";
      masses_.mb=Flavour(kf_b).Mass();
      masses_.mt=Flavour(kf_t).Mass();
      masses_.mw=Flavour(kf_Wplus).Mass();
      double gf(1.0/sqrt(2.0)/std::abs(sqr(model->ComplexConstant("cvev"))));
      couplings_.gw=sqrt(8*sqr(Flavour(kf_Wplus).Mass())*gf/sqrt(2));
      return true;
    }

    PHASIC::Process_Base *InitializeProcess
    (const PHASIC::Process_Info &pi, bool add) { return NULL; }
    bool PerformTests() { return true; }
    bool NewLibraries() { return false; }

    void SetClusterDefinitions(PDF::Cluster_Definitions_Base *const defs) {}

    ATOOLS::Cluster_Amplitude *ClusterConfiguration
    (PHASIC::Process_Base *const proc,const size_t &mode,const double &kt2)
    { return NULL; }

  }; // end of class Wbb_Interface
 
  DECLARE_GETTER(Wbb_Interface_Getter,"Wbb",ME_Generator_Base,ME_Generator_Key);
  ME_Generator_Base *Wbb_Interface_Getter::operator()(const ME_Generator_Key &key) const
  {
    return new Wbb_Interface();
  }
  void Wbb_Interface_Getter::PrintInfo(std::ostream &str,const size_t width) const
  { 
    str<<"Interface to the Wbb loop ME generator"; 
  }

  class Wbb_Process: public PHASIC::Virtual_ME2_Base {
  protected:
    double *p_p, *p_res;
  public:
    Wbb_Process(const PHASIC::Process_Info& pi,
		const ATOOLS::Flavour_Vector& flavs):
      Virtual_ME2_Base(pi,flavs)
    {
      p_p = new double[4*(WBB_NMX+1)];
      p_res = new double[4];
      // m_mode=1;
      rpa->gen.AddCitation
	(1,"NLO ME for Wbb from \\cite{FebresCordero:2006sj}.");
    }

    ~Wbb_Process()
    {
      delete [] p_p;
      delete [] p_res;
    }

    void Calc(const ATOOLS::Vec4D_Vector &p)
    {
      scale_.xmur=sqrt(m_mur2);
      int s(m_flavs[0].IsAnti());
      for (size_t n(0);n<p.size();++n) GetMom(p_p,n,p[(s&&n<2)?1-n:n]);
      m2_virt_wbb_(p_p,p_res);
      m_res.Finite()=p_res[0];
      m_res.IR()=p_res[1];
      m_res.IR2()=p_res[2];
      m_born=p_res[3]*sqr(4.0*M_PI*(*MODEL::as)(m_mur2));
    }

    double Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom)
    {
      return 4.0*M_PI;
    }

  };// end of class Wbb_Process

  DECLARE_VIRTUALME2_GETTER(Wbb_Process_Getter,"Wbb_Process")
  Virtual_ME2_Base *Wbb_Process_Getter::operator()(const Process_Info &pi) const
  {
    if (pi.m_loopgenerator!="Wbb") return NULL;
    if (MODEL::s_model->Name()!=std::string("SM")) return NULL;
    if (pi.m_oew!=1 || pi.m_fi.m_nloewtype!=nlo_type::lo)return NULL;
    if (!(pi.m_fi.m_nloqcdtype&nlo_type::loop)) return NULL;
    Flavour_Vector fl(pi.ExtractFlavours());
    DEBUG_VAR(fl);
    if (fl[0].Strong() && fl[1].Strong() &&
	fl[2].Kfcode()==kf_Wplus &&
	fl[3].Kfcode()==kf_b && fl[4]==fl[3].Bar() &&
	Flavour(kf_b).Mass()) {
      msg_Info()<<"!";
      return new Wbb_Process(pi,fl);
    }
    return NULL;
  }

}// end of namespace WBB
