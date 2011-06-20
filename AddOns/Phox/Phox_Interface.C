#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "MODEL/Main/Running_AlphaQED.H"

namespace PHOX {

  class Phox: public PHASIC::Virtual_ME2_Base {
  private:
    MODEL::Running_AlphaS   * p_as;
    MODEL::Running_AlphaQED * p_aqed;
    double                    m_normcorr, m_cplfac, m_qcharge4;

    bool DiPhoxME2(const ATOOLS::Vec4D_Vector& mom,const double & q2);
  public:
    Phox(const PHASIC::Process_Info& pi,
	 const ATOOLS::Flavour_Vector& flavs);
    ~Phox();
    void   Calc(const ATOOLS::Vec4D_Vector& momenta);
    double Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& momenta);
  };
}


#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace PHOX;
using namespace PHASIC;
using namespace ATOOLS;

Phox::Phox(const int & pID,const Process_Info& pi,
	   const Flavour_Vector& flavs):
  Virtual_ME2_Base(pi,flavs), 
  p_as((MODEL::Running_AlphaS *)
       MODEL::s_model->GetScalarFunction(std::string("alpha_S"))),
  p_aqed((MODEL::Running_AlphaQED *)
       MODEL::s_model->GetScalarFunction(std::string("alpha_QED"))),
  m_normcorr(4.*(3.0-1.0)*(3.0+1.0)/(3.0*3.0)),
  m_qcharge4(pow(flavs[0].Charge(),4))
{
  rpa.gen.AddCitation
    (1,"The NLO matrix elements have been taken from PHOX \\cite{}.");
}

Phox::~Phox()
{
}

bool Phox::DiPhoxME2(const Vec4D_Vector &mom,const double & q2) {
  double s = 2.*mom[0]*mom[1];
  double t = 2.*mom[0]*mom[2];
  double u = 2.*mom[1]*mom[2];
  double prefactor = m_cplfac/(t*u);

  m_res.finite() = M_PI*M_PI*t*t + M_PI*M_PI*u*u + log(s/q2)*
    ((-2.*s*s- 2.*t*t)*log(-(t/q2)) + (-2.*s*s - 2.*u*u)*log(-(u/q2)) - 
     2.*s*s - 2.*t*t - 4.*t*u - 2.*u*u) + (s*s + t*t)*pow(log(-(t/q2)),2) + 
    s*s*pow(log(-(u/q2)),2) + 2.*s*s*pow(log(s/q2),2) + 
    3.*t*t*log(-(u/q2)) + u*(2.*t + 3.*u)*log(-(t/q2)) + 
    2.*t*u*log(-(u/q2)) + u*u*pow(log(-(u/q2)),2) + s*s - 4.*t*t - 4.*u*u;

  m_res.IR1() = 2.0*s*s-u*u-t*t+2.0*log(s/q2)*u*u+2.0*log(s/q2)*t*t;
  m_res.IR2() = -2.*(u*u+t*t);
  
  m_res.IR1()    *= prefactor;
  m_res.IR2()    *= prefactor;
  m_res.finite() *= prefactor;
}

void Phox::Calc(const Vec4D_Vector &p)
{
  m_cplfac  = (*p_as)(m_mur2)*(*p_aqed)(m_mur2)*m_qcharge4;
  // can check ME by checking independence of q2.
  double q2 = 1.; 
  if (!DiPhoxME2(p,q2)) {
    msg_Error()<<"Error in "<<METHOD<<":\n"
	       <<"   Could not evaluate virtual ME^2 of qqb -> gamma gamma.\n";
  }
}

double Phox::Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom)
{
  return 4.*M_PI;
}

extern "C" { void chooser_(); }

DECLARE_VIRTUALME2_GETTER(Phox_Getter,"Phox")
Virtual_ME2_Base *Phox_Getter::operator()(const Process_Info &pi) const
{
  DEBUG_FUNC("");
  if (pi.m_loopgenerator!="PHOX")        return NULL;
  if (pi.m_oew!=2)                       return NULL;
  if (pi.m_fi.m_nloewtype!=nlo_type::lo) return NULL;
  if (pi.m_fi.m_nloqcdtype&nlo_type::loop) {
    Flavour_Vector fl(pi.ExtractFlavours());
    if (fl.size()!=4 ||
	!(fl[0].IsQuark() && fl[1]==fl[0].Bar() &&
	  fl[2].IsPhoton() && fl[3].IsPhoton())) {
      msg_Error()<<"Error in "<<METHOD<<":\n"
		 <<"   Tried to initialse a phox interface with flavours";
      for (int i=0;i<fl.size();i++) msg_Error()<<" "<<fl[i];
      msg_Error()<<".\n"<<"   Return 'NULL' and hope for the best.\n";
      return NULL;
    }
    if (!(MODEL::s_model->Name()==std::string("SM") ||
	  MODEL::s_model->Name()==std::string("THDM"))) {
      msg_Error()<<"Error in "<<METHOD<<":\n"
		 <<"   Tried to initialse a phox interface with model "
		 <<MODEL::s_model->Name()<<".\n   Return 'NULL' and hope.\n";
      return NULL;
    }
    
    return new Phox(pi,fl);
  }
  return NULL;
}
