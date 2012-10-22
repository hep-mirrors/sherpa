#ifndef AddOns_HAWK_HAWK_Interface_H
#define AddOns_HAWK_HAWK_Interface_H

#include "AddOns/HAWK/HAWK_Wrapper.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "MODEL/Main/Model_Base.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"

namespace HAWK {

  class HAWK_Interface: public PHASIC::ME_Generator_Base {
  public :

    // constructor
    HAWK_Interface();

    // destructor
    ~HAWK_Interface();

    // member functions
    bool Initialize(const std::string &path,const std::string &file,
		    MODEL::Model_Base *const model,
		    BEAM::Beam_Spectra_Handler *const beam,
		    PDF::ISR_Handler *const isr);
    PHASIC::Process_Base * InitializeProcess
    (const PHASIC::Process_Info &pi, bool add);
    bool PerformTests();

    void SetClusterDefinitions(PDF::Cluster_Definitions_Base *const defs);

    ATOOLS::Cluster_Amplitude *ClusterConfiguration
    (PHASIC::Process_Base *const proc,const size_t &mode,const double &kt2);

  }; 

  class HAWK_Process: public PHASIC::Virtual_ME2_Base {
  private:
    MODEL::Running_AlphaS * p_as;
    double                * p_p;
    double                * p_m2i, * p_m2i0, * p_m2if, * p_m2if0;
    void CallHAWK(const int & i,const int & j,const int & k,const int & l);
  public:
    HAWK_Process(const PHASIC::Process_Info& pi,
		 const ATOOLS::Flavour_Vector& flavs);
    ~HAWK_Process();
    void Calc(const ATOOLS::Vec4D_Vector& momenta);
    double Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom);
  }; // end of class HAWK_Interface
}

extern "C" { 
  void mat2_(double * p,
	     double * m2i0,double * m2if0,
	     double * m2i,double * m2if);
}

using namespace HAWK;
using namespace PHASIC; 
using namespace ATOOLS;
using namespace std;

HAWK_Interface::HAWK_Interface(): 
  ME_Generator_Base("HAWK")
{
  msg_Out()<<METHOD<<".\n";
}

HAWK_Interface::~HAWK_Interface() 
{
}

bool HAWK_Interface::Initialize
(const string &path,const string &file,MODEL::Model_Base *const model,
 BEAM::Beam_Spectra_Handler *const beam,PDF::ISR_Handler *const isrhandler)
{
  param_.pi      = M_PI;
  param_.alpha   = model->ScalarFunction(string("alpha_QED"),
					 sqr(Flavour(kf_Z).Mass()));
  param_.alpha0  = model->ScalarConstant(string("alpha_QED(0)"));
  param_.alphaZ  = param_.alpha;
  param_.el      = sqrt(4.*param_.pi*param_.alpha);
  param_.GF      = model->ScalarConstant(string("GF"));
  param_.alphas  = model->ScalarConstant(string("alpha_S(MZ)"));
  param_.sw2     = model->ScalarConstant(string("sin2_thetaW"));
  param_.sw      = sqrt(param_.sw2);
  param_.cw2     = 1.-param_.sw2;
  param_.cw      = sqrt(param_.cw2);
  double wfac(sqrt(1.+sqr(Flavour(kf_Wplus).Width()/
				  Flavour(kf_Wplus).Mass())));
  param_.mw      = Flavour(kf_Wplus).Mass()/wfac;
  param_.mw2     = sqr(param_.mw);
  param_.gw      = Flavour(kf_Wplus).Width()/wfac;
  double zfac(sqrt(1.+sqr(Flavour(kf_Z).Width()/
				  Flavour(kf_Z).Mass())));
  param_.mz      = Flavour(kf_Z).Mass()/zfac;
  param_.mz2     = sqr(param_.mz);
  param_.gz      = Flavour(kf_Z).Width()/zfac;
  param_.mh      = Flavour(kf_h0).Mass();
  param_.mh2     = sqr(param_.mh);
  for (int i=0;i<3;i++) {
    for (int j=0;j<3;j++) {
      param_.v[3*i+j]  = model->ComplexMatrixElement(string("CKM"),i,j);
      param_.cv[3*i+j] = conj(model->ComplexMatrixElement(string("CKM"),i,j));
    }
  }
  param_.ml[0]  = Flavour(kf_e).Mass();
  param_.ml[1]  = Flavour(kf_mu).Mass();
  param_.ml[2]  = Flavour(kf_tau).Mass();
  param_.mqp[0] = Flavour(kf_u).HadMass();
  param_.mqp[1] = Flavour(kf_c).HadMass();
  param_.mqp[2] = Flavour(kf_t).Mass();
  param_.mqm[0] = Flavour(kf_d).HadMass();;
  param_.mqm[1] = Flavour(kf_s).HadMass();;
  param_.mqm[2] = Flavour(kf_b).Mass();;
  for (int i(0);i<3;i++) {
    param_.ml2[i]  = sqr(param_.ml[i]);
    param_.mqp2[i] = sqr(param_.mqp[i]);
    param_.mqm2[i] = sqr(param_.mqm[i]);
  }
  rcoptions_.qborn       = 1;
  rcoptions_.qw          = 1;
  rcoptions_.qz          = 1;
  rcoptions_.qschan      = 1;
  rcoptions_.qtchan      = 1;
  rcoptions_.qch2        = 1;
  rcoptions_.qchint      = 1;
  rcoptions_.qbini       = 0;// have to check for b's in 93.
  rcoptions_.qbfin       = 0;// have to check for b's in 93.
  rcoptions_.qwidth      = 1;
  rcoptions_.qfact       = 0;
  rcoptions_.qbos        = 0;
  rcoptions_.qferm       = 0;
  rcoptions_.qsoft       = 1;
  rcoptions_.qhh2        = 1;
  rcoptions_.qqcddiag    = 0;
  rcoptions_.qqcdnondiag = 0;
  rcoptions_.qqcdgsplit  = 0;
  rcoptions_.qqcdggfus   = 0;
  rcoptions_.qcp         = 0;

  xparam_.xmw2 = cparam_.cmw2 = Complex(param_.mw2,-param_.mw*param_.gw);
  xparam_.xmz2 = cparam_.cmz2 = Complex(param_.mz2,-param_.mz*param_.gz);
  xparam_.xcw2 = cparam_.ccw2 = cparam_.cmw2/cparam_.cmz2;
  xparam_.xsw2 = cparam_.csw2 = 1.-cparam_.ccw2;
  xparam_.xmw  = cparam_.cmw  = sqrt(cparam_.cmw2);
  xparam_.xmz  = cparam_.cmz  = sqrt(cparam_.cmz2);
  xparam_.xcw  = cparam_.ccw  = sqrt(cparam_.ccw2);
  xparam_.xsw  = cparam_.csw  = sqrt(cparam_.csw2);

  xparam_.xmh  = param_.mh;
  xparam_.xmh2 = param_.mh2;
  for (int i=0;i<3;i++) {
    xparam_.xml[i]   = param_.ml[i];
    xparam_.xml2[i]  = param_.ml2[i];
    xparam_.xmqp[i]  = param_.mqp[i];
    xparam_.xmqp2[i] = param_.mqp2[i];
    xparam_.xmqm[i]  = param_.mqm[i];
    xparam_.xmqm2[i] = param_.mqm2[i];
  }
  xparam_.zero = Complex(0.,0.);
  qf_.qf[0]  = qf_.qf[2] = qf_.qu = 2./3.;
  qf_.qf[1]  = qf_.qf[3] = qf_.qd = -1./3.;
  qf_.ql     = -1.;
  qf_.qn     = 0.;
  qf_.mu     = param_.mqm[0];
  qf_.mu2    = sqr(qf_.mu);
  qf_.md     = param_.mqp[0];
  qf_.md2    = sqr(qf_.md);
  qf_.guu[0] = -qf_.qu*xparam_.xsw/xparam_.xcw; 
  qf_.guu[2] = -qf_.qu*xparam_.xsw/xparam_.xcw+0.5/(xparam_.xsw*xparam_.xcw); 
  qf_.gdd[0] = -qf_.qd*xparam_.xsw/xparam_.xcw; 
  qf_.gdd[2] = -qf_.qd*xparam_.xsw/xparam_.xcw-0.5/(xparam_.xsw*xparam_.xcw); 
  qf_.guu[1] = qf_.gdd[1] = xparam_.zero;

  for (int i=0;i<5;i++) {
    fmass_.rmf2[i] = (double(i)+1.)*1.e-20;
    fmass_.cmf2[i] = Complex((double(i)+1.)*1.e-20,0.0);
  }

  gevfb_.gevfb = ATOOLS::rpa->Picobarn()*1000.;
  
  uv_.mudim2  = sqr(param_.mw);
  ir_.lambda  = param_.mw; 
  ir_.lambda2 = sqr(ir_.lambda); 
  ir_.deltas  = 0.1; 

  constsub_.constff =  1.0-sqr(M_PI)/3.;
  constsub_.constfi =  1.0-sqr(M_PI)/2.;
  constsub_.constif = -1.5+sqr(M_PI)/6.;
  constsub_.constii =  1.5-sqr(M_PI)/3.;

  hh2_.chh2 = 62.0308*sqr(param_.GF*param_.mh2/sqr(4.*M_PI)/sqrt(2.));
  return true;
}

Process_Base *HAWK_Interface::
InitializeProcess(const Process_Info &pi, bool add)
{
  return NULL;
}

bool HAWK_Interface::PerformTests()
{
  return true;
}
  
void HAWK_Interface::SetClusterDefinitions
(PDF::Cluster_Definitions_Base *const defs)
{
}

Cluster_Amplitude *HAWK_Interface::ClusterConfiguration
(Process_Base *const proc,const size_t &mode,const double &kt2)
{
  return NULL;
}

HAWK_Process::HAWK_Process(const Process_Info& pi,
			   const Flavour_Vector& flavs) :
  Virtual_ME2_Base(pi,flavs), 
  p_as((MODEL::Running_AlphaS *)
       MODEL::s_model->GetScalarFunction(std::string("alpha_S")))
{
  rpa->gen.AddCitation
    (1,"The NLO matrix elements have been taken from HAWK \\cite{}.");

  p_p = new double[24];
  for (size_t i=0;i<24;i++) p_p[i] = 0.;
  p_m2i   = new double[sqr(2*HAWK_NF+1)];
  p_m2i0  = new double[sqr(2*HAWK_NF+1)];
  p_m2if  = new double[sqr(sqr(2*HAWK_NF+1))];
  p_m2if0 = new double[sqr(sqr(2*HAWK_NF+1))];
}

HAWK_Process::~HAWK_Process() { }


void HAWK_Process::CallHAWK(const int & i,const int & j,
			    const int & k,const int & l) {
  mat2_(p_p,p_m2i0,p_m2if0,p_m2i,p_m2if);

  msg_Out()<<"Born & loop level for {"<<i<<" "<<j<<"} --> {"<<k<<" "<<l<<"}: "
	   <<p_m2if0[mr2(i,j,k,l)]<<" "<<p_m2if[mr2(i,j,k,l)]<<".\n";
}

void HAWK_Process::Calc(const Vec4D_Vector &p)
{
  for (size_t n(0);n<p.size();++n) GetMom(p_p,n,p[n]);
  long int i(m_flavs[0]), j(m_flavs[1]);
  long int k(i), l(j);
  msg_Out()<<"----------- "<<METHOD<<" -----------\n"
	   <<"flavs = {"<<i<<", "<<j<<"} -> {"<<k<<", "<<l<<"}"
	   <<" with "<<p.size()<<" vectors.\n";
  for (size_t n=0;n<p.size();n++) {
    msg_Out()<<" p["<<n<<"] = (";
    for (size_t i=0;i<3;i++) msg_Out()<<p_p[mp(n,i)]<<",";
    msg_Out()<<p_p[mp(n,3)]<<") from "<<p[n]<<" "
	     <<"("<<sqrt(dabs(p[n].Abs2()))<<")\n";
  }
  CallHAWK(i,j,k,l);

  m_res.Finite() = p_m2if[mr2(i,j,k,l)];
  m_res.IR()     = p_m2if0[mr2(i,j,k,l)];
  m_res.IR2()    = p_m2if0[mr2(i,j,k,l)];

  msg_Debugging()<<METHOD<<" yields "<<m_res.Finite()
		 <<" + 1/eps * "<<m_res.IR()
		 <<" + 1/eps^2 * "<<m_res.IR2()
		 <<" ...  .\n";
}

double HAWK_Process::Eps_Scheme_Factor(const Vec4D_Vector& mom)
{
  return 4.*M_PI;
}

namespace PHASIC {

  DECLARE_GETTER(HAWK_Interface_Getter,"HAWK",
		 ME_Generator_Base,ME_Generator_Key);

  ME_Generator_Base *HAWK_Interface_Getter::
  operator()(const ME_Generator_Key &key) const
  {
    return new HAWK_Interface();
  }

  void HAWK_Interface_Getter::
  PrintInfo(ostream &str,const size_t width) const
  { 
    str<<"Interface to the HAWK loop ME generator"; 
  }
}

DECLARE_VIRTUALME2_GETTER(HAWK_ME_Getter,"HAWK_ME")
Virtual_ME2_Base * HAWK_ME_Getter::operator()(const Process_Info &pi) const
{
  msg_Out()<<"\n\n"<<"**** In "<<METHOD<<":\n"
	   <<"**** loop generator = "<<pi.m_loopgenerator<<"\n"
	   <<"**** model          = "<<MODEL::s_model->Name()<<"\n";
  if (pi.m_loopgenerator!="HAWK")                       return NULL;
  if (MODEL::s_model->Name()!=std::string("SM") ||
      MODEL::s_model->ScalarConstant("Yukawa_b")>0. ||
      !Flavour(kf_h0).IsOn())                           return NULL;
  if (pi.m_fi.m_nloewtype!=nlo_type::lo)                return NULL;
  if (!(pi.m_fi.m_nloqcdtype&nlo_type::loop))           return NULL;
  Flavour_Vector fl(pi.ExtractFlavours());
  //if (fl[0].IsGluon() || fl[1].IsGluon())               return NULL;
  //if (pi.m_fi.m_ps.size()!=1)                           return NULL;
  //Flavour flh(pi.m_fi.m_ps[0].m_fl[0]);
  //if (!flh==Flavour(kf_h0))                             return NULL;
  msg_Out()<<" size = "<<pi.m_fi.m_ps.size()<<".\n";
  return new HAWK_Process(pi,fl);
}

#endif
