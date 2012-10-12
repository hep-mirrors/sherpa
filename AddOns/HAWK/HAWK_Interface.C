#ifndef AddOns_HAWK_HAWK_Interface_H
#define AddOns_HAWK_HAWK_Interface_H

#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Process/ME_Generator_Base.H"

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

  }; // end of class HAWK_Interface
 
} // end of namespace HAWK

#endif

#include "MODEL/Main/Model_Base.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "AddOns/HAWK/HAWK_Wrapper.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"

using namespace HAWK;
using namespace PHASIC; 
using namespace ATOOLS;
using namespace std;

HAWK_Interface::HAWK_Interface(): 
  ME_Generator_Base("HAWK")
{
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
  param_.alphasZ = model->ScalarConstant(string("alpha_S(MZ)"));
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
      param_v[3*i+j]  = model->ComplexMatrixElement(string("CKM"),i,j);
      param_cv[3*i+j] = model->ComplexMatrixElement(string("CKM"),i,j).conj();
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
  for (i(0);i<3;i++) {
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
  qf_.guu[0] = -qu*xsw/xcw; 
  qf_.guu[2] = -qu*xsw/xcw+0.5/(xsw*xcw); 
  qf_.gdd[0] = -qd*xsw/xcw; 
  qf_.gdd[2] = -qd*xsw/xcw-0.5/(xsw*xcw); 
  qf_.guu[1] = qf_.gdd[1] = xparam_.zero;

  for (i(0);i<5;i++) {
    fmass_.rmf2[i] = (double(i)+1.)*1.e-20;
    fmass_.cmf2[i] = Complex((double(i)+1.)*1.e-20,0.0);
  }

  gevfb_.gevfb = rpa.gen.PicoBarn()*1000.;
  
  uv_.mudim   = sqr(param_.mw);
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

Process_Base *HAWK_Interface::InitializeProcess(const Process_Info &pi, bool add)
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
