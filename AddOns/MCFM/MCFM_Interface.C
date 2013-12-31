#ifndef AddOns_MCFM_MCFM_Interface_H
#define AddOns_MCFM_MCFM_Interface_H

#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Process/ME_Generator_Base.H"

namespace MCFM {

  class MCFM_Interface: public PHASIC::ME_Generator_Base {
  public :

    // constructor
    MCFM_Interface();

    // destructor
    ~MCFM_Interface();

    // member functions
    bool Initialize(const std::string &path,const std::string &file,
		    MODEL::Model_Base *const model,
		    BEAM::Beam_Spectra_Handler *const beam,
		    PDF::ISR_Handler *const isr);
    PHASIC::Process_Base *InitializeProcess(const PHASIC::Process_Info &pi, bool add);
    int  PerformTests();
    bool NewLibraries();

  }; // end of class MCFM_Interface
 
} // end of namespace MCFM

#endif

#include "MODEL/Main/Model_Base.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "AddOns/MCFM/MCFM_Wrapper.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"

using namespace MCFM;
using namespace PHASIC; 
using namespace ATOOLS;

MCFM_Interface::MCFM_Interface(): 
  ME_Generator_Base("MCFM")
{
}

MCFM_Interface::~MCFM_Interface() 
{
}

bool MCFM_Interface::Initialize
(const std::string &path,const std::string &file,MODEL::Model_Base *const model,
 BEAM::Beam_Spectra_Handler *const beam,PDF::ISR_Handler *const isrhandler)
{
  std::ifstream procfile((rpa->gen.Variable("SHERPA_CPP_PATH")+"/process.DAT").c_str());
  if (!procfile.good())
    THROW(fatal_error,"MCFM's 'process.DAT' is missing. Consider copying it to this directory.");
  msg_Info()<<METHOD<<"(): {\n";
  nproc_.nproc=-1;
  // masses and widths
  nflav_.nflav=Flavour(kf_jet).Size()/2;
  msg_Info()<<"  n_f = "<<nflav_.nflav<<"\n";
  masses_.md=Flavour(kf_d).Mass();
  masses_.mu=Flavour(kf_u).Mass();
  masses_.ms=Flavour(kf_s).Mass();
  masses_.mc=Flavour(kf_c).Mass();
  masses_.mb=Flavour(kf_b).Mass();
  masses_.mt=Flavour(kf_t).Mass();
  masses_.mel=Flavour(kf_e).Mass();
  masses_.mmu=Flavour(kf_mu).Mass();
  masses_.mtau=Flavour(kf_tau).Mass();
  masses_.hmass=Flavour(kf_h0).Mass();
  masses_.hwidth=Flavour(kf_h0).Width();
  masses_.wmass=Flavour(kf_Wplus).Mass();
  masses_.wwidth=Flavour(kf_Wplus).Width();
  masses_.zmass=Flavour(kf_Z).Mass();
  masses_.zwidth=Flavour(kf_Z).Width();
  masses_.twidth=Flavour(kf_t).Width();
  masses_.tauwidth=Flavour(kf_tau).Width();
  masses_.mtausq=sqr(masses_.mtau);
  masses_.mcsq=sqr(Flavour(kf_c).Mass(true));
  masses_.mbsq=sqr(Flavour(kf_b).Mass(true));
  breit_.n2=breit_.n3=0;
  breit_.mass2=Flavour(kf_t).Mass();
  breit_.width2=Flavour(kf_t).Width();
  breit_.mass3=breit_.width3=0.;
  // ew params 
  ewscheme_.ewscheme = 3;
  ewinput_.aemmz_inp = model->ScalarFunction(std::string("alpha_QED"));
  ewinput_.gf_inp    = model->ScalarConstant(std::string("GF"));
  ewinput_.xw_inp    = model->ScalarConstant(std::string("sin2_thetaW"));
  ewinput_.wmass_inp = model->ScalarConstant(std::string("MW"));
  ewinput_.zmass_inp = model->ScalarConstant(std::string("MZ"));
  // ckm elements
  // must check the syntax for off-diag elements.
  cabib_.Vud=model->ComplexMatrixElement(std::string("CKM"),0,0).real();
  cabib_.Vus=model->ComplexMatrixElement(std::string("CKM"),0,1).real();
  cabib_.Vub=model->ComplexMatrixElement(std::string("CKM"),0,2).real();
  cabib_.Vcd=model->ComplexMatrixElement(std::string("CKM"),1,0).real();
  cabib_.Vcs=model->ComplexMatrixElement(std::string("CKM"),1,1).real();
  cabib_.Vcb=model->ComplexMatrixElement(std::string("CKM"),1,2).real();
  //msg_Out()<<"Check this:"
  //	   <<cabib_.Vud<<" "<<cabib_.Vus<<" "<<cabib_.Vub<<std::endl
  //	   <<"           "
  //	   <<cabib_.Vcd<<" "<<cabib_.Vcs<<" "<<cabib_.Vcb<<std::endl;
  // set couplings
  scale_.scale       = ewinput_.zmass_inp;
  scale_.musq        = sqr(scale_.scale);
  nlooprun_.nlooprun = MODEL::as->Order()+1;
  couple_.amz        = model->ScalarFunction(std::string("alpha_S"));

  if (model->Name()==std::string("SM+AGC")){
    anomcoup_.delg1_z  = model->ScalarConstant(std::string("g1_Z"))-1;
    anomcoup_.delg1_g  = model->ScalarConstant(std::string("g1_gamma"))-1;
    anomcoup_.lambda_g = model->ScalarConstant(std::string("lambda_gamma"));
    anomcoup_.lambda_z = model->ScalarConstant(std::string("lambda_Z"));
    anomcoup_.h1Z      = model->ScalarConstant(std::string("h1_Z"));
    anomcoup_.h2z      = model->ScalarConstant(std::string("h2_Z"));
    anomcoup_.h3z      = model->ScalarConstant(std::string("h3_Z"));
    anomcoup_.h4z      = model->ScalarConstant(std::string("h4_Z"));
    anomcoup_.h1gam    = model->ScalarConstant(std::string("h1_gamma"));
    anomcoup_.h2gam    = model->ScalarConstant(std::string("h2_gamma"));
    anomcoup_.h3gam    = model->ScalarConstant(std::string("h3_gamma"));
    anomcoup_.h4gam    = model->ScalarConstant(std::string("h4_gamma"));
    anomcoup_.delk_g   = model->ScalarConstant(std::string("kappa_gamma"))-1;
    anomcoup_.delk_z   = model->ScalarConstant(std::string("kappa_Z"))-1;
    anomcoup_.tevscale = model->ScalarConstant(std::string("UNITARIZATION_SCALE"))/1000;
    }
  
  if (!zerowidth_.zerowidth) limits_.bbsqmin = 1.; 

  qcdcouple_.as      = model->ScalarFunction(std::string("alpha_S"));
  qcdcouple_.gsq     = 4.*M_PI*qcdcouple_.as;
  qcdcouple_.ason2pi = qcdcouple_.as/(2.*M_PI);
  qcdcouple_.ason4pi = qcdcouple_.as/(4.*M_PI);
  std::string dummy  = std::string("mstw8lo");
  dummy.copy(pdlabel_.pdlabel,7);
  limits_.wsqmin = 1.e-6;
  limits_.wsqmax = 1.e99;

  msg_Info()<<"}\n";

  verbose_.verbose = true;
  return true;
}

Process_Base *MCFM_Interface::InitializeProcess(const Process_Info &pi, bool add)
{
  return NULL;
}

int MCFM_Interface::PerformTests()
{
  return 1;
}

bool MCFM_Interface::NewLibraries()
{
  return false;
}

DECLARE_GETTER(MCFM_Interface,"MCFM",ME_Generator_Base,ME_Generator_Key);

ME_Generator_Base *ATOOLS::Getter
<ME_Generator_Base,ME_Generator_Key,MCFM_Interface>::
operator()(const ME_Generator_Key &key) const
{
  return new MCFM_Interface();
}

void ATOOLS::Getter<ME_Generator_Base,ME_Generator_Key,MCFM_Interface>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"Interface to the MCFM loop ME generator"; 
}
