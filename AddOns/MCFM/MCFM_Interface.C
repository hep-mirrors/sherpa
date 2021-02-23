#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include "MODEL/Main/Model_Base.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "AddOns/MCFM/MCFM_Wrapper.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"

using namespace PHASIC; 
using namespace ATOOLS;

namespace MCFM {

  class MCFM_Interface: public PHASIC::ME_Generator_Base {
  public :

    MCFM_Interface(): ME_Generator_Base("MCFM") {}

    bool Initialize(const std::string &path,const std::string &file,
		    MODEL::Model_Base *const model,
		    BEAM::Beam_Spectra_Handler *const beam,
		    PDF::ISR_Handler *const isr)
    {
      DEBUG_FUNC("");
      std::string pdname(rpa->gen.Variable("SHERPA_CPP_PATH")+"/process.DAT");
      if (!FileExists(pdname))
	Copy(MCFM_PATH+std::string("Bin/process.DAT"),pdname);
      nproc_.nproc=-1;
      // masses and widths
      nflav_.nflav=Flavour(kf_jet).Size()/2;
      msg_Debugging()<<"n_f = "<<nflav_.nflav<<"\n";
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
      ewscheme_.ewscheme=0;
      ewinput_.aemmz_inp=model->ScalarConstant("alpha_QED");
      ewinput_.gf_inp=1.0/sqrt(2.0)/std::abs(sqr(model->ComplexConstant("cvev")));
      ewinput_.xw_inp=std::abs(model->ComplexConstant("csin2_thetaW"));
      ewinput_.wmass_inp=Flavour(kf_Wplus).Mass();
      ewinput_.zmass_inp=Flavour(kf_Z).Mass();
      // ckm elements
      cabib_.Vud=model->ComplexConstant("CKM_0_0").real();
      cabib_.Vus=model->ComplexConstant("CKM_0_1").real();
      cabib_.Vub=model->ComplexConstant("CKM_0_2").real();
      cabib_.Vcd=model->ComplexConstant("CKM_1_0").real();
      cabib_.Vcs=model->ComplexConstant("CKM_1_1").real();
      cabib_.Vcb=model->ComplexConstant("CKM_1_2").real();
      // scales and strong coupling
      mcfmscale_.scale=ewinput_.zmass_inp;
      mcfmscale_.musq=sqr(mcfmscale_.scale);
      nlooprun_.nlooprun=MODEL::as->Order()+1;
      couple_.amz=model->ScalarConstant("alpha_S");
      if (!zerowidth_.zerowidth) limits_.bbsqmin = 1.;
      qcdcouple_.as=model->ScalarConstant("alpha_S");
      qcdcouple_.gsq=4.*M_PI*qcdcouple_.as;
      qcdcouple_.ason2pi=qcdcouple_.as/(2.*M_PI);
      qcdcouple_.ason4pi=qcdcouple_.as/(4.*M_PI);
      std::string dummy("mstw8lo");
      dummy.copy(pdlabel_.pdlabel,255);
      limits_.wsqmin=1.e-6;
      limits_.wsqmax=1.e99;
      verbose_.verbose=true;
      return true;
    }

    PHASIC::Process_Base *InitializeProcess
    (const PHASIC::Process_Info &pi, bool add) { return NULL; }

    int  PerformTests() { return 1; }
    bool NewLibraries() { return false; }

  }; // end of class MCFM_Interface

} // end of namespace MCFM

using namespace MCFM;

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
