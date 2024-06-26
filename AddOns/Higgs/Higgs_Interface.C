#ifndef AddOns_Higgs_Higgs_Interface_H
#define AddOns_Higgs_Higgs_Interface_H

#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "Higgs_Tree.H"
#include "Higgs_Virtual.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"

namespace HIGGS {

  class Higgs_Interface: public PHASIC::ME_Generator_Base {
  private:

    MODEL::Model_Base *p_model;

    void RegisterDefaults() const;

  public :

    // constructor
    Higgs_Interface();

    // destructor
    ~Higgs_Interface();

    // member functions
    bool Initialize(MODEL::Model_Base *const model,
		    BEAM::Beam_Spectra_Handler *const beam,
		    PDF::ISR_Handler *const isr,
		    YFS::YFS_Handler *const yfs);
    PHASIC::Process_Base *InitializeProcess(const PHASIC::Process_Info &pi, bool add);
    int  PerformTests();
    bool NewLibraries();

  }; // end of class Higgs_Interface

} // end of namespace WHITEHAT

#endif

#include "Wrappers.H"
#include "MODEL/Main/Model_Base.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "ATOOLS/Org/Settings.H"
#include "ATOOLS/Org/Message.H"

using namespace HIGGS;
using namespace PHASIC;
using namespace ATOOLS;

Higgs_Interface::Higgs_Interface():
  ME_Generator_Base("Higgs")
{
  RegisterDefaults();
}

Higgs_Interface::~Higgs_Interface()
{
}

void Higgs_Interface::RegisterDefaults() const
{
  Settings& s = Settings::GetMainSettings();
  s["HIGGS_INTERFERENCE_ONLY"].SetDefault(0);
  s["HIGGS_INTERFERENCE_MODE"].SetDefault(7);
  s["HIGGS_INTERFERENCE_SPIN"].SetDefault(0);
  s["HIGGS_INTERFERENCE_KAPPAG"].SetDefault(1.0);
  s["HIGGS_INTERFERENCE_KAPPAQ"].SetDefault(1.0);
  s["HIGGS_ON_SHELL"].SetDefault(false);
}

bool Higgs_Interface::Initialize(MODEL::Model_Base *const model,
                                 BEAM::Beam_Spectra_Handler *const beam,
                                 PDF::ISR_Handler *const isrhandler,
                                 YFS::YFS_Handler *const yfshandler)
{
  p_model=model;
  Higgs_Tree::SetModel(p_model);
  Higgs_Virtual::SetModel(p_model);
  // Higgs_Virtual::SetModel(p_model);
  s_mt=p_model->GetScalarFunction("m"+Flavour(kf_t).IDName());
  s_mb=p_model->GetScalarFunction("m"+Flavour(kf_b).IDName());
  s_mc=p_model->GetScalarFunction("m"+Flavour(kf_c).IDName());
  G_F=1.0/sqrt(2.0)/std::abs(sqr(model->ComplexConstant("cvev")));
  DEBUG_VAR(model->ComplexConstant("cvev")<<" "<<G_F);
  s2W=std::abs(model->ComplexConstant(std::string("csin2_thetaW")));
  c2W=1.0-s2W;
  sW=sqrt(s2W);
  cW=sqrt(c2W);
  // double ehc_scale2=p_model->ScalarConstant("EHC_SCALE2");
  // DEBUG_VAR(sqrt(ehc_scale2));
  // alpha0=p_model->ScalarFunction("alpha_QED",ehc_scale2);
  // alpha0=p_model->GetInteractionModel()->ScalarConstant("alpha_QED");
  // DEBUG_VAR(1.0/alpha0);
  m_u=Flavour(kf_u).Mass(true);
  m_d=Flavour(kf_d).Mass(true);
  m_s=Flavour(kf_s).Mass(true);
  m_e=Flavour(kf_e).Yuk();
  m_mu=Flavour(kf_mu).Yuk();
  m_tau=Flavour(kf_tau).Yuk();
  m_W=Flavour(kf_Wplus).Mass();
  m_Z=Flavour(kf_Z).Mass();
  sumQsq = ( N_f>=6. ? 3. : 2.) * 4./9. + ( N_f>=5. ? 3. : 2.) * 1./9.;
  DEBUG_VAR(sumQsq<<" "<<1.0/9.0*(3.0+2.0*4.0)<<" "<<N_f);
  sumQ4 = 2. * 16./81. + 3. * 1./81.;
  sumQ6 = 2. * 64./729. + 3. * 1./729.;
  e_u=Flavour(kf_u).Charge();
  e_c=Flavour(kf_c).Charge();
  e_t=Flavour(kf_t).Charge();
  e_d=Flavour(kf_d).Charge();
  e_s=Flavour(kf_s).Charge();
  e_b=Flavour(kf_b).Charge();
  I_3u=Flavour(kf_u).IsoWeak();
  I_3c=Flavour(kf_c).IsoWeak();
  I_3t=Flavour(kf_t).IsoWeak();
  I_3d=Flavour(kf_d).IsoWeak();
  I_3s=Flavour(kf_s).IsoWeak();
  I_3b=Flavour(kf_b).IsoWeak();
  return true;
}

Process_Base *Higgs_Interface::InitializeProcess(const Process_Info &pi, bool add)
{
  return NULL;
}

int Higgs_Interface::PerformTests()
{
  return 1;
}
  
bool Higgs_Interface::NewLibraries()
{
  return false;
}

DECLARE_GETTER(Higgs_Interface,"Higgs",ME_Generator_Base,ME_Generator_Key);

ME_Generator_Base *ATOOLS::Getter<ME_Generator_Base,ME_Generator_Key,
				  Higgs_Interface>::
operator()(const ME_Generator_Key &key) const
{
  msg_Info()<<"#####################################################\n"
	    <<"##                                                 ##\n"
	    <<"##  gg->yy real & virtual corrections computed by  ##\n"
	    <<"##  Z. Bern, L. J. Dixon, C. Schmidt, Y. Li        ##\n"
	    <<"##  Please cite  Phys.Rev. D66 (2002) 074018       ##\n"
	    <<"##               Phys.Rev.Lett. 111 (2013) 111802  ##\n"
	    <<"##                                                 ##\n"
	    <<"#####################################################\n";
  rpa->gen.AddCitation(1,"The Higgs library is described in \\cite{}.");
  return new Higgs_Interface();
}

void ATOOLS::Getter<ME_Generator_Base,ME_Generator_Key,Higgs_Interface>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"Interface to the Higgs loop ME generator"; 
}
