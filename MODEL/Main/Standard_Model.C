#include "MODEL/Main/Standard_Model.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "MODEL/Main/Strong_Coupling.H"
#include "MODEL/Main/Running_Fermion_Mass.H"
#include "PDF/Main/ISR_Handler.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include <iomanip>

using namespace MODEL;
using namespace ATOOLS;
using namespace std;

DECLARE_GETTER(Standard_Model_Getter,
	       "SM",Model_Base,Model_Arguments);

Model_Base *Standard_Model_Getter::
operator()(const Model_Arguments &args) const
{
  return new Standard_Model(args.m_path,args.m_file,args.m_elementary);
}

void Standard_Model_Getter::
PrintInfo(ostream &str,const size_t width) const
{ 
  str<<"The Standard Model\n";
  str<<setw(width+4)<<" "<<"{\n"
     <<setw(width+7)<<" "<<"parameter specification [keyword=value]\n"
     <<setw(width+7)<<" "<<"- EW_SCHEME (values 0,1,2,3, EW input schemes, see documentation)\n"
     <<setw(width+7)<<" "<<"- WIDTH_SCHEME (Fixed or CMS, see documentation)\n"
     <<setw(width+7)<<" "<<"- ALPHAS(MZ) (strong coupling at MZ)\n"
     <<setw(width+7)<<" "<<"- ALPHAS(default) (fixed strong coupling)\n"
     <<setw(width+7)<<" "<<"- ORDER_ALPHAS (0,1,2 -> 1, 2, 3-loop running)\n"
     <<setw(width+7)<<" "<<"- 1/ALPHAQED(0) (alpha QED Thompson limit)\n"
     <<setw(width+7)<<" "<<"- 1/ALPHAQED(default) (fixed alpha QED)\n"
     <<setw(width+7)<<" "<<"- SIN2THETAW (weak mixing angle)\n"
     <<setw(width+7)<<" "<<"- VEV (Higgs vev)\n"
     <<setw(width+7)<<" "<<"- LAMBDA (Higgs quartic coupling)\n"
     <<setw(width+7)<<" "<<"- CKMORDER (0,1,2,3 - order of CKM expansion in Cabibbo angle)\n"
     <<setw(width+7)<<" "<<"- CABIBBO (Cabibbo angle in Wolfenstein parameterization)\n"
     <<setw(width+7)<<" "<<"- A (Wolfenstein A)\n"
     <<setw(width+7)<<" "<<"- RHO (Wolfenstein Rho)\n"
     <<setw(width+7)<<" "<<"- ETA (Wolfenstein Eta)\n"
     <<setw(width+4)<<" "<<"}";
  str<<"Infrared continuation of alphaS:\n";
  str<<setw(width+4)<<" "<<"{\n"
     <<setw(width+7)<<" "<<"- AS_FORM (values 0,1,2,3,10, see documentation)\n"
     <<setw(width+7)<<" "<<"- Q2_AS (corresponding infrared parameter, see documentation)\n"
     <<setw(width+4)<<" "<<"}";
}

DECLARE_GETTER(Standard_Model_Top_Getter,
	       "SM+Top",Model_Base,Model_Arguments);

Model_Base *Standard_Model_Top_Getter::
operator()(const Model_Arguments &args) const
{
  return new Standard_Model(args.m_path,args.m_file,args.m_elementary,1);
}

void Standard_Model_Top_Getter::
PrintInfo(ostream &str,const size_t width) const
{ 
  str<<"The Standard Model + non-standard top couplings\n";
  str<<setw(width+4)<<" "<<"{\n"
     <<setw(width+7)<<" "<<"parameter specification [keyword=value]\n"
     <<setw(width+7)<<" "<<"- EW_SCHEME (values 0,1,2,3, EW input schemes, see documentation)\n"
     <<setw(width+7)<<" "<<"- WIDTH_SCHEME (Fixed or CMS, see documentation)\n"
     <<setw(width+7)<<" "<<"- ALPHAS(MZ) (strong coupling at MZ)\n"
     <<setw(width+7)<<" "<<"- ALPHAS(default) (fixed strong coupling)\n"
     <<setw(width+7)<<" "<<"- ORDER_ALPHAS (0,1,2 -> 1, 2, 3-loop running)\n"
     <<setw(width+7)<<" "<<"- 1/ALPHAQED(0) (alpha QED Thompson limit)\n"
     <<setw(width+7)<<" "<<"- 1/ALPHAQED(default) (fixed alpha QED)\n"
     <<setw(width+7)<<" "<<"- SIN2THETAW (weak mixing angle)\n"
     <<setw(width+7)<<" "<<"- VEV (Higgs vev)\n"
     <<setw(width+7)<<" "<<"- LAMBDA (Higgs quartic coupling)\n"
     <<setw(width+7)<<" "<<"- CKMORDER (0,1,2,3 - order of CKM expansion in Cabibbo angle)\n"
     <<setw(width+7)<<" "<<"- CABIBBO (Cabibbo angle in Wolfenstein parameterization)\n"
     <<setw(width+7)<<" "<<"- A (Wolfenstein A)\n"
     <<setw(width+7)<<" "<<"- RHO (Wolfenstein Rho)\n"
     <<setw(width+7)<<" "<<"- ETA (Wolfenstein Eta)\n"
     <<setw(width+7)<<" "<<"- KAPPA_{TBW} (rel. factor of tbW w.r.t. SM)\n"
     <<setw(width+7)<<" "<<"- THETA_{TBW} (angle of tbW in PL-PR w.r.t. SM)\n"
     <<setw(width+4)<<" "<<"}";
}


DECLARE_GETTER(Standard_Model_Gen4_Getter,"SM+4thGen",
	       Model_Base,Model_Arguments);

Model_Base *Standard_Model_Gen4_Getter::
operator()(const Model_Arguments &args) const
{
  return new Standard_Model(args.m_path,args.m_file,args.m_elementary,2);
}

void Standard_Model_Gen4_Getter::
PrintInfo(ostream &str,const size_t width) const
{ 
  str<<"The Standard Model + 4th generation\n";
  str<<setw(width+4)<<" "<<"{\n"
     <<setw(width+7)<<" "<<"parameter specification [keyword=value]\n"
     <<setw(width+7)<<" "<<"- EW_SCHEME (values 0,1,2,3, EW input schemes, see documentation)\n"
     <<setw(width+7)<<" "<<"- WIDTH_SCHEME (Fixed or CMS, see documentation)\n"
     <<setw(width+7)<<" "<<"- ALPHAS(MZ)          (strong coupling at MZ)\n"
     <<setw(width+7)<<" "<<"- ALPHAS(default)     (fixed strong coupling)\n"
     <<setw(width+7)<<" "<<"- ORDER_ALPHAS        (0,1,2 -> 1, 2, 3-loop running)\n"
     <<setw(width+7)<<" "<<"- 1/ALPHAQED(0)       (alpha QED Thompson limit)\n"
     <<setw(width+7)<<" "<<"- 1/ALPHAQED(default) (fixed alpha QED)\n"
     <<setw(width+7)<<" "<<"- SIN2THETAW          (weak mixing angle)\n"
     <<setw(width+7)<<" "<<"- VEV    (Higgs vev)\n"
     <<setw(width+7)<<" "<<"- LAMBDA (Higgs quartic coupling)\n"
     <<setw(width+7)<<" "<<"masses, widths, etc:\n"
     <<setw(width+7)<<" "<<"- MASS[7,8,17,18], WIDTH[7,8,17,18], etc.\n"
     <<setw(width+7)<<" "<<"3x3 CKM matrix:\n"
     <<setw(width+7)<<" "<<"- CKMORDER (0,1,2,3 - order of CKM expansion in Cabibbo angle)\n"
     <<setw(width+7)<<" "<<"- CABIBBO  (Cabibbo angle in Wolfenstein parameterization)\n"
     <<setw(width+7)<<" "<<"- A        (Wolfenstein A)\n"
     <<setw(width+7)<<" "<<"- RHO      (Wolfenstein Rho)\n"
     <<setw(width+7)<<" "<<"- ETA      (Wolfenstein Eta)\n"
     <<setw(width+7)<<" "<<"4th generation extension (Phys.Lett.B192:441,1987):\n"
     <<setw(width+7)<<" "<<"- A_14     (quark mixing angle a_14)\n"
     <<setw(width+7)<<" "<<"- A_24     (quark mixing angle a_24)\n"
     <<setw(width+7)<<" "<<"- A_34     (quark mixing angle a_34)\n"
     <<setw(width+7)<<" "<<"- PHI_2    (quark mixing angle phi_2)\n"
     <<setw(width+7)<<" "<<"- PHI_3    (quark mixing angle phi_3)\n"
     <<setw(width+7)<<" "<<"4th generation lepton mixing:\n"
     <<setw(width+7)<<" "<<"- THETA_L14   (lepton mixing angle theta_L,14)\n"
     <<setw(width+7)<<" "<<"- THETA_L24   (lepton mixing angle theta_L,24)\n"
     <<setw(width+7)<<" "<<"- THETA_L34   (lepton mixing angle theta_L,34)\n"
     <<setw(width+7)<<" "<<"- PHI_L2    (lepton mixing angle phi_L,2)\n"
     <<setw(width+7)<<" "<<"- PHI_L3    (lepton mixing angle phi_L,3)\n"
     <<setw(width+7)<<" "<<"possibly different Yukawa couplings for 4th generation:\n"
     <<setw(width+7)<<" "<<"- YUKAWA_U4   (Yukawa coupling of 4th generation up-type quark)\n"
     <<setw(width+7)<<" "<<"- YUKAWA_D4   (Yukawa coupling of 4th generation down-type quark)\n"
     <<setw(width+7)<<" "<<"- YUKAWA_L4   (Yukawa coupling of 4th generation lepton)\n"
     <<setw(width+7)<<" "<<"- YUKAWA_Nu4  (Yukawa coupling of 4th generation neutrino)\n"
     <<setw(width+7)<<" "<<"output of mixing matrices [1=on,0=off(default)]:\n"
     <<setw(width+7)<<" "<<"- OUTPUT_MIXING  (Print the matrices for lepton and quark mixing)\n"
     <<setw(width+4)<<" "<<"}";
}


Standard_Model::Standard_Model(string _dir,string _file,
			       bool _elementary,int _trivialextension) :
  Model_Base(_dir,_file,_elementary), m_trivialextension(_trivialextension)
{
  ParticleInit();
  CustomContainerInit();
  if (m_elementary) {
    ATOOLS::OutputParticles(msg->Info());
    ATOOLS::OutputContainers(msg->Info());
  }
}

bool Standard_Model::ModelInit(const PDF::ISR_Handler_Map& isr)
{
  if (m_elementary) 
    msg_Info()<<"Initialize the Standard Model from "<<m_dir
	      <<" / "<<m_file<<endl;
  m_name      = string("SM");
  p_numbers          = new ScalarNumbersMap();
  p_constants        = new ScalarConstantsMap();
  p_complexconstants = new ComplexConstantsMap();
  p_functions        = new ScalarFunctionsMap();
  p_matrices         = new ComplexMatricesMap();
  
  (*p_numbers)["Extension"] = m_trivialextension;

  FillSpectrum(isr);

  return true;
}

void Standard_Model::ParticleInit() {
  //add SM particles
  //kf_code,mass,width,charge,icharge,strong,spin,majorana,take,stable,massive,idname,tex_name
  s_kftable[1] = new  Particle_Info(1,0.01,.0,-1,-1,3,1,0,1,1,0,"d","d");
  s_kftable[2] = new  Particle_Info(2,0.005,.0,2,1,3,1,0,1,1,0,"u","u");
  s_kftable[3] = new  Particle_Info(3,0.2,.0,-1,-1,3,1,0,1,1,0,"s","s");
  s_kftable[4] = new  Particle_Info(4,1.42,.0,2,1,3,1,0,1,1,0,"c","c");
  s_kftable[5] = new  Particle_Info(5,4.8,.0,-1,-1,3,1,0,1,1,0,"b","b");
  s_kftable[6] = new  Particle_Info(6,175.,1.5,2,1,3,1,0,1,1,1,"t","t");
  s_kftable[11] = new Particle_Info(11,0.000511,.0,-3,-1,0,1,0,1,1,0,"e-","e^-");
  s_kftable[12] = new Particle_Info(12,.0,.0,0,1,0,1,0,1,1,0,"nu_e","\\nu_e");
  s_kftable[13] = new Particle_Info(13,.105,.0,-3,-1,0,1,0,1,1,0,"mu-","\\mu^-");
  s_kftable[14] = new Particle_Info(14,.0,.0,0,1,0,1,0,1,1,0,"nu_mu","\\nu_\\mu");
  s_kftable[15] = new Particle_Info(15,1.777,2.36E-12,-3,-1,0,1,0,1,0,0,"tau-","\\tau^-");
  s_kftable[16] = new Particle_Info(16,.0,.0,0,1,0,1,0,1,1,0,"nu_tau","\\nu_\\tau");
  s_kftable[21] = new Particle_Info(21,.0,.0,0,0,8,2,-1,1,1,0,"G","g");
  s_kftable[22] = new Particle_Info(22,.0,.0,0,0,0,2,-1,1,1,0,"P","\\gamma");
  s_kftable[23] = new Particle_Info(23,91.188,2.49,0,0,0,2,-1,1,1,1,"Z","Z");
  s_kftable[24] = new Particle_Info(24,80.419,2.06,3,0,0,2,0,1,1,1,"W+","W^+");
  s_kftable[25] = new Particle_Info(25,120.,0.0037,0,0,0,0,-1,1,1,1,"h0","h_0");
  if (m_trivialextension==2) {
    s_kftable[7]  = new Particle_Info(7,500.,38.0,-1,-1,3,1,0,1,0,1,"D_4","D_4");
    s_kftable[8]  = new Particle_Info(8,500.,38.2,2,1,3,1,0,1,0,1,"U_4","U_4");
    s_kftable[17] = new Particle_Info(17,500.,38.3,-3,-1,0,1,0,1,0,1,"L_4-","L_4^-");
    s_kftable[18] = new Particle_Info(18,500.,38.3,0,1,0,1,0,1,0,1,"Nu_4","\\nu_4");
  }

  // pseudoparticles for comix
  s_kftable[921] = new Particle_Info(921,0.0,0.0,0,0,8,4,-1,1,1,0,"G4","G_4",1);
  s_kftable[923] = new Particle_Info(923,91.188,2.49,0,0,0,4,-1,1,1,1,"Z4","Z_4",1);
  s_kftable[924] = new Particle_Info(924,80.419,2.06,3,0,0,4,0,1,1,1,"W+4","W^+_4",1);
  s_kftable[925] = new Particle_Info(925,120.0,0.0037,0,0,0,0,0,1,1,1,"h04","h_{04}",1);

  ReadParticleData();

  // add pseudo particles and containers
  s_kftable[kf_none] = new
    Particle_Info(kf_none,-1,0,0,0,0,0,-1,0,1,0,"no_particle","no particle",1,1);
  s_kftable[kf_resummed] = new
    Particle_Info(kf_resummed,0.,0.,0,0,1,2,1,1,1,0,"r","resummed",0,1);
  s_kftable[kf_bjet] = new
    Particle_Info(kf_bjet,0.,0.,0,0,1,2,0,1,1,0,"bj","bjet",1,1);

  s_kftable[kf_fermion] = new
    Particle_Info(kf_fermion,0.,0., 0,0,0,1,0,1,1,0,"fermion","fermion",1,1);
  s_kftable[kf_jet] = new
    Particle_Info(kf_jet,0.,0.,0,0,1, 2,1,1,1,0,"j","jet",1,1);
  s_kftable[kf_quark] = new
    Particle_Info(kf_quark,0.,0.,0, 0,1,1,0,1,1,0,"Q","Quark",1,1);
  s_kftable[kf_lepton] = new
    Particle_Info(kf_lepton,0.,0.,-3,-1,0,1,0,1,1,0,"lepton","lepton",1,1);
  s_kftable[kf_neutrino] = new
    Particle_Info(kf_neutrino,0.,0.,0,1,0, 1,0,1,1,0,"neutrino","neutrino",1,1);
  s_kftable[kf_fermion]->Clear();
  s_kftable[kf_jet]->Clear();
  s_kftable[kf_resummed]->Clear();
  s_kftable[kf_resummed]->m_resummed=true;
  s_kftable[kf_quark]->Clear();
  s_kftable[kf_lepton]->Clear();
  s_kftable[kf_neutrino]->Clear();
  double jet_mass_threshold=p_dataread->GetValue<double>("JET_MASS_THRESHOLD", 500.0);
  for (int i=1;i<7+(m_trivialextension==2?2:0);i++) {
    Flavour addit((kf_code)i);
    if ((addit.Mass()==0.0 || !addit.IsMassive()) && addit.IsOn()) {
      if (addit.Mass(true)<=jet_mass_threshold) {
        s_kftable[kf_jet]->Add(addit);
        s_kftable[kf_jet]->Add(addit.Bar());
        s_kftable[kf_quark]->Add(addit);
        s_kftable[kf_quark]->Add(addit.Bar());
        s_kftable[kf_fermion]->Add(addit);
        s_kftable[kf_fermion]->Add(addit.Bar());
      }
      else {
        PRINT_INFO("Ignoring "<<addit<<" due to JET_MASS_THRESHOLD.");
      }
    }
  }
  s_kftable[kf_jet]->Add(Flavour(kf_gluon));
  s_kftable[kf_jet]->SetResummed();
  for (int i=11;i<17+(m_trivialextension==2?2:0);i+=2) {
    Flavour addit((kf_code)i);
    if ((addit.Mass()==0.0 || !addit.IsMassive()) && addit.IsOn()) {
      s_kftable[kf_lepton]->Add(addit);
      s_kftable[kf_lepton]->Add(addit.Bar());
      s_kftable[kf_fermion]->Add(addit);
      s_kftable[kf_fermion]->Add(addit.Bar());
    }
  }
  for (int i=12;i<17+(m_trivialextension==2?2:0);i+=2) {
    Flavour addit((kf_code)i);
    if ((addit.Mass()==0.0) && addit.IsOn()) {
      s_kftable[kf_neutrino]->Add(addit);
      s_kftable[kf_neutrino]->Add(addit.Bar());
      s_kftable[kf_fermion]->Add(addit);
      s_kftable[kf_fermion]->Add(addit.Bar());
    }
  }
}

void Standard_Model::FillSpectrum(const PDF::ISR_Handler_Map& isr)
{
  p_dataread->RereadInFile();
  FixEWParameters();  
  FixCKM();
  p_constants->insert(make_pair(string("Yukawa_e"), 
				p_dataread->GetValue<double>("YUKAWA_E",0.)));
  p_constants->insert(make_pair(string("Yukawa_mu"), 
				p_dataread->GetValue<double>("YUKAWA_MU",0.)));
  p_constants->insert(make_pair(string("Yukawa_tau"), 
				p_dataread->GetValue<double>("YUKAWA_TAU",Flavour(kf_tau).Mass(true))));
  p_constants->insert(make_pair(string("Yukawa_d"), 
				p_dataread->GetValue<double>("YUKAWA_D",0.)));
  p_constants->insert(make_pair(string("Yukawa_u"), 
				p_dataread->GetValue<double>("YUKAWA_U",0.)));
  p_constants->insert(make_pair(string("Yukawa_s"), 
				p_dataread->GetValue<double>("YUKAWA_S",0.)));
  p_constants->insert(make_pair(string("Yukawa_c"), 
				p_dataread->GetValue<double>("YUKAWA_C",0.)));
  p_constants->insert(make_pair(string("Yukawa_b"), 
				p_dataread->GetValue<double>("YUKAWA_B",Flavour(kf_b).Mass(true))));
  p_constants->insert(make_pair(string("Yukawa_t"), 
				p_dataread->GetValue<double>("YUKAWA_T",Flavour(kf_t).Mass(true))));
  ATOOLS::Flavour(kf_e).SetYuk(ScalarConstant("Yukawa_e"));
  ATOOLS::Flavour(kf_mu).SetYuk(ScalarConstant("Yukawa_mu"));
  ATOOLS::Flavour(kf_tau).SetYuk(ScalarConstant("Yukawa_tau"));
  ATOOLS::Flavour(kf_u).SetYuk(ScalarConstant("Yukawa_u"));
  ATOOLS::Flavour(kf_d).SetYuk(ScalarConstant("Yukawa_d"));
  ATOOLS::Flavour(kf_c).SetYuk(ScalarConstant("Yukawa_c"));
  ATOOLS::Flavour(kf_s).SetYuk(ScalarConstant("Yukawa_s"));
  ATOOLS::Flavour(kf_t).SetYuk(ScalarConstant("Yukawa_t"));
  ATOOLS::Flavour(kf_b).SetYuk(ScalarConstant("Yukawa_b"));

  // Extra coupling parameters for non-Standard tbW coupling
  if (m_trivialextension==1) {
    p_constants->insert(make_pair(string("tbW_RelFactor"),
				  m_trivialextension==1?
				  p_dataread->GetValue<double>("KAPPA_{TBW}",1.):1.));
    p_constants->insert(make_pair(string("tbW_Angle"),
				  m_trivialextension==1?
				  p_dataread->GetValue<double>("THETA_{TBW}",0.):0.));
  }
  if (m_trivialextension==2) { 
    p_constants->insert(make_pair(string("Yukawa_D4"), 
				  p_dataread->GetValue<double>("YUKAWA_D4",Flavour(kf_D4).Mass(true))));
    p_constants->insert(make_pair(string("Yukawa_U4"), 
				  p_dataread->GetValue<double>("YUKAWA_U4",Flavour(kf_U4).Mass(true))));
    p_constants->insert(make_pair(string("Yukawa_L4"), 
				  p_dataread->GetValue<double>("YUKAWA_L4",Flavour(kf_L4).Mass(true))));
    p_constants->insert(make_pair(string("Yukawa_Nu4"), 
				  p_dataread->GetValue<double>("YUKAWA_Nu4",Flavour(kf_Nu4).Mass(true))));
  }

  int    order_alphaS	= p_dataread->GetValue<int>("ORDER_ALPHAS",1);
  int    th_alphaS	= p_dataread->GetValue<int>("THRESHOLD_ALPHAS",1);
  double alphaS         = p_dataread->GetValue<double>("ALPHAS(MZ)",0.118);
  double alphaS_default = p_dataread->GetValue<double>("ALPHAS(default)",alphaS);
  double MZ2            = sqr((*p_constants)[string("MZ")]);

  as = new Running_AlphaS(alphaS,MZ2,order_alphaS,th_alphaS,isr);
  as->SetDefault(alphaS_default);
  p_constants->insert(make_pair(string("alpha_S(MZ)"),alphaS));
  p_functions->insert(make_pair(string("alpha_S"),as));


  double Q2aS      = p_dataread->GetValue<double>("Q2_AS",1.);
  string asf  = p_dataread->GetValue<string>("As_Form",
						       string("smooth"));
  asform::code as_form(asform::smooth);
  if (asf==string("constant"))    as_form = asform::constant;
  else if (asf==string("frozen")) as_form = asform::frozen;
  else if (asf==string("smooth")) as_form = asform::smooth;
  else if (asf==string("IR0"))    as_form = asform::IR0;
  else if (asf==string("GDH"))    as_form = asform::GDH_inspired;
  Strong_Coupling * strong_cpl(new Strong_Coupling(as,as_form,Q2aS));
  p_functions->insert(make_pair(string("strong_cpl"),strong_cpl));

  Running_Fermion_Mass * md   = 
    new Running_Fermion_Mass(Flavour(kf_d),
			     ScalarConstant(string("Yukawa_d")),as);
  Running_Fermion_Mass * mu   = 
    new Running_Fermion_Mass(Flavour(kf_u),
			     ScalarConstant(string("Yukawa_u")),as);
  Running_Fermion_Mass * ms   = 
    new Running_Fermion_Mass(Flavour(kf_s),
			     ScalarConstant(string("Yukawa_s")),as);
  Running_Fermion_Mass * mc   = 
    new Running_Fermion_Mass(Flavour(kf_c),
			     ScalarConstant(string("Yukawa_c")),as);
  Running_Fermion_Mass * mb   = 
    new Running_Fermion_Mass(Flavour(kf_b),
			     ScalarConstant(string("Yukawa_b")),as);
  Running_Fermion_Mass * mt   = 
    new Running_Fermion_Mass(Flavour(kf_t),
			     ScalarConstant(string("Yukawa_t")),as);
  Running_Fermion_Mass * me   = 
    new Running_Fermion_Mass(Flavour(kf_e),
			     ScalarConstant(string("Yukawa_e")),as);
  Running_Fermion_Mass * mmu  = 
    new Running_Fermion_Mass(Flavour(kf_mu),
			     ScalarConstant(string("Yukawa_mu")),as);
  Running_Fermion_Mass * mtau = 
    new Running_Fermion_Mass(Flavour(kf_tau),
			     ScalarConstant(string("Yukawa_tau")),as);
  p_functions->insert(make_pair(string("m")+
				string(Flavour(kf_d).IDName()),md));
  p_functions->insert(make_pair(string("m")+
				string(Flavour(kf_u).IDName()),mu));
  p_functions->insert(make_pair(string("m")+
				string(Flavour(kf_s).IDName()),ms));
  p_functions->insert(make_pair(string("m")+
				string(Flavour(kf_c).IDName()),mc));
  p_functions->insert(make_pair(string("m")+
				string(Flavour(kf_b).IDName()),mb));
  p_functions->insert(make_pair(string("m")+
				string(Flavour(kf_t).IDName()),mt));
  p_functions->insert(make_pair(string("m")+
				string(Flavour(kf_e).IDName()),me));
  p_functions->insert(make_pair(string("m")+
				string(Flavour(kf_mu).IDName()),mmu));
  p_functions->insert(make_pair(string("m")+
				string(Flavour(kf_tau).IDName()),mtau));

  if (m_trivialextension==2) {
    Running_Fermion_Mass * mD4  = 
      new Running_Fermion_Mass(Flavour(kf_D4),
			       ScalarConstant(string("Yukawa_D4")),as);
    Running_Fermion_Mass * mU4  = 
      new Running_Fermion_Mass(Flavour(kf_U4),
			       ScalarConstant(string("Yukawa_U4")),as);
    Running_Fermion_Mass * mL4  = 
      new Running_Fermion_Mass(Flavour(kf_L4),
			       ScalarConstant(string("Yukawa_L4")),as);
    Running_Fermion_Mass * mNu4  = 
      new Running_Fermion_Mass(Flavour(kf_Nu4),
			       ScalarConstant(string("Yukawa_Nu4")),as);
    p_functions->insert(make_pair(string("m")+
				  string(Flavour(kf_D4).IDName()),mD4));
    p_functions->insert(make_pair(string("m")+
				  string(Flavour(kf_U4).IDName()),mU4));
    p_functions->insert(make_pair(string("m")+
				  string(Flavour(kf_L4).IDName()),mL4));
    p_functions->insert(make_pair(string("m")+
				  string(Flavour(kf_Nu4).IDName()),mNu4));
  }
}

void Standard_Model::FixEWParameters() {
  double MW,MZ,MH,GW,GZ,GH,alphaQED0,sin2thetaW,cos2thetaW,vev,lambdaH,GF;
  Complex csin2thetaW, ccos2thetaW, cvev,clambdaH, I(0.0, 1.0);
  string widthscheme = p_dataread->GetValue<string>("WIDTH_SCHEME","Fixed");
  p_numbers->insert(make_pair(string("WidthScheme"), widthscheme=="CMS"));
  int ewscheme = p_dataread->GetValue<int>("EW_SCHEME",0);
  // TODO: Shouldn't we make 1 the default?

  switch (ewscheme) {
  case 0:
    // all SM parameters given explicitly
    MW         = Flavour(kf_Wplus).Mass();
    MZ         = Flavour(kf_Z).Mass();
    MH         = Flavour(kf_h0).Mass();
    GW         = Flavour(kf_Wplus).Width();
    GZ         = Flavour(kf_Z).Width();
    GH         = Flavour(kf_h0).Width();
    alphaQED0  = 1./p_dataread->GetValue<double>("1/ALPHAQED(0)",137.03599976);
    aqed       = new Running_AlphaQED(alphaQED0);
    aqed->SetDefault(1./p_dataread->GetValue<double>("1/ALPHAQED(default)",
                                                     1./(*aqed)(sqr(MZ))));
    sin2thetaW = p_dataread->GetValue<double>("SIN2THETAW",0.23);
    cos2thetaW = 1.-sin2thetaW;
    vev        = p_dataread->GetValue<double>("VEV",246.);
    lambdaH    = p_dataread->GetValue<double>("LAMBDA",0.47591);
    GF         = p_dataread->GetValue<double>("GF",1.16639e-5);
    if (widthscheme=="CMS") {
      THROW(not_implemented, "CMS not implemented for EW_SCHEME=0");
    }
    break;
  case 1:
    // SM parameters given by alphaQED0, M_W, M_Z, M_H
    MW         = Flavour(kf_Wplus).Mass();
    MZ         = Flavour(kf_Z).Mass();
    MH         = Flavour(kf_h0).Mass();
    GW         = Flavour(kf_Wplus).Width();
    GZ         = Flavour(kf_Z).Width();
    GH         = Flavour(kf_h0).Width();
    alphaQED0   = 1./p_dataread->GetValue<double>("1/ALPHAQED(0)",137.03599976);
    aqed       = new Running_AlphaQED(alphaQED0);
    aqed->SetDefault(1./p_dataread->GetValue<double>("1/ALPHAQED(default)",
                                                     1./(*aqed)(sqr(MZ))));
    cos2thetaW = sqr(MW/MZ);
    sin2thetaW = 1.-cos2thetaW;
    vev        = 2.*MW*sqrt(sin2thetaW/(4.*M_PI*aqed->AqedFixed()));
    lambdaH    = 2.*sqr(MH/vev);
    GF         = p_dataread->GetValue<double>("GF",1.16639e-5);
    if (widthscheme=="CMS") {
      Complex muW2(MW*(MW-I*GW)), muZ2(MZ*(MZ-I*GZ)), muH2(MH*(MH-I*GH));
      ccos2thetaW = muW2/muZ2;
      csin2thetaW = 1.-ccos2thetaW;
      cvev        = 2.*sqrt(muW2*csin2thetaW/(4.*M_PI*aqed->AqedFixed()));
      clambdaH    = 2.*muH2/(cvev*cvev);
    }
    break;
  case 2:
    // SM parameters given by alphaQED0, sinthetaW, v, M_H
    alphaQED0   = 1./p_dataread->GetValue<double>("1/ALPHAQED(0)",137.03599976);
    vev        = p_dataread->GetValue<double>("VEV",246.);
    sin2thetaW = p_dataread->GetValue<double>("SIN2THETAW",0.23);
    cos2thetaW = 1.-sin2thetaW;
    MW         = vev/2.*sqrt((4.*M_PI*aqed->AqedFixed())/sin2thetaW);
    MZ         = vev/2.*sqrt((4.*M_PI*aqed->AqedFixed())*(1/sin2thetaW+1/cos2thetaW));
    MH         = Flavour(kf_h0).Mass();
    GW         = Flavour(kf_Wplus).Width();
    GZ         = Flavour(kf_Z).Width();
    GH         = Flavour(kf_h0).Width();
    lambdaH    = 2.*sqr(MH/vev);
    aqed       = new Running_AlphaQED(alphaQED0);
    aqed->SetDefault(1./p_dataread->GetValue<double>("1/ALPHAQED(default)",
                                                     1./(*aqed)(sqr(MZ))));
    GF         = p_dataread->GetValue<double>("GF",1.16639e-5);
    if (widthscheme=="CMS") {
      THROW(not_implemented, "CMS not implemented for EW_SCHEME=2");
    }
    break;
  case 3:
    //gmu scheme
    GF         = p_dataread->GetValue<double>("GF",1.16639e-5);
    MW         = Flavour(kf_Wplus).Mass();
    MZ         = Flavour(kf_Z).Mass();
    MH         = Flavour(kf_h0).Mass();
    GW         = Flavour(kf_Wplus).Width();
    GZ         = Flavour(kf_Z).Width();
    GH         = Flavour(kf_h0).Width();
    sin2thetaW = 1.-sqr(MW/MZ);
    cos2thetaW = 1.-sin2thetaW;
    vev        = 1./(pow(2.,0.25)*sqrt(GF));
    lambdaH    = 2.*sqr(MH/vev); 
    alphaQED0  = 1./p_dataread->GetValue<double>("1/ALPHAQED(0)",137.03599976);
    aqed       = new Running_AlphaQED(alphaQED0);
    aqed->SetDefault(sqrt(2.)*GF/M_PI*sqr(MW)*sin2thetaW);
    if (widthscheme=="CMS") {
      Complex muW2(MW*(MW-I*GW)), muZ2(MZ*(MZ-I*GZ)), muH2(MH*(MH-I*GH));
      ccos2thetaW = muW2/muZ2;
      csin2thetaW = 1.-ccos2thetaW;
      //aqed->SetDefault(Complex(sqrt(2.)*GF/M_PI*muW2*csin2thetaW).real());
      aqed->SetDefault(sqrt(2.)*GF/M_PI*MW*MW*(1-MW*MW/MZ/MZ));
      cvev        = 1./(pow(2.,0.25)*sqrt(GF));
      clambdaH    = 2.*muH2/(cvev*cvev);
      break;
    }
    break;
  default:
    THROW(not_implemented, "Unknown EW_SCHEME="+ToString(ewscheme));
    break;
  }

  Flavour(kf_Wplus).SetMass(MW);
  Flavour(kf_Z).SetMass(MZ);
  Flavour(kf_h0).SetMass(MH);

  p_functions->insert(make_pair(string("alpha_QED"),aqed));
  p_constants->insert(make_pair(string("alpha_QED(0)"),alphaQED0));
  p_constants->insert(make_pair(string("sin2_thetaW"), sin2thetaW));
  p_constants->insert(make_pair(string("cos2_thetaW"), cos2thetaW));
  p_constants->insert(make_pair(string("vev"),         vev));
  p_constants->insert(make_pair(string("MW"),          MW));
  p_constants->insert(make_pair(string("MZ"),          MZ));
  p_constants->insert(make_pair(string("MH"),          MH));
  p_constants->insert(make_pair(string("lambdaH"),     lambdaH));
  p_constants->insert(make_pair(string("GF"),          GF));
  p_numbers->insert(make_pair(string("HIGGS_PARITY"),
			      p_dataread->GetValue<int>("HIGGS_PARITY",1)));

  if (widthscheme=="CMS") {
    p_complexconstants->insert(make_pair(string("ccos2_thetaW"), ccos2thetaW));
    p_complexconstants->insert(make_pair(string("csin2_thetaW"), csin2thetaW));
    p_complexconstants->insert(make_pair(string("cvev"), cvev));
    p_complexconstants->insert(make_pair(string("clambdaH"), clambdaH));
  }
}

void Standard_Model::FixCKM() {
  CMatrix CKM(3);

  for (int i=0;i<3;i++) {
    for (int j=i;j<3;j++) CKM[i][j] = CKM[j][i] = Complex(0.,0.);
    CKM[i][i] = Complex(1.,0.);
  }
  
  double Cabibbo=0.0,A=.8,rho,eta;
  m_ckmorder     = p_dataread->GetValue<int>("CKMORDER",0);  
  if (m_ckmorder>0) {
    Cabibbo    = p_dataread->GetValue<double>("CABIBBO",0.2272);
    CKM[0][0] += sqr(Cabibbo)/2. * Complex(-1.,0.);
    CKM[1][1] += sqr(Cabibbo)/2. * Complex(-1.,0.);
    CKM[0][1] += Cabibbo * Complex( 1.,0.);
    CKM[1][0] += Cabibbo * Complex(-1.,0.);
  }
  if (m_ckmorder>1) {
    A          = p_dataread->GetValue<double>("A",0.818);
    CKM[1][2] += A*sqr(Cabibbo)  * Complex( 1.,0.);
    CKM[2][1] += A*sqr(Cabibbo)  * Complex(-1.,0.);
  }
  if (m_ckmorder>2) {
    eta        = p_dataread->GetValue<double>("ETA",0.349);
    rho        = p_dataread->GetValue<double>("RHO",0.227);
    CKM[0][2] += A*pow(Cabibbo,3) * Complex(rho,-eta);
    CKM[2][0] += A*pow(Cabibbo,3) * Complex(1.-rho,-eta);
  }

  if (m_trivialextension!=2) {
    p_matrices->insert(make_pair(string("CKM"),CKM));

    CMatrix L_CKM(3);   
    for (int i=0;i<3;i++) {
      for (int j=0;j<3;j++) {
	if (i!=j) L_CKM[i][j] = Complex(0.,0.);
	else      L_CKM[i][j] = Complex(1.,0.);
      }
    }
    p_matrices->insert(make_pair(string("L_CKM"),L_CKM));
  }
  else {
    CMatrix V4(4);   
    for (int i=0;i<3;i++) {
      for (int j=0;j<3;j++) {
        V4[i][j] = CKM[i][j];
        V4[3][j] = Complex(0.,0.);
      }
      V4[i][3] = Complex(0.,0.); 
    }
    V4[3][3] = Complex(1.,0.);
    double a14  = p_dataread->GetValue<double>("A_14",0.);
    double a24  = p_dataread->GetValue<double>("A_24",0.);
    double a34  = p_dataread->GetValue<double>("A_34",0.);
    double phi2 = p_dataread->GetValue<double>("PHI_2",0.);
    double phi3 = p_dataread->GetValue<double>("PHI_3",0.);
    Complex I   = Complex(0.,1.);
    CMatrix CKM4(4),Phase14(4),Phase24(4),Phase34(4);
    for (int i=0;i<4;i++) {
      Phase14[i][i] = Phase24[i][i] = Phase34[i][i] = Complex(1.,0.);
      for (int j=i+1;j<4;j++) {
        Phase14[i][j] = Phase14[j][i] =
        Phase24[i][j] = Phase24[j][i] =
        Phase34[i][j] = Phase34[j][i] = Complex(0.,0.);
      }
    }
    // build matrix Phase14
    Phase14[0][0] = Phase14[3][3] = cos(a14);
    Phase14[0][3] = sin(a14)*exp(-I*phi3);
    Phase14[3][0] = -sin(a14)*exp(I*phi3);
    // build matrix Phase24
    Phase24[1][1] = Phase24[3][3] = cos(a24);
    Phase24[1][3] = sin(a24)*exp(-I*phi2);
    Phase24[3][1] = -sin(a24)*exp(I*phi2);
    // build matrix Phase34
    Phase34[2][2] = Phase34[3][3] = cos(a34);
    Phase34[2][3] = sin(a34);
    Phase34[3][2] = -sin(a34);

    // build new CKM4 = CKM*P34*P24*P14
    for (unsigned int i=0;i<4;++i) {
      for (unsigned int j=0;j<4;++j) {
        Complex sum(0.,0.);
        for (unsigned int k=0;k<4;++k) {
          for (unsigned int l=0;l<4;++l) {
            for (unsigned int m=0;m<4;++m) {
              sum += V4[i][k]*Phase34[k][l]*Phase24[l][m]*Phase14[m][j];
            }
          }
        }
        CKM4[i][j] = sum;
      }
    }

    double t14   = p_dataread->GetValue<double>("THETA_L14",0.);
    double t24   = p_dataread->GetValue<double>("THETA_L24",0.);
    double t34   = p_dataread->GetValue<double>("THETA_L34",0.);
    double phiL2 = p_dataread->GetValue<double>("PHI_L2",0.);
    double phiL3 = p_dataread->GetValue<double>("PHI_L3",0.);
    CMatrix L_V4(4),L_CKM4(4),L_Phase14(4),L_Phase24(4),L_Phase34(4);;   
    for (int i=0;i<4;i++) {
      for (int j=0;j<4;j++) {
        if (i!=j) L_V4[i][j] = Complex(0.,0.);
        else      L_V4[i][j] = Complex(1.,0.);
      }
    }
    for (int i=0;i<4;i++) {
      L_Phase14[i][i] = L_Phase24[i][i] = L_Phase34[i][i] = Complex(1.,0.);
      for (int j=i+1;j<4;j++) {
        L_Phase14[i][j] = L_Phase14[j][i] =
        L_Phase24[i][j] = L_Phase24[j][i] =
        L_Phase34[i][j] = L_Phase34[j][i] = Complex(0.,0.);
      }
    }
    // build matrix L_Phase14
    L_Phase14[0][0] = L_Phase14[3][3] = cos(t14);
    L_Phase14[0][3] = sin(t14)*exp(-I*phiL3);
    L_Phase14[3][0] = -sin(t14)*exp(I*phiL3);
    // build matrix L_Phase24
    L_Phase24[1][1] = L_Phase24[3][3] = cos(t24);
    L_Phase24[1][3] = sin(t24)*exp(-I*phiL2);
    L_Phase24[3][1] = -sin(t24)*exp(I*phiL2);
    // build matrix L_Phase34
    L_Phase34[2][2] = L_Phase34[3][3] = cos(t34);
    L_Phase34[2][3] = sin(t34);
    L_Phase34[3][2] = -sin(t34);

    // build new L_CKM4 = L_CKM*L_P34*L_P24*L_P14
    for (unsigned int i=0;i<4;++i) {
      for (unsigned int j=0;j<4;++j) {
        Complex sum(0.,0.);
        for (unsigned int k=0;k<4;++k) {
          for (unsigned int l=0;l<4;++l) {
            for (unsigned int m=0;m<4;++m) {
              sum += L_V4[i][k]*L_Phase34[k][l]*L_Phase24[l][m]*L_Phase14[m][j];
            }
          }
        }
        L_CKM4[i][j] = sum;
      }
    }

    p_matrices->insert(make_pair(string("CKM"),CKM4));
    p_matrices->insert(make_pair(string("L_CKM"),L_CKM4));

    bool output = p_dataread->GetValue<int>("OUTPUT_MIXING",0);
    if (output) {
      unsigned int os(25);
      msg_Out()<<"quark mixing matrix:\n";
      for (int i=0;i<4;++i)
        msg_Out()<<setw(os)<<CKM4[i][0]<<setw(os)<<CKM4[i][1]
                 <<setw(os)<<CKM4[i][2]<<setw(os)<<CKM4[i][3]
                 <<"\n";
      msg_Out()<<"\n";
      msg_Out()<<"lepton mixing matrix:\n";
      for (int i=0;i<4;++i)
        msg_Out()<<setw(os)<<L_CKM4[i][0]<<setw(os)<<L_CKM4[i][1]
                 <<setw(os)<<L_CKM4[i][2]<<setw(os)<<L_CKM4[i][3]
                 <<"\n";
      msg_Out()<<"\n";
    }
  }
}

bool Standard_Model::CheckFlavours(int nin, int nout, Flavour* flavs)
{
  // baryon number
  double bnum(0.);
  for (int i=0;i<nin;i++) bnum-=flavs[i].BaryonNumber();
  for (int i=nin;i<nin+nout;i++) bnum+=flavs[i].BaryonNumber();
  if (!IsZero(bnum)) return false;
  // lepton numbers
  int lnum[3]={0,0,0};
  for (int i=0;i<nin;i++) if (flavs[i].IsLepton()) 
    lnum[flavs[i].LeptonFamily()-1]-=flavs[i].LeptonNumber();
  for (int i=nin;i<nin+nout;i++) if (flavs[i].IsLepton()) 
    lnum[flavs[i].LeptonFamily()-1]+=flavs[i].LeptonNumber();
  for (int i=0;i<3;i++) if (lnum[i]!=0) return false;
  return true;
}
