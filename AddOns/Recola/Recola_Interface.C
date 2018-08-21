#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include "MODEL/Main/Model_Base.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Process/External_ME_Args.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Library_Loader.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include <algorithm>
#include <sys/stat.h>

#include "Recola_Interface.H"
//#include "FMATRIX.h"
using namespace PHASIC;
using namespace MODEL;
using namespace ATOOLS;
using namespace std;

namespace Recola {

  std::string    Recola_Interface::s_recolaprefix     = std::string("");
  double         Recola_Interface::s_light_fermion_threshold=0.1;
  unsigned int   Recola_Interface::s_recolaProcIndex = 0;
  bool           Recola_Interface::s_processesGenerated = false;
  int            Recola_Interface::s_default_flav = 0;
  double         Recola_Interface::s_default_alphaqcd = 0;
  double         Recola_Interface::s_default_scale = 0;
  int            Recola_Interface::s_ewscheme = 3;

  std::vector<double> Recola_Interface::s_pdfmass(6);
  
  std::map<int,PHASIC::Process_Info> m_procmap;
  std::map<int, bool> m_interference;
  int m_doint;

  std::string Recola_Interface::particle2Recola(const int p){
    if(p==1)  return "d";
    if(p==-1) return "d~";
    if(p==2)  return "u";
    if(p==-2) return "u~";
    if(p==3)  return "s";
    if(p==-3) return "s~";
    if(p==4)  return "c";
    if(p==-4) return "c~";
    if(p==5)  return "b";
    if(p==-5) return "b~";
    if(p==6)  return "t";
    if(p==-6) return "t~";

    if(p==11) return "e-";
    if(p==-11)return "e+";
    if(p==12) return "nu_e";
    if(p==-12)return "nu_e~";

    if(p==13) return "mu-";
    if(p==-13)return "mu+";
    if(p==14) return "nu_mu";
    if(p==-14)return "nu_mu~";

    if(p==15) return "tau-";
    if(p==-15)return "tau+";
    if(p==16) return "nu_tau";
    if(p==-16)return "nu_tau~";

    if(p==21) return "g";
    if(p==22) return "A";
    if(p==23) return "Z";
    if(p==24) return "W+";
    if(p==-24)return "W-";
    if(p==25) return "H";
    
    THROW(fatal_error, "Unknown particle id "+ToString(p));
  }

  std::string Recola_Interface::particle2Recola(const std::string p){
    if(p=="d")     return "d";
    if(p=="db")    return "d~";
    if(p=="u")     return "u";
    if(p=="ub")    return "u~";
    if(p=="s")     return "s";
    if(p=="sb")    return "s~";
    if(p=="c")     return "c";
    if(p=="cb")    return "c~";
    if(p=="b")     return "b";
    if(p=="bb")    return "b~";
    if(p=="t")     return "t";
    if(p=="tb")    return "t~";

    if(p=="e-")    return "e-";
    if(p=="e+")    return "e+";
    if(p=="nue")     return "nu_e";
    if(p=="nueb")   return "nu_e~";

    if(p=="mu-")   return "mu-";
    if(p=="mu+")   return "mu+";
    if(p=="numu")   return "nu_mu";
    if(p=="numub")  return "nu_mu~";


    if(p=="tau-")  return "tau-";
    if(p=="tau+")  return "tau+";
    if(p=="nutau")  return "nu_tau";
    if(p=="nutaub") return "nu_tau~";

    if(p=="G")     return "g";
    if(p=="P")     return "A";
    if(p=="Z")     return "Z";
    if(p=="W+")    return "W+";
    if(p=="W-")    return "W-";
    if(p=="h0")    return "H";

    THROW(fatal_error, "Unknown particle id "+ToString(p));
  }

  std::string Recola_Interface::process2Recola(const Flavour_Vector& fl)
  {
    std::string process = particle2Recola(fl[0].IDName())
      + " " + particle2Recola(fl[1].IDName()) + " -> ";
    for(size_t i=2; i<fl.size(); ++i)
      process += particle2Recola(fl[i].IDName())+" ";
    return process;
  }

  bool Recola_Interface::Initialize(const string &path,const string &file,
				    MODEL::Model_Base *const model,
				    BEAM::Beam_Spectra_Handler *const beam,
				    PDF::ISR_Handler *const isr)
  {
    // find RECOLA installation prefix with several overwrite options
    struct stat st;
    Data_Reader reader(" ",";","#","=");

    s_recolaprefix = rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Recola";

    if(stat(s_recolaprefix.c_str(),&st) != 0) s_recolaprefix = RECOLA_PREFIX;
    s_recolaprefix = reader.GetValue<string>("RECOLA_PREFIX", s_recolaprefix);
    msg_Info()<<"Initialising Recola generator from "<<s_recolaprefix<<endl;

    if(MODEL::s_model->Name() != "SM")
      THROW(not_implemented, "ONLY Standard Model so far supported in RECOLA");

    // load library dynamically,
    // keep this at the beginning of the function!
    s_loader->AddPath(s_recolaprefix+"/lib");
    if (!s_loader->LoadLibrary("recola")) 
      THROW(fatal_error, "Failed to load librecola.");

    int recolaVerbosity=reader.GetValue<int>("RECOLA_VERBOSITY",0);
    if(recolaVerbosity<0 || recolaVerbosity >2)
      THROW(fatal_error, "Invalid value for RECOLA_VERBOSITY");
    set_print_level_squared_amplitude_rcl(recolaVerbosity);
    set_print_level_amplitude_rcl(recolaVerbosity);
    set_print_level_correlations_rcl(recolaVerbosity);

    string recolaOutput = reader.GetValue<string>("RECOLA_OUTPUT","*");
    set_output_file_rcl(recolaOutput.c_str());

    int recolaOnShellZW = reader.GetValue<int>("RECOLA_ONSHELLZW",0);

    // set particle masses/widths
    if(recolaOnShellZW != 0){
      set_onshell_mass_z_rcl(Flavour(kf_Z).Mass(),Flavour(kf_Z).Width());
      set_onshell_mass_w_rcl(Flavour(kf_Wplus).Mass(),Flavour(kf_Wplus).Width());
    }
    else{
      set_pole_mass_z_rcl(Flavour(kf_Z).Mass(),Flavour(kf_Z).Width());
      set_pole_mass_w_rcl(Flavour(kf_Wplus).Mass(),Flavour(kf_Wplus).Width());
    }
    set_pole_mass_h_rcl(Flavour(kf_h0).Mass(),Flavour(kf_h0).Width());
    set_pole_mass_electron_rcl(Flavour(kf_e).Mass());
    set_pole_mass_muon_rcl(Flavour(kf_mu).Mass(),Flavour(kf_mu).Width());
    set_pole_mass_tau_rcl(Flavour(kf_tau).Mass(),Flavour(kf_tau).Width());
    set_pole_mass_up_rcl(Flavour(kf_u).Mass());
    set_pole_mass_down_rcl(Flavour(kf_d).Mass());
    set_pole_mass_strange_rcl(Flavour(kf_s).Mass());
    set_pole_mass_charm_rcl(Flavour(kf_c).Mass(),Flavour(kf_c).Width());
    set_pole_mass_bottom_rcl(Flavour(kf_b).Mass(),Flavour(kf_b).Width());
    set_pole_mass_top_rcl(Flavour(kf_t).Mass(),Flavour(kf_t).Width());
    s_light_fermion_threshold = reader.GetValue<double>("RECOLA_LIGHT_FERMION_THRESHOLD",0);
    set_light_fermions_rcl(s_light_fermion_threshold);
    set_delta_ir_rcl(0.0,M_PI*M_PI/6.0); // adapts the conventions from COLLIER to Catani-Seymour
    
    PDF::PDF_Base *pdf=isr->PDF(0);
    int pdfnf=pdf->ASInfo().m_flavs.size();
    s_default_alphaqcd=pdf->ASInfo().m_asmz;
    s_default_scale=pdf->ASInfo().m_mz2;
    s_default_flav=pdfnf;
    if (pdfnf>10) pdfnf-=10;
    if (pdfnf==-1) pdfnf=6;

    double cmass=pdf->ASInfo().m_flavs[3].m_mass;
    double bmass=pdf->ASInfo().m_flavs[4].m_mass;
    double tmass=pdf->ASInfo().m_flavs[5].m_mass;

    cmass=reader.GetValue<double>("RECOLA_AS_RUN_MASS_C", cmass);
    cmass=reader.GetValue<double>("RECOLA_AS_REN_MASS_C", cmass);
    bmass=reader.GetValue<double>("RECOLA_AS_RUN_MASS_B", bmass);
    bmass=reader.GetValue<double>("RECOLA_AS_REN_MASS_B", bmass);
    tmass=reader.GetValue<double>("RECOLA_AS_RUN_MASS_T", tmass);
    tmass=reader.GetValue<double>("RECOLA_AS_REN_MASS_T", tmass);

    for (int i=0; i<3; i++){
      if (i<pdfnf)
	s_pdfmass[i]=pdf->ASInfo().m_flavs[i].m_thres;
    }
    s_pdfmass[3]=cmass;
    s_pdfmass[4]=bmass;
    s_pdfmass[5]=tmass;

    set_alphas_masses_rcl(cmass,bmass,tmass,Flavour(kf_c).Width(),Flavour(kf_b).Width(),Flavour(kf_t).Width()); 

    // Init EW scheme
    s_ewscheme = ToType<int>(rpa->gen.Variable("EW_SCHEME"),1);
    double aqed = MODEL::aqed->Default();
    switch (s_ewscheme)
      {
      case 3:
	use_gfermi_scheme_and_set_alpha_rcl(aqed);
	break;
      case 2:
	//Currently not implemented in Recola, GFermi is used instead!
        use_alphaz_scheme_rcl();
        break;
      case 1:
	use_alpha0_scheme_rcl();
	break;
      default:
	THROW(fatal_error, "Unsupported EW scheme");
      }

    // Set IR and UV scales 
    double IRscale=reader.GetValue<double>("IR_SCALE",100);
    double UVscale=reader.GetValue<double>("UV_SCALE",100);
    set_mu_ir_rcl(IRscale);
    set_mu_uv_rcl(UVscale);

    // Used to be done in both Recola_Born and Recola_Virtual.
    // This just sets s_default_flav.
    // Can for sure be simplified and made more transparent
    int nlight=0;
    int fixed=reader.GetValue<int>("RECOLA_FIXED_FLAVS",Recola_Interface::GetDefaultFlav()+10);
    int default_flavscheme(fixed);
    if (default_flavscheme==16) default_flavscheme=-1;
    if (fixed>0 && fixed<10){
      nlight=fixed;
    }
    else{
      if (default_flavscheme>10){
	nlight=Recola_Interface::PDFnf(sqr(Flavour(kf_Z).Mass()),default_flavscheme-10);
      }
      if (default_flavscheme==-1)
	nlight=-1;
      if (default_flavscheme==-2 || default_flavscheme==0){
	if (Flavour(kf_c).Mass()!=0)
	  nlight=3;
	else if (Flavour(kf_b).Mass()!=0)
	  nlight=4;
	else if (Flavour(kf_t).Mass()!=0)
	  nlight=5;
	else {
	  msg_Out()<<"WARNING: 6 light flavours detected.\n";
	  nlight=6;
	}
      }
    }
    s_default_flav=nlight;
    
    // Set AlphaS
    set_alphas_rcl(MODEL::as->Default(),
		   Flavour(kf_Z).Mass(),
		   s_default_flav);

    return true;
  }

  int Recola_Interface::RegisterBorn(const External_ME_Args& args,
				     const int& amptype)
  {
    /* Assign index to process */
    increaseProcIndex(); const int& procIndex(getProcIndex());
    string procstring(process2Recola(args.Flavours()));
    define_process_rcl(procIndex,
		       procstring.c_str(),
		       "NLO");

    /* TODO: set properly */
    m_interference[procIndex]=false;

    /* Set coupling orders */
    unselect_all_gs_powers_BornAmpl_rcl(procIndex);
    unselect_all_gs_powers_LoopAmpl_rcl(procIndex);
    if (amptype==12)
      select_gs_power_LoopAmpl_rcl(procIndex,args.m_orders[0]);
    else
      select_gs_power_BornAmpl_rcl(procIndex,args.m_orders[0]);
       
    return procIndex;
  }
      
  int Recola_Interface::RegisterProcess(const Process_Info& pi,
					int amptype)
  {

    increaseProcIndex();
    msg_Debugging()<<"Recola_Interface::RegisterProcess called\n";
    int procIndex(getProcIndex());
    msg_Debugging()<<"ProcIndex = " <<procIndex <<"\n"; 
    msg_Debugging()<<"process string = "<<process2Recola(pi)<<"\n";
    m_procmap[procIndex]=pi;
    if (!pi.m_nlomode && amptype!=12) {
      msg_Debugging() << "no NLO mode detected!\n";
      return 0;
    }
    // set procIndex to map with flavours
    define_process_rcl(procIndex,process2Recola(pi).c_str(),"NLO");
    unselect_all_gs_powers_BornAmpl_rcl(procIndex);
    unselect_all_gs_powers_LoopAmpl_rcl(procIndex);
    m_interference[procIndex]=false;
    Data_Reader reader(" ",";","#","=");
    reader.AddIgnore("[");
    reader.AddIgnore("]");
    m_doint=reader.GetValue<int>("RECOLA_INTERFERENCE",1);
    int quarkcount(0), gluoncount(0);
    int tempQCD(pi.m_maxcpl[0]), tempEW(pi.m_maxcpl[1]);
    if(pi.m_fi.m_nlotype==nlo_type::loop){
      if (m_doint){
	Flavour_Vector inflavs(pi.m_ii.GetExternal());
	Flavour_Vector outflavs(pi.m_fi.GetExternal());
	for (int i=0; i<inflavs.size(); i++){
	  if (inflavs[i].IsQuark())
	    quarkcount++;
	  else if (inflavs[i].IsGluon())
	    gluoncount++;
	}
	for (int i=0; i<outflavs.size(); i++){
	  if (outflavs[i].IsQuark())
	    quarkcount++;
	  else if (outflavs[i].IsGluon())
	    gluoncount++;
	}
	tempQCD-=gluoncount;
	if ((pi.m_fi.m_nlocpl[1]==1) && (quarkcount>=4) && (pi.m_maxcpl[0]>=2)){
	  m_interference[procIndex]=true;
	}
	if ((pi.m_fi.m_nlocpl[0]==1) && (quarkcount>=4) && (pi.m_maxcpl[1]>=2)){	  
	  m_interference[procIndex]=true;
	  tempQCD-=1;
	}
	tempEW=quarkcount-2-tempQCD;
      }

      if (m_interference[procIndex]){
	int maxBqcd, minBqcd;
	
	if ((tempQCD+2*pi.m_fi.m_nlocpl[0])>tempEW)
	  maxBqcd=pi.m_maxcpl[0]+tempEW-pi.m_fi.m_nlocpl[0];
	else
	  maxBqcd=pi.m_maxcpl[0]+tempQCD+pi.m_fi.m_nlocpl[0];
	if ((tempEW+2*pi.m_fi.m_nlocpl[1])>tempQCD)
	  minBqcd=pi.m_maxcpl[0]-tempQCD-pi.m_fi.m_nlocpl[0];
	else
	  minBqcd=pi.m_maxcpl[0]-tempEW-pi.m_fi.m_nlocpl[0]-2*pi.m_fi.m_nlocpl[1];

	while (quarkcount>=2 && maxBqcd>=minBqcd){
	  select_gs_power_LoopAmpl_rcl(procIndex,2.*pi.m_maxcpl[0]-maxBqcd);
	  select_gs_power_BornAmpl_rcl(procIndex,maxBqcd);
	  quarkcount-=2;
	  maxBqcd-=2;
	}
      }
      
      else{
	select_gs_power_BornAmpl_rcl(procIndex,pi.m_maxcpl[0]-pi.m_fi.m_nlocpl[0]);
	select_gs_power_LoopAmpl_rcl(procIndex,pi.m_maxcpl[0]+pi.m_fi.m_nlocpl[0]);
      }
    }
    else if (amptype==12) {
      select_gs_power_LoopAmpl_rcl(procIndex,pi.m_maxcpl[0]);
    }
    msg_Debugging()<<"procIndex "<<procIndex<<" returned\n";  
    return procIndex;
  }
  
  void Recola_Interface::EvaluateLoop(int id, const Vec4D_Vector& momenta, double& bornres, METOOLS::DivArrD& virt)
  {
    const int NN = momenta.size();
    double fpp[NN][4];
    
    for (int i=0; i<NN; i++){
      for (int mu=0; mu<4; mu++){
	fpp[i][mu] = momenta[i][mu];
      }
    }
    double fA2[2]={0.0};
    
    bool momcheck(0);
    //    compute_process_rcl(id,fpp,NN,"NLO",fA2);
    compute_process_rcl(id,fpp,"NLO",fA2,momcheck); // Change discussed in meeting. Mathieu 12/04/2017
    int procIndex(id);
    PHASIC::Process_Info pi(m_procmap[id]);
    
    // TODO: had to change get_squared_amplitude_rcl to
    // get_squared_amplitude_r1_rcl, is this correct?
    if (m_interference[procIndex]){
      get_squared_amplitude_rcl(id,pi.m_maxcpl[0]-pi.m_fi.m_nlocpl[0],"LO",fA2[0]);
      get_squared_amplitude_rcl(id,pi.m_maxcpl[0],"NLO",fA2[1]);
    }
 
    bornres = fA2[0];
    virt.Finite()=fA2[1];
  }
  
  int Recola_Interface::PDFnf(double scale, int maxn){
    int nf(0);
    for (int i=0; i<=maxn; i++){
      nf=i;
      if (sqrt(scale)<s_pdfmass[i])
	break;
    }
    return nf;
  }

  int Recola_Interface::PerformTests()
  {
    return 1;
  }
  
  void Recola_Interface::PrepareTerminate()
  {
    
  }


}

using namespace Recola;

DECLARE_GETTER(Recola_Interface,"Recola",ME_Generator_Base,ME_Generator_Key);

ME_Generator_Base *ATOOLS::Getter<ME_Generator_Base,ME_Generator_Key,
                                  Recola_Interface>::
operator()(const ME_Generator_Key &key) const
{
  return new Recola::Recola_Interface();
}

void ATOOLS::Getter<ME_Generator_Base,ME_Generator_Key,Recola_Interface>::
PrintInfo(ostream &str,const size_t width) const
{
  str<<"Interface to the Recola loop ME generator";
}
