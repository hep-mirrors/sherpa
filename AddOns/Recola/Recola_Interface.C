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

using namespace PHASIC;
using namespace MODEL;
using namespace ATOOLS;
using namespace std;


namespace Recola {

  std::string    Recola_Interface::s_recolaprefix = std::string("");
  bool           Recola_Interface::s_ignore_model = false;
  bool           Recola_Interface::s_exit_on_error= true;
  double         Recola_Interface::s_light_fermion_threshold=0.1;
  size_t         Recola_Interface::s_recolaProcIndex = 0;
  bool           Recola_Interface::s_processesGenerated = false;
  size_t         Recola_Interface::s_getPDF_default = 0;
  size_t         Recola_Interface::s_default_flav = 6;
  size_t         Recola_Interface::s_fixed_flav = 0;
  double         Recola_Interface::s_default_alphaqcd = 0;
  double         Recola_Interface::s_default_scale = 0;
  std::vector<double> Recola_Interface::s_pdfmass(6);
  bool           Recola_Interface::s_compute_poles = false;
  size_t         Recola_Interface::s_vmode = 0;
  int            Recola_Interface::s_ewscheme = 3;
  int            Recola_Interface::s_amptype = 1;
  bool           Recola_Interface::s_mass_reg = 0;
  double         Recola_Interface::s_photon_mass = 0.1;
  bool           Recola_Interface::s_check_mass = 0;
  bool           Recola_Interface::s_self_energy_off = 0;
  bool           Recola_Interface::s_use_decay = 0;
  bool           Recola_Interface::s_qed_only = 0;
  bool           Recola_Interface::s_ew_only = 0;

  
  std::map<size_t,PHASIC::Process_Info> Recola_Interface::s_procmap;
  std::map<size_t,ATOOLS::asscontrib::type> Recola_Interface::s_asscontribs;
  
  std::map<size_t, bool> Recola_Interface::s_interference;
  size_t Recola_Interface::s_doint = 1;

  Recola_Interface::Recola_Interface() :
        ME_Generator_Base("Recola") {RegisterDefaults();}

  Recola_Interface::~Recola_Interface() {}



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
    if(p=="ve")     return "nu_e";
    if(p=="veb")   return "nu_e~";

    if(p=="mu-")   return "mu-";
    if(p=="mu+")   return "mu+";
    if(p=="vmu")   return "nu_mu";
    if(p=="vmub")  return "nu_mu~";


    if(p=="tau-")  return "tau-";
    if(p=="tau+")  return "tau+";
    if(p=="vtau")  return "nu_tau";
    if(p=="vtaub") return "nu_tau~";

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
  

  void Recola_Interface::RegisterDefaults() const
  {
    Scoped_Settings s{ Settings::GetMainSettings()["Recola"] };
    Scoped_Settings spol{ Settings::GetMainSettings()["POLARIZATIONS"] };
    s["VERBOSITY"].SetDefault(0);
    s["IGNORE_MODEL"].SetDefault(0);
    s["EXIT_ON_ERROR"].SetDefault(1);
    s["USE_I_IN_EWAPPROX"].SetDefault(false);
    s["GETPDF_DEFAULT"].SetDefault(0);
    s["IR_SCALE"].SetDefault(100);
    s["UV_SCALE"].SetDefault(100);
    s["LIGHT_FERMION_THRESHOLD"].SetDefault(1e-20);
    s["OUTPUT"].SetDefault("*");
    s["INTERFERENCE"].SetDefault(1);
    s["ONSHELLZW"].SetDefault(0);
    s["COMPUTE_POLES"].SetDefault(0);
    s["COLLIER_CACHE"].SetDefault(-1);
    s["VMODE"].SetDefault(0);
    s["AMPTYPE"].SetDefault(1);
    s["PHOTON_MASS"].SetDefault(0.1);
    s["USE_DECAY"].SetDefault(0);
    s["MASS_REG"].SetDefault(false);
    s["NO_SELF_ENERGY"].SetDefault(false);
    s["COMPLEX_MASS_SCHEME"].SetDefault(true);
    s["QED_ONLY"].SetDefault(false);
    s["EW_ONLY"].SetDefault(false);
    // find RECOLA installation prefix with several overwrite options
    char *var=NULL;
    s_recolaprefix = rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Recola";
    s_recolaprefix = string(((var=getenv("RECOLA_PREFIX"))==NULL ? s_recolaprefix : var));
    struct stat st;
    if(stat(s_recolaprefix.c_str(), &st) != 0)
        s_recolaprefix = RECOLA_PREFIX;
    s["PREFIX"].SetDefault(s_recolaprefix);
    s_recolaprefix = s["PREFIX"].Get<string>();
  }



  bool Recola_Interface::Initialize(MODEL::Model_Base *const model,
            BEAM::Beam_Spectra_Handler *const beam,
            PDF::ISR_Handler *const isr,
            YFS::YFS_Handler *const yfs)
  {
    msg_Info()<<"Initialising Recola generator from "<<s_recolaprefix<<endl;
    Scoped_Settings s{ Settings::GetMainSettings()["Recola"] };

    //This check could be added if a different model from the SM wants to be used
    s_ignore_model = s["IGNORE_MODEL"].Get<bool>();
    if (s_ignore_model) {
      msg_Info()<<METHOD<<"(): Recola will use the "
                        <<"Standard Model even if you set a "
                        <<"different model without warning."
                        <<std::endl;
    }

    // VERBOSITY
    int recolaVerbosity=s["VERBOSITY"].Get<int>();
    if(recolaVerbosity<0 || recolaVerbosity >2)
      THROW(fatal_error, "Invalid value for RECOLA_VERBOSITY");
    set_print_level_squared_amplitude_rcl(recolaVerbosity);
    set_print_level_amplitude_rcl(recolaVerbosity);
    set_print_level_correlations_rcl(recolaVerbosity);

    string recolaOutput = s["OUTPUT"].Get<std::string>();
    s_amptype           = s["AMPTYPE"].Get<int>();
    set_output_file_rcl(recolaOutput.c_str());
    s_vmode = s["VMODE"].Get<int>();
    msg_Tracking()<<METHOD<<"(): Set V-mode to "<<s_vmode<<endl;
    s_photon_mass = s["PHOTON_MASS"].Get<double>();
    s_use_decay   = s["USE_DECAY"].Get<bool>();
    s_qed_only    = s["QED_ONLY"].Get<bool>();
    s_ew_only     = s["EW_ONLY"].Get<bool>();

    // if(s_qed_only) set_qed_only_rcl();
    // if(s_ew_only) set_ew_only_rcl();
    if(s_qed_only==1 && s_ew_only==1) THROW(fatal_error, "Only one of 'QED_ONLY' and 'EW_ONLY' can be enabled at a time");
    if (s_vmode&2) THROW(fatal_error,"Inclusion of I operator not implemented.");

    // Compute poles
    int cp(0);
    s_compute_poles = s["COMPUTE_POLES"].Get<int>();
    if (s_compute_poles) {
      msg_Info()<<METHOD<<"(): Instructing Recola to compute poles."<<std::endl;
    }


    if(MODEL::s_model->Name() != "SM")
      THROW(not_implemented, "ONLY Standard Model so far supported in RECOLA");
    


    bool recolaOnShellZW = s["ONSHELLZW"].Get<bool>();
    bool recolaComplexMass = s["COMPLEX_MASS_SCHEME"].Get<bool>();
    if(!recolaComplexMass) set_on_shell_scheme_rcl();
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
    s_light_fermion_threshold = s["LIGHT_FERMION_THRESHOLD"].Get<double>();
    set_light_fermions_rcl(s_light_fermion_threshold);
    set_delta_ir_rcl(0.0,M_PI*M_PI/6.0); // adapts the conventions from COLLIER to Catani-Seymour
    
    PDF::PDF_Base *pdf=isr->PDF(0);
    auto pdfnf = -1;
    auto cmass = 0.0;
    auto bmass = 0.0;
    auto tmass = 0.0;
    bool hadronic_beam1 = beam->GetBeam(0)->Beam().IsHadron(); 
    bool hadronic_beam2 = beam->GetBeam(1)->Beam().IsHadron(); 
    if(hadronic_beam1!=hadronic_beam2) THROW(not_implemented,"Recola interface cannot handle DIS yet.");
    bool hadronic_beam = hadronic_beam1;

    if (hadronic_beam) {
      pdfnf=pdf->ASInfo().m_flavs.size();
      s_default_alphaqcd=pdf->ASInfo().m_asmz;
      s_default_scale=pdf->ASInfo().m_mz2;
      s_default_flav=pdfnf;
      if (pdfnf>10) pdfnf-=10;
      if (pdfnf==-1) pdfnf=6;
      cmass=pdf->ASInfo().m_flavs[3].m_mass;
      bmass=pdf->ASInfo().m_flavs[4].m_mass;
      tmass=pdf->ASInfo().m_flavs[5].m_mass;
      if(pdf->ASInfo().m_flavs.size()<6) tmass=Flavour(kf_t).Mass();
    } 
      else {
      pdfnf = MODEL::as->Nf(1.e20);
      s_default_alphaqcd=MODEL::as->AsMZ();
      s_default_scale=Flavour{kf_Z}.Mass();
      s_default_flav=pdfnf;
      const auto thresholds = MODEL::as->Thresholds(0.0, 1e20);
      tmass = sqrt(thresholds[thresholds.size() - 1]);
      bmass = sqrt(thresholds[thresholds.size() - 2]);
      cmass = sqrt(thresholds[thresholds.size() - 3]);
    }


    if (hadronic_beam) {
      for (int i=0; i<3; i++){
        if (i<pdfnf) s_pdfmass[i]=pdf->ASInfo().m_flavs[i].m_thres;
      }
    } 
    else {
      for (int i{0}; i < 3; ++i) s_pdfmass[i] = Flavour{i+1}.Mass(1);
    }

    s_pdfmass[3]=cmass;
    s_pdfmass[4]=bmass;
    s_pdfmass[5]=tmass;
    set_alphas_masses_rcl(cmass,bmass,tmass,
                          Flavour(kf_c).Width(),Flavour(kf_b).Width(),
                          Flavour(kf_t).Width()); 
    
    s_fixed_flav=Recola_Interface::GetDefaultFlav()+10;
    s_mass_reg = s["MASS_REG"].Get<bool>();
    s_check_mass = s["MASS_REG"].Get<bool>();
    s_self_energy_off = s["NO_SELF_ENERGY"].Get<bool>();
    if(s_self_energy_off) switchon_resonant_selfenergies_rcl();
    
    if(s_mass_reg){
      s_photon_mass = s["PHOTON_MASS"].Get<double>();
      if(s_photon_mass != yfs->m_photonMass){
        msg_Error()<<"Mismatch between YFS and Recola photon mass"
                   <<"\n mass in YFS = "<<yfs->m_photonMass
                   <<"\n mass in RECOLA = "<<s_photon_mass<<endl;
        THROW(fatal_error,"Mismatch in photon mass regulator");
    }
    
    set_dynamic_settings_rcl(1);
    use_mass_reg_soft_rcl(s_photon_mass);
  }
    return true;
  }



  // This function is specific for LO or Born processes since they read 
  // settings from External_ME_Args... 
  int Recola_Interface::RegisterProcess(const External_ME_Args& args,
             const int& amptype)
  {
    DEBUG_FUNC("");
    increaseProcIndex();
    msg_Debugging()<<"Recola_Interface::RegisterProcess called\n";
    int procIndex(getProcIndex());
    msg_Debugging()<<"ProcIndex = " <<procIndex <<"\n"; 
    msg_Debugging()<<"process string = "<<process2Recola(args.Flavours())<<"\n";
    
      
    string procstring(process2Recola(args.Flavours()));
    if(amptype==1) define_process_rcl(procIndex, procstring.c_str(), "LO");
    else define_process_rcl(procIndex, procstring.c_str(), "NLO");
    
    
    Scoped_Settings s{ Settings::GetMainSettings()["Recola"] };

    /* Set coupling orders */
    unselect_all_gs_powers_BornAmpl_rcl(procIndex);
    unselect_all_gs_powers_LoopAmpl_rcl(procIndex);
    if (amptype==12)
    {
      // set collier caching level
      int cc=s["COLLIER_CACHE"].Get<int>();
      if (cc>=0) split_collier_cache_rcl(procIndex,cc);
      select_gs_power_LoopAmpl_rcl(procIndex,args.m_orders[0]);
    }
    else
      select_gs_power_BornAmpl_rcl(procIndex,args.m_orders[0]);
       
    return procIndex;
  }
      
  size_t Recola_Interface::RegisterProcess(const Process_Info& pi,
          int amptype)
  {
    DEBUG_FUNC("");
    std::string decayprocess, myprocess;
    if(s_use_decay) {
      std::string decayprocess = particle2Recola(pi.m_ii.m_ps[0].m_fl.IDName())
                    + " " + particle2Recola(pi.m_ii.m_ps[1].m_fl.IDName()) + " -> ";
      for(auto dec: pi.m_fi.m_ps ){
        PRINT_VAR(dec.m_ps.size());
        if(dec.m_ps.size() != 0){
          decayprocess+= particle2Recola(dec.m_fl.IDName());
          decayprocess+="(";
          for (int i = 0; i < dec.m_ps.size(); ++i)
          {
            decayprocess+= particle2Recola(dec.m_ps[i].m_fl);
            if(i==0) decayprocess+=" ";
          }
          decayprocess+=")";
        }
        else{
          decayprocess+= particle2Recola(dec.m_fl.IDName())+" ";
        }
        myprocess=decayprocess;
      }
    }
    increaseProcIndex();
    msg_Debugging()<<"Recola_Interface::RegisterProcess called\n";
    int procIndex(getProcIndex());
    msg_Debugging()<<"ProcIndex = " <<procIndex <<"\n"; 
    msg_Debugging()<<"process string = "<<process2Recola(pi)<<"\n";
    
    // set procIndex to map with flavours
    s_procmap[procIndex]=pi;
    if (pi.m_nlomode && amptype==12) {
      msg_Debugging() << "Recola cannot do looploop NLO!\n";
      return 0;
    }

    // define process in Recola, at this stage always 'NLO'
    // if(s_use_decay) define_process_rcl(procIndex,myprocess.c_str(),"NLO");
    // else define_process_rcl(procIndex,process2Recola(pi).c_str(),"NLO");
    if(s_use_decay){
      if(amptype==-1) define_process_rcl(procIndex, myprocess.c_str(), "LO");
      else define_process_rcl(procIndex,myprocess.c_str(),"NLO");
    }
    else{
      if(s_use_decay) define_process_rcl(procIndex,myprocess.c_str(),"NLO");
      else define_process_rcl(procIndex,process2Recola(pi).c_str(),"NLO");
    }
    s_interference[procIndex]=false;
    Scoped_Settings s{ Settings::GetMainSettings()["Recola"] };
    
    s_doint=s["INTERFERENCE"].Get<int>();
    
    // set collier caching level
    int cc=s["COLLIER_CACHE"].Get<int>();
    if (cc>=0) split_collier_cache_rcl(procIndex,cc);

    // find out whether we need multiple orders or not
    s_asscontribs[procIndex]=pi.m_fi.m_asscontribs;
    
    // if we only need specific orders, set them
    if (s_asscontribs[procIndex]==asscontrib::none) {
      // unset all powers of the amplitude
      unselect_all_gs_powers_BornAmpl_rcl(procIndex);
      unselect_all_gs_powers_LoopAmpl_rcl(procIndex);
      
      int quarkcount(0), gluoncount(0);
      int tempQCD(pi.m_maxcpl[0]), tempEW(pi.m_maxcpl[1]);
      if(pi.m_fi.m_nlotype==nlo_type::loop){
        
        // Check whether for this process any interference 
        // diagram is present
        if (s_doint){
          Flavour_Vector inflavs(pi.m_ii.GetExternal());
          Flavour_Vector outflavs(pi.m_fi.GetExternal());
            for (int i=0; i<inflavs.size(); i++){
              if (inflavs[i].IsQuark()) quarkcount++;
              else if (inflavs[i].IsGluon()) gluoncount++;
            }

          for (int i=0; i<outflavs.size(); i++){
            if (outflavs[i].IsQuark()) quarkcount++;
            else if (outflavs[i].IsGluon()) gluoncount++;
          }
          tempQCD-=gluoncount;
          if ((pi.m_fi.m_nlocpl[1]==1) && (quarkcount>=4) && (pi.m_maxcpl[0]>=2)){
            s_interference[procIndex]=true;
          }
          if ((pi.m_fi.m_nlocpl[0]==1) && (quarkcount>=4) && (pi.m_maxcpl[1]>=2)){    
            s_interference[procIndex]=true;
            tempQCD-=1;
          }
          tempEW=quarkcount-2-tempQCD;
        }

        // If interference is present set orders properly
        if (s_interference[procIndex]){
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

        // If there is no interference set orders normally checking for 
        // EW or QCD NLO
        else{
          // now set the requested powers for the amplitude
          if (pi.m_fi.m_nlocpl[0]==1) {
            double borngspower=pi.m_maxcpl[0]-pi.m_fi.m_nlocpl[0];
            double loopgspower=pi.m_maxcpl[0]+pi.m_fi.m_nlocpl[0];
            
            msg_Debugging()<<"QCD Tree gs-power = "<<borngspower<<std::endl
                         <<"    Loop gs-power = "<<loopgspower<<std::endl;
            
            select_gs_power_BornAmpl_rcl(procIndex,borngspower);
            select_gs_power_LoopAmpl_rcl(procIndex,loopgspower);
          }
        
          else if (pi.m_fi.m_nlocpl[1]==1) {
            msg_Debugging()<<"EW Tree gs-power = "<<pi.m_maxcpl[0]<<std::endl
                           <<"   Loop gs-power = "<<pi.m_maxcpl[0]<<std::endl;
            select_gs_power_BornAmpl_rcl(procIndex,pi.m_maxcpl[0]);
            select_gs_power_LoopAmpl_rcl(procIndex,pi.m_maxcpl[0]);
          }
        }
      }
      else {
        msg_Debugging()<<"Born gs-power = "<<pi.m_maxcpl[0]<<std::endl;
        if(amptype==12)
          select_gs_power_LoopAmpl_rcl(procIndex,pi.m_maxcpl[0]);
        else
          select_gs_power_BornAmpl_rcl(procIndex,pi.m_maxcpl[0]);
      }
    }
    
    // If ass contributions are needed just initialize with every power
    else {
      std::cout<<"initialized with every power\n";
      msg_Debugging()<<"Initialise Tree and Loop with all gs-powers"<<std::endl;
    }
    msg_Debugging()<<"procIndex "<<procIndex<<" returned\n";
    return procIndex;
  }
  

  void Recola_Interface::GenerateProcesses(const double& alpha, const double& alphas,
                                           const double& muIR, const double& muUV,
                                           const double& muR2)
  {
    DEBUG_FUNC("");
    
    Settings& s = Settings::GetMainSettings();
    ew_scheme::code ewscheme = s["EW_SCHEME"].Get<ew_scheme::code>();
    ew_scheme::code ewrenscheme = s["EW_REN_SCHEME"].Get<ew_scheme::code>();
    if (ewscheme!=ewrenscheme) THROW(fatal_error,"Inconsistent input scheme.");
    switch (ewscheme) {
      case 1:
        // use_alpha0_scheme_and_set_alpha_rcl(AlphaQED());
        use_alpha0_scheme_rcl(alpha);
        break;
      case 2:
        // use_alphaz_scheme_and_set_alpha_rcl(AlphaQED());
        use_alphaz_scheme_rcl(alpha);
        break;
      case 3:
        use_gfermi_scheme_and_set_alpha_rcl(alpha);
        break;
      default:
        msg_Error()<<"The EW scheme "<<ewscheme<<" is not available with the "
                   <<"Sherpa+Recola interface. Valid options are:\n"
                   <<"  1) alpha(0)\n"
                   <<"  2) alpha(M_Z)\n"
                   <<"  3) GFermi"<<std::endl;
        THROW(fatal_error,"Unknown EW_SCHEME setting.");
    }


    int nlight=0;
    set_mu_ir_rcl(muIR);
    set_mu_uv_rcl(muUV);
    size_t fflav(Recola_Interface::GetFixedFlav());
    double alpha_mat;
    int default_flavscheme(fflav);
    if (default_flavscheme==16) default_flavscheme=-1;
    if (fflav>0 && fflav<10) nlight=fflav;
    else {
      if (default_flavscheme>10)
        nlight=Recola_Interface::PDFnf(muR2,default_flavscheme-10);
      if (default_flavscheme==-1)
        nlight=-1;
      if (default_flavscheme==-2 || default_flavscheme==0) {
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
    if (nlight==0) {
      msg_Error()<<METHOD<<"(): Cannot determine number of flavours."<<std::endl;
    }
    if (nlight>6) {
      msg_Error()<<METHOD<<"(): Too many light flavours: "<<nlight
                         <<",  maximum is 6"<<std::endl;
    }
    double default_alphaQCD=Recola_Interface::GetDefaultAlphaQCD();
    double default_scale=Recola_Interface::GetDefaultScale();
    set_alphas_rcl(default_alphaQCD,sqrt(default_scale),nlight);
    msg_Debugging()<<"use \\alpha_s="<<alphas<<" at \\mu_R="<<sqrt(muR2)
                   <<std::endl;

    msg_Out()<<"Processes in Recola are being generated..." << endl;
    generate_processes_rcl();
    Recola_Interface::setProcGenerationTrue();
    msg_Out()<<"Process generation in Recola completed..." << endl;
  }

  void Recola_Interface::EvaluateLoop(int id, const Vec4D_Vector& momenta, double& bornres, 
                                      METOOLS::DivArrD& virt, std::vector<double> &asscontribs)
  {
    vector<double> pp(4*momenta.size());

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
    
    PHASIC::Process_Info pi(s_procmap[id]);
    
    /*
    if (s_interference[id]){
      get_squared_amplitude_rcl(id,pi.m_maxcpl[0]-pi.m_fi.m_nlocpl[0],"LO",fA2[0]);
      get_squared_amplitude_rcl(id,pi.m_maxcpl[0],"NLO",fA2[1]);
    }
    */

    double boqcd = pi.m_maxcpl[0]-pi.m_fi.m_nlocpl[0];
    double voqcd = pi.m_maxcpl[0];
 

    msg_Debugging()<<"Getting Born ..."<<std::endl;
    get_squared_amplitude_rcl(id,boqcd,"LO",fA2[0]);
    msg_Debugging()<<"... B="<<fA2[0]<<std::endl;
    msg_Debugging()<<"Getting V ..."<<std::endl;
    get_squared_amplitude_rcl(id,voqcd,"NLO",fA2[1]);
    msg_Debugging()<<"... V="<<fA2[1]<<std::endl;
    bornres = fA2[0];
    virt.Finite()=fA2[1];


   if (s_asscontribs[id]) {
      if (s_asscontribs[id]&asscontrib::EW) {
        if (!(asscontribs.size()>0)) THROW(fatal_error,"Inconsistent state.");
        msg_Debugging()<<"Getting V_EW ..."<<std::endl;
        get_squared_amplitude_rcl(id,boqcd,"NLO",asscontribs[0]);
        msg_Debugging()<<"... V_EW="<<asscontribs[0]<<std::endl;
      }
      if (s_asscontribs[id]&asscontrib::LO1) {
        if (!(asscontribs.size()>1)) THROW(fatal_error,"Inconsistent state.");
        msg_Debugging()<<"Getting BsubLO1 ..."<<std::endl;
        if (boqcd>=1) get_squared_amplitude_rcl(id,boqcd-1,"LO",asscontribs[1]);
        msg_Debugging()<<"... BsubLO1="<<asscontribs[1]<<std::endl;
      }
      if (s_asscontribs[id]&asscontrib::LO2) {
        if (!(asscontribs.size()>2)) THROW(fatal_error,"Inconsistent state.");
        msg_Debugging()<<"Getting BsubLO2 ..."<<std::endl;
        if (boqcd>=2) get_squared_amplitude_rcl(id,boqcd-2,"LO",asscontribs[2]);
        msg_Debugging()<<"... BsubLO2="<<asscontribs[2]<<std::endl;
      }
      if (s_asscontribs[id]&asscontrib::LO3) {
        if (!(asscontribs.size()>3)) THROW(fatal_error,"Inconsistent state.");
        msg_Debugging()<<"Getting BsubLO3 ..."<<std::endl;
        if (boqcd>=3) get_squared_amplitude_rcl(id,boqcd-3,"LO",asscontribs[3]);
        msg_Debugging()<<"... BsubLO3="<<asscontribs[3]<<std::endl;
      }
    }

    if (s_compute_poles) {
      double fE12[2]={0.};
        double fE22[2]={0.};
      // compute 1/eps^2
      set_delta_ir_rcl(0.,100.+M_PI*M_PI/6.);
      compute_process_rcl(id,fpp,"NLO",fE22);
      msg_Info()<<"B = "<<fE22[0]<<" <==> "<<fA2[0]<<std::endl;
      // compute 1/eps
      set_delta_uv_rcl(100.);
      compute_process_rcl(id,fpp,"NLO",fE12);
      msg_Info()<<"B = "<<fE12[0]<<" <==> "<<fA2[0]<<std::endl;
      // reset
      set_delta_ir_rcl(0.,M_PI*M_PI/6.);
      set_delta_uv_rcl(0.);
    }
  }
  

  void Recola_Interface::EvaluateReal(int id, const Vec4D_Vector& momenta, METOOLS::DivArrD& real, int amptype)
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
    int procIndex(id);
    PHASIC::Process_Info pi(s_procmap[id]);
    
    /*if (s_interference[procIndex]){
      get_squared_amplitude_rcl(id,pi.m_maxcpl[0],"LO",fA2[0]);
      }*/
    if(amptype==12)
    {
      compute_process_rcl(id,fpp,"NLO",fA2,momcheck);
      real.Finite() = fA2[1];

    }
    else if (amptype==1)
    {
      compute_process_rcl(id,fpp,"LO",fA2,momcheck);
      real.Finite() = fA2[0];
    } 
  }

  void Recola_Interface::EvaluateBorn(int id, const Vec4D_Vector& momenta, double& bornres, int amptype)
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
    int procIndex(id);
    PHASIC::Process_Info pi(s_procmap[id]);
    
    /*if (s_interference[procIndex]){
      get_squared_amplitude_rcl(id,pi.m_maxcpl[0],"LO",fA2[0]);
      }*/
    if(amptype==12)
    {
      compute_process_rcl(id,fpp,"NLO",fA2,momcheck);
      bornres = fA2[1];
    }
    else if (amptype==1)
    {
      compute_process_rcl(id,fpp,"LO",fA2,momcheck);
      bornres = fA2[0];
    } 
  }
  
  size_t Recola_Interface::PDFnf(double scale, size_t maxn){
    size_t nf(0);
    double sqrtscale(sqrt(scale));
    for (size_t i(0);i<=maxn;++i) {
      nf=i;
      if (sqrtscale<s_pdfmass[i]) break;
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
