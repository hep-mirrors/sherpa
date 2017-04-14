#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include "MODEL/Main/Model_Base.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Process/Process_Base.H"
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
  bool           Recola_Interface::s_ignore_model = false;
  bool           Recola_Interface::s_exit_on_error= true;
  double         Recola_Interface::s_light_fermion_threshold=0.1;
  unsigned int   Recola_Interface::s_recolaProcIndex = 0;
  bool           Recola_Interface::s_processesGenerated = false;
  double         Recola_Interface::s_default_alpha = -1.;
  double         Recola_Interface::s_default_scale = -1.;
  int            Recola_Interface::s_default_flav = 0.;
  int            Recola_Interface::s_getPDF_default = 0;
  int            Recola_Interface::s_fixed = 0;
  std::vector<double> Recola_Interface::s_pdfmass(6);
  
  std::map<int,PHASIC::Process_Info> m_procmap;

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

    std::cout << "ERROR: particle id not found in particle2Recola" << std::endl;
    std::cout << "WARNING: Goldstone modes 'p0', 'p+', 'p-' not implemented" << std::endl;
    std::cout << "consider particle list in ATOOLS/Phys/Flavour_Tags.H" << std::endl;
    exit(1);
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

    std::cout << "ERROR: particle id="<<p<< " not found in particle2Recola" << std::endl;
    std::cout << "WARNING: Goldstone modes 'p0', 'p+', 'p-' not implemented" << std::endl;
    std::cout << "consider particle list in ATOOLS/Phys/Flavour_Tags.H" << std::endl;
    exit(1);
  }


  std::string Recola_Interface::process2Recola(const Process_Info &pi){
    Flavour_Vector fl=pi.ExtractFlavours();
    std::string process = particle2Recola(fl[0].IDName())
      + " " + particle2Recola(fl[1].IDName()) + " -> ";
    for(size_t i=2; i<fl.size(); ++i){
      process += particle2Recola(fl[i].IDName())+" ";
    }
    return process;
  }

  Recola_Interface::~Recola_Interface()
  {

  }

  bool Recola_Interface::Initialize(const string &path,const string &file,
				    MODEL::Model_Base *const model,
				    BEAM::Beam_Spectra_Handler *const beam,
				    PDF::ISR_Handler *const isr)
  {
    // find RECOLA installation prefix with several overwrite options
    s_default_alpha=s_default_scale=-1.;
    s_default_flav=0;
    struct stat st;
    Data_Reader reader(" ",";","#","=");
    s_ignore_model = reader.GetValue<int>("RECOLA_IGNORE_MODEL",0); 

    s_exit_on_error = reader.GetValue<int>("RECOLA_EXIT_ON_ERROR",1); 
    if (s_ignore_model) msg_Info()<<METHOD<<"(): Recola will use the "
                                  <<"Standard Model even if you set a "
                                  <<"different model without warning."
                                  <<std::endl;

    s_recolaprefix = rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Recola";
    s_getPDF_default = reader.GetValue<int>("RECOLA_GETPDF_DEFAULT",0);

    if(stat(s_recolaprefix.c_str(),&st) != 0) s_recolaprefix = RECOLA_PREFIX;
    s_recolaprefix = reader.GetValue<string>("RECOLA_PREFIX", s_recolaprefix);
    msg_Info()<<"Initialising Recola generator from "<<s_recolaprefix<<endl;

    if(MODEL::s_model->Name() != "SM"){
      THROW(fatal_error, "ONLY Standard Model so far supported in RECOLA");
    }

    // load library dynamically
    s_loader->AddPath(s_recolaprefix);
    if (!s_loader->LoadLibrary("recola")) THROW(fatal_error, "Failed to load librecola.");

    int recolaVerbosity=0;
    recolaVerbosity = reader.GetValue<int>("RECOLA_VERBOSITY",recolaVerbosity);
    if(recolaVerbosity<0 || recolaVerbosity >2){
      cout << "no valid Value for RECOLA_VERBOSITY"<< endl;
      cout << "Verbosity set to 'silent'" << endl;
      recolaVerbosity = 0;
    }
    set_print_level_squared_amplitude_rcl(recolaVerbosity);
    set_print_level_amplitude_rcl(recolaVerbosity);
    set_print_level_correlations_rcl(recolaVerbosity);

    string recolaOutput = reader.GetValue<string>("RECOLA_OUTPUT","*");
    set_output_file_rcl(recolaOutput);

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
    
    PDF::PDF_Base *pdf(NULL);
    pdf=isr->PDF(0);
    s_default_alpha=pdf->GetDefaultAlpha();
    s_default_scale=pdf->GetDefaultScale();
    int pdfnf=pdf->GetFlavourScheme();

    if (pdfnf>10) pdfnf-=10;
    if (pdfnf==-1) pdfnf=6;

    double cmass=pdfnf>3?pdf->ASInfo().m_flavs[3].m_thres:Flavour(kf_c).Mass();
    double bmass=pdfnf>4?pdf->ASInfo().m_flavs[4].m_thres:Flavour(kf_b).Mass();
    double tmass=pdfnf>5?pdf->ASInfo().m_flavs[5].m_thres:Flavour(kf_t).Mass();

    cmass=reader.GetValue<double>("RECOLA_AS_RUN_MASS_C", cmass);
    bmass=reader.GetValue<double>("RECOLA_AS_RUN_MASS_B", bmass);
    tmass=reader.GetValue<double>("RECOLA_AS_RUN_MASS_T", tmass);

    for (int i=0; i<3; i++){
      if (i<pdfnf)
	s_pdfmass[i]=pdf->ASInfo().m_flavs[i].m_thres;
    }
    s_pdfmass[3]=cmass;
    s_pdfmass[4]=bmass;
    s_pdfmass[5]=tmass;

    set_alphas_masses_rcl(cmass,bmass,tmass,Flavour(kf_c).Width(),Flavour(kf_b).Width(),Flavour(kf_t).Width()); 
    return true;
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
    define_process_rcl(procIndex,process2Recola(pi),"NLO");
    unselect_all_gs_powers_BornAmpl_rcl(procIndex);
    unselect_all_gs_powers_LoopAmpl_rcl(procIndex);
    Data_Reader reader(" ",";","#","=");
    reader.AddIgnore("[");
    reader.AddIgnore("]");
    if(pi.m_fi.m_nloqcdtype==nlo_type::loop){
      select_gs_power_BornAmpl_rcl(procIndex,pi.m_maxcpl[0]-(pi.m_fi.m_nloqcdtype==nlo_type::loop));
      select_gs_power_LoopAmpl_rcl(procIndex,pi.m_maxcpl[0]+(pi.m_fi.m_nloqcdtype==nlo_type::loop));
    }
    else if (amptype==12) {
      select_gs_power_LoopAmpl_rcl(procIndex,pi.m_maxcpl[0]);
    }
    msg_Debugging()<<"procIndex "<<procIndex<<" returned\n";  
    return procIndex;
  }
  
  void Recola_Interface::EvaluateLoop(int id, const Vec4D_Vector& momenta, double& bornres, METOOLS::DivArrD& virt)
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
    
//    compute_process_rcl(id,fpp,NN,"NLO",fA2);
    compute_process_rcl(id,fpp,"NLO",fA2); // Change discussed in meeting. Mathieu 12/04/2017
    int procIndex(id);
    PHASIC::Process_Info pi(m_procmap[id]);
    
    get_squared_amplitude_rcl(id,pi.m_maxcpl[0]-(pi.m_fi.m_nloqcdtype==nlo_type::loop),"LO",fA2[0]);
    get_squared_amplitude_rcl(id,pi.m_maxcpl[0],"NLO",fA2[1]);
 
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
