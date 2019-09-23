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
using namespace Recola;
using namespace PHASIC;
using namespace MODEL;
using namespace ATOOLS;
using namespace std;

std::string    Recola_Interface::s_recolaprefix = std::string("");
bool           Recola_Interface::s_ignore_model = false;
bool           Recola_Interface::s_exit_on_error = true;
bool           Recola_Interface::s_use_iop_in_ewapprox = false;
double         Recola_Interface::s_light_fermion_threshold = 0.1;
size_t         Recola_Interface::s_recolaProcIndex = 0;
bool           Recola_Interface::s_processesGenerated = false;
double         Recola_Interface::s_default_alphaqcd = -1.;
double         Recola_Interface::s_default_scale = -1.;
int            Recola_Interface::s_default_flav = 0.;
int            Recola_Interface::s_getPDF_default = 0;
int            Recola_Interface::s_fixed_flav = 0;
double         Recola_Interface::s_ir_scale = 100.;
double         Recola_Interface::s_uv_scale = 100.;
int            Recola_Interface::s_collier_cache = -1;
size_t         Recola_Interface::s_doint = 1;
std::vector<double> Recola_Interface::s_pdfmass(6);

std::map<size_t,PHASIC::Process_Info> Recola_Interface::s_procmap;
std::map<size_t,PHASIC::asscontrib::type> Recola_Interface::s_asscontribs;

std::ostream & Recola::operator<<(std::ostream & s,
                                  const amptype::code &at)
{
  if      (at==amptype::none)      s<<"none";
  else if (at==amptype::treetree)  s<<"|tree|^2";
  else if (at==amptype::treeloop)  s<<"2Re(tree*loop)";
  else if (at==amptype::looploop)  s<<"|loop|^2";
  return s;
}

std::string Recola_Interface::particle2Recola(const int p)
{
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

  msg_Error()<<"ERROR: particle id="<<p<<" not found in particle2Recola"<<std::endl;
  msg_Error()<<"WARNING: Goldstone modes 'p0', 'p+', 'p-' not implemented"<<std::endl;
  msg_Error()<<"consider particle list in ATOOLS/Phys/Flavour_Tags.H"<<std::endl;
  THROW(fatal_error,"Particle ID could not be translated.");
}

std::string Recola_Interface::particle2Recola(const std::string p)
{
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

  msg_Error()<<"ERROR: particle id="<<p<<" not found in particle2Recola"<<std::endl;
  msg_Error()<<"WARNING: Goldstone modes 'p0', 'p+', 'p-' not implemented"<<std::endl;
  msg_Error()<<"consider particle list in ATOOLS/Phys/Flavour_Tags.H"<<std::endl;
  THROW(fatal_error,"Particle ID could not be translated.");
}


std::string Recola_Interface::process2Recola(const Process_Info &pi)
{
  Flavour_Vector fl=pi.ExtractFlavours();
  std::string process = particle2Recola(fl[0].IDName())
    + " " + particle2Recola(fl[1].IDName()) + " -> ";
  for(size_t i=2; i<fl.size(); ++i){
    process += particle2Recola(fl[i].IDName())+" ";
  }
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
  s_ignore_model = reader.GetValue<int>("RECOLA_IGNORE_MODEL",0);

  s_exit_on_error = reader.GetValue<int>("RECOLA_EXIT_ON_ERROR",1);
  if (s_ignore_model) msg_Info()<<METHOD<<"(): Recola will use the "
                                <<"Standard Model even if you set a "
                                <<"different model without warning."
                                <<std::endl;

  s_exit_on_error = reader.GetValue<size_t>("RECOLA_EXIT_ON_ERROR",1);
  s_use_iop_in_ewapprox = reader.GetValue<size_t>("RECOLA_USE_I_IN_EWAPPROX",0);

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
  if (!s_loader->LoadLibrary("recola"))
    THROW(fatal_error, "Failed to load librecola.");

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

  // set fixed IR and UV scales
  s_ir_scale=reader.GetValue<int>("RECOLA_IR_SCALE",100.);
  s_uv_scale=reader.GetValue<int>("RECOLA_UV_SCALE",100.);

  // set particle masses/widths
  // boson masses
  int recolaOnShellZW = reader.GetValue<int>("RECOLA_ONSHELLZW",0);
  if(recolaOnShellZW != 0) {
    set_onshell_mass_z_rcl(Flavour(kf_Z).Mass(),Flavour(kf_Z).Width());
    set_onshell_mass_w_rcl(Flavour(kf_Wplus).Mass(),Flavour(kf_Wplus).Width());
  }
  else {
    set_pole_mass_z_rcl(Flavour(kf_Z).Mass(),Flavour(kf_Z).Width());
    set_pole_mass_w_rcl(Flavour(kf_Wplus).Mass(),Flavour(kf_Wplus).Width());
  }
  set_pole_mass_h_rcl(Flavour(kf_h0).Mass(),Flavour(kf_h0).Width());
  // lepton masses
  set_pole_mass_electron_rcl(Flavour(kf_e).Mass());
  set_pole_mass_muon_rcl(Flavour(kf_mu).Mass(),Flavour(kf_mu).Width());
  set_pole_mass_tau_rcl(Flavour(kf_tau).Mass(),Flavour(kf_tau).Width());
  // quark masses
  set_pole_mass_up_rcl(Flavour(kf_u).Mass());
  set_pole_mass_down_rcl(Flavour(kf_d).Mass());
  set_pole_mass_strange_rcl(Flavour(kf_s).Mass());
  set_pole_mass_charm_rcl(Flavour(kf_c).Mass(),Flavour(kf_c).Width());
  set_pole_mass_bottom_rcl(Flavour(kf_b).Mass(),Flavour(kf_b).Width());
  set_pole_mass_top_rcl(Flavour(kf_t).Mass(),Flavour(kf_t).Width());
  // light fermion threshold
  s_light_fermion_threshold=reader.GetValue<double>("RECOLA_LIGHT_FERMION_THRESHOLD",0);
  set_light_fermions_rcl(s_light_fermion_threshold);

  // adapt the conventions from COLLIER to Catani-Seymour
  set_delta_ir_rcl(0.0,M_PI*M_PI/6.0);

  // pdf consistency
  PDF::PDF_Base *pdf(isr->PDF(0));
  s_default_alphaqcd=pdf->ASInfo().m_asmz;
  s_default_scale=pdf->ASInfo().m_mz2;
  int pdfnf=pdf->ASInfo().m_flavs.size();
  s_default_flav=pdfnf;
  s_fixed_flav=reader.GetValue<int>("RECOLA_FIXED_FLAVS",s_default_flav+10);
  if (pdfnf>10) pdfnf-=10;
  if (pdfnf==-1) pdfnf=6;

  double cmass(0), bmass(0), tmass(0);
  if (pdf->ASInfo().m_allflavs.size()==0) {
    cmass=Flavour(kf_c).Mass();
    bmass=Flavour(kf_b).Mass();
    tmass=Flavour(kf_t).Mass();
  }
  else {
    cmass=pdf->ASInfo().m_allflavs[3].m_mass;
    bmass=pdf->ASInfo().m_allflavs[4].m_mass;
    tmass=pdf->ASInfo().m_allflavs[5].m_mass;
  }
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
  set_alphas_masses_rcl(cmass,bmass,tmass,
                        Flavour(kf_c).Width(),Flavour(kf_b).Width(),
                        Flavour(kf_t).Width());
  s_collier_cache=reader.GetValue<int>("RECOLA_COLLIER_CACHE",-1);
  s_doint=reader.GetValue<size_t>("RECOLA_INTERFERENCE",1);
  return true;
}

int Recola_Interface::RegisterProcess(const Process_Info& pi,
                                      const amptype::code& at)
{
  // check whether valid
  if (pi.m_maxcpl[0]<0 || pi.m_maxcpl[1]<0) {
    msg_Error()<<"Couplings are "<<pi.m_maxcpl<<", not initialising process"
                   <<std::endl;
    return -1;
  }
  DEBUG_FUNC(pi.m_maxcpl);
  // init process
  increaseProcIndex();
  msg_Debugging()<<"Recola_Interface::RegisterProcess called\n";
  int procIndex(getProcIndex());
  msg_Debugging()<<"ProcIndex = " <<procIndex <<"\n";
  msg_Debugging()<<"process string = "<<process2Recola(pi)<<"\n";
  // set procIndex to map with flavours
  s_procmap[procIndex]=pi;
  if (!pi.m_nlomode && at!=amptype::looploop) {
    msg_Error() << "no NLO mode detected!\n";
    return 0;
  }
  // define process in Recola, at this stage always 'NLO'
  define_process_rcl(procIndex,process2Recola(pi),"NLO");
  // set collier cache
  if (s_collier_cache>=0) split_collier_cache_rcl(procIndex,s_collier_cache);
  // find out whether we need multiple orders or not
  s_asscontribs[procIndex]=pi.m_fi.m_asscontribs;
  // if we only need specific orders, set them
  if (s_asscontribs[procIndex]==asscontrib::none) {
    // unset all powers of the amplitude
    unselect_all_gs_powers_BornAmpl_rcl(procIndex);
    unselect_all_gs_powers_LoopAmpl_rcl(procIndex);
    // now set the requested powers for the amplitude
    if (pi.m_fi.m_nloqcdtype==nlo_type::loop) {
      double borngspower(pi.m_maxcpl[0]-(pi.m_fi.m_nloqcdtype==nlo_type::loop));
      double loopgspower(pi.m_maxcpl[0]+(pi.m_fi.m_nloqcdtype==nlo_type::loop));
      msg_Debugging()<<"QCD Tree gs-power = "<<borngspower<<std::endl
                     <<"    Loop gs-power = "<<loopgspower<<std::endl;
      select_gs_power_BornAmpl_rcl(procIndex,borngspower);
      select_gs_power_LoopAmpl_rcl(procIndex,loopgspower);
    }
    if (pi.m_fi.m_nloewtype==nlo_type::loop) {
      msg_Debugging()<<"EW Tree gs-power = "<<pi.m_maxcpl[0]<<std::endl
                     <<"   Loop gs-power = "<<pi.m_maxcpl[0]<<std::endl;
      select_gs_power_BornAmpl_rcl(procIndex,pi.m_maxcpl[0]);
      select_gs_power_LoopAmpl_rcl(procIndex,pi.m_maxcpl[0]);
    }
    else if (at==amptype::treetree) {
      msg_Debugging()<<"Born gs-power = "<<pi.m_maxcpl[0]<<std::endl;
      select_gs_power_LoopAmpl_rcl(procIndex,pi.m_maxcpl[0]);
    }
  }
  else {
    msg_Debugging()<<"Initialise Tree and Loop with all gs-powers"<<std::endl;
  }
  msg_Debugging()<<"procIndex "<<procIndex<<" returned\n";
  //
  return procIndex;
}

void Recola_Interface::EvaluateProcess(int id, const Vec4D_Vector& momenta,
                                       const size_t& voqcd, const size_t& boqcd,
                                       METOOLS::DivArrD& Vqcd, double& B,
                                       std::vector<double> &asscontribs)
{
  DEBUG_FUNC("Voqcd="<<voqcd<<", Boqcd="<<boqcd);
  vector<double> pp(4*momenta.size());

  const int NN = momenta.size();
  double fpp[NN][4];

  for (int i=0; i<NN; i++){
    for (int mu=0; mu<4; mu++){
      fpp[i][mu] = momenta[i][mu];
    }
  }
  double fA2[2]={0.};

  compute_process_rcl(id,fpp,"NLO",fA2);
  msg_Debugging()<<"Getting Born ..."<<std::endl;
  get_squared_amplitude_rcl(id,boqcd,"LO",fA2[0]);
  msg_Debugging()<<"... B="<<fA2[0]<<std::endl;
  msg_Debugging()<<"Getting V ..."<<std::endl;
  get_squared_amplitude_rcl(id,voqcd,"NLO",fA2[1]);
  msg_Debugging()<<"... V="<<fA2[1]<<std::endl;

  B = fA2[0];
  Vqcd.Finite()=fA2[1];

  if (s_asscontribs[id]) {
    if (s_asscontribs[id]&asscontrib::EW) {
      if (!asscontribs.size()>0) THROW(fatal_error,"Inconsistent state.");
      msg_Debugging()<<"Getting V_EW ..."<<std::endl;
      get_squared_amplitude_rcl(id,boqcd,"NLO",asscontribs[0]);
      msg_Debugging()<<"... V_EW="<<asscontribs[0]<<std::endl;
    }
    if (s_asscontribs[id]&asscontrib::LO1) {
      if (!asscontribs.size()>1) THROW(fatal_error,"Inconsistent state.");
      msg_Debugging()<<"Getting BsubLO1 ..."<<std::endl;
      if (boqcd>=1) get_squared_amplitude_rcl(id,boqcd-1,"LO",asscontribs[1]);
      msg_Debugging()<<"... BsubLO1="<<asscontribs[1]<<std::endl;
    }
    if (s_asscontribs[id]&asscontrib::LO2) {
      if (!asscontribs.size()>2) THROW(fatal_error,"Inconsistent state.");
      msg_Debugging()<<"Getting BsubLO2 ..."<<std::endl;
      if (boqcd>=2) get_squared_amplitude_rcl(id,boqcd-2,"LO",asscontribs[2]);
      msg_Debugging()<<"... BsubLO2="<<asscontribs[2]<<std::endl;
    }
    if (s_asscontribs[id]&asscontrib::LO3) {
      if (!asscontribs.size()>3) THROW(fatal_error,"Inconsistent state.");
      msg_Debugging()<<"Getting BsubLO3 ..."<<std::endl;
      if (boqcd>=3) get_squared_amplitude_rcl(id,boqcd-3,"LO",asscontribs[3]);
      msg_Debugging()<<"... BsubLO3="<<asscontribs[3]<<std::endl;
    }
  }
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
