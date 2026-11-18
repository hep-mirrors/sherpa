#include "AddOns/GoSam/GoSam_Interface.H"
#include "AddOns/GoSam/GoSam_Communicator.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Library_Loader.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Phys/KF_Table.H"

#include <fstream>

using namespace GoSam;
using namespace PHASIC;
using namespace ATOOLS;
using namespace std;


string GoSam_Interface::s_contractName = "";
string GoSam_Interface::s_type = "";
bool GoSam_Interface::s_newLibraries = false;
bool GoSam_Interface::s_contract = false;

GoSam_Interface::GoSam_Interface() : PHASIC::ME_Generator_Base("GoSam"),
  m_drmode(-1), m_widthscheme(""),
  m_contractBaseName("") {
  RegisterDefaults();
}

void GoSam_Interface::RegisterDefaults() {
    auto s = Settings::GetMainSettings()["GOSAM"];
    s["IRregularisation"].SetDefault("CDR");
    s["ContractFileName"].SetDefault("SherpaGoSam");
    auto pdata = s["ParticleCodeMap"];
}

bool GoSam_Interface::Initialize(MODEL::Model_Base *const model,
                                 BEAM::Beam_Spectra_Handler *const beam,
                                 PDF::ISR_Handler *const isr,
                                 YFS::YFS_Handler *const yfs) {
  auto s = Settings::GetMainSettings()["GOSAM"];
  m_contractBaseName = s["ContractFileName"].Get<string>();
  string irr = s["IRregularisation"].Get<string>();
  if(irr == "CDR") m_drmode = 0;
  else if(irr == "DRED") m_drmode = 1;

  if(m_drmode <  0) THROW(fatal_error,"Unknown or unset IRregularisation.");
  /* if(m_drmode == 1) THROW(fatal_error,"Currently only CDR supported as IRregularisation."); */

  if (model->ScalarNumber(std::string("WidthScheme")))
    m_widthscheme = "ComplexMass";
  else
    m_widthscheme = "FixedWidth";
  if(m_widthscheme == "") THROW(fatal_error,"Unknown or unset WidthScheme.");

  s_contract = false;
  s_contractName = rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/GoSam/" + m_contractBaseName;
  string cname = s_contractName + ".olc";
  if (access( cname.c_str(), F_OK ) != -1)  {
    s_contract = true;
    s_contractName = cname;
  }
  else {
    s_contractName += ".olp";
  }
  s_type = s["BornType"].SetDefault("Tree").Get<string>();
  if(s_type != "Tree" && s_type != "LoopInduced") THROW(fatal_error,"Unknown Born type.");
  if (!s_contract) {
    MakeDir(rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/GoSam");
    GoSam_Communicator lhfile(s_contractName);
    if (lhfile.FileStatus()==0) {
      lhfile.AddParameter("InterfaceVersion    BLHA2");
      lhfile.AddParameter("CorrectionType      QCD");
      if (m_drmode==1) lhfile.AddParameter("IRregularisation    DRED");
      else             lhfile.AddParameter("IRregularisation    CDR");
      lhfile.AddParameter("WidthScheme         "+m_widthscheme);
      if ( s["SetEWScheme"].SetDefault(true).Get<bool>() ) {
        string ewscheme = Settings::GetMainSettings()["EW_SCHEME"].Get<string>();
        lhfile.AddParameter("EWScheme            "+ewscheme);
      }
      /* Should be quarks only here?*/
      int tmparr[] = { kf_c, kf_b, kf_t, kf_Wplus, kf_Z, kf_h0};
      vector<int> pdgids (tmparr, tmparr + sizeof(tmparr) / sizeof(tmparr[0]) );
      std::string massive_string="";
      for (size_t i=0; i<pdgids.size(); ++i)  if (Flavour(pdgids[i]).Mass()>0.0)   massive_string+=(" "+ToString(pdgids[i]));
      if (massive_string.length()) lhfile.AddParameter("MassiveParticles "+massive_string);

      lhfile.AddParameter("");
      lhfile.AddParameter("# process list");
    }
    string gosamSettings = "";
    for(auto& key: s["Settings"].GetKeys()) {
      string value = s["Settings"][key].SetDefault("").Get<string>();
      if(value != "") {
        gosamSettings += key+"="+value+"\n";
      }
    }
    ofstream ofile;
    string gosamSettingsName = rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/GoSam/"+"gosam.in";
    ofile.open(gosamSettingsName.c_str(),ios::trunc);
    ofile<<gosamSettings;
    ofile.close();
    msg_Out()<<"Wrote parameters from runcard to "<<gosamSettingsName<<":\n"<<gosamSettings;
  }
  else {
    s_loader->LoadLibrary(rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/GoSam/build","golem_olp");
    int check = 0;
    int errmode = 1;


    bool autosettings = s["AutoSetParameters"].SetDefault(true).Get<bool>();
    if(autosettings) {
      auto pdata = s["ParticleCodeMap"];
      for(kf_code kf: {kf_Z,kf_Wplus,kf_h0,kf_e,kf_mu,kf_tau,kf_s,kf_c,kf_b,kf_t}) {
        SetMass(kf,pdata[ToString(kf)].SetDefault("").Get<string>());
      }
      for(kf_code kf: {kf_Z,kf_Wplus}) {
        SetWidth(kf,pdata[ToString(kf)].SetDefault("").Get<string>());
      }
      SetParameter("alpha",MODEL::s_model->ScalarConstant("alpha_QED"),errmode);
    }
    for(auto& key: s["Parameters"].GetKeys()) {
      double value = s["Parameters"][key].SetDefault(0).Get<double>();
      SetParameter(key,value);
    }
    
    /// Start
    GoSam::OLP_Start(s_contractName.c_str(),&check);
    string reportname = rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/GoSam/settingsreport.txt";
    Remove(reportname);
    OLP_PrintParameter(reportname.c_str());
  }
  return true;
}

int GoSam_Interface::RegisterProcess(const Process_Info& pi,const Flavour_Vector& flavs,
                                     bool active) {
  if (!active) return -1;
  size_t pn = flavs.size();
  int proc_id = -1;

  GoSam_Communicator lhfile(s_contractName);
  if (!s_contract) {
    if(lhfile.CheckProcess(2,pn-2,flavs)==-1) {
      lhfile.AddParameter("AmplitudeType            Loop");
      lhfile.AddParameter("CouplingPower QCD        "+ToString(pi.m_maxcpl[0]-1));
      lhfile.AddParameter("CouplingPower QED        "+ToString(pi.m_maxcpl[1]));
      lhfile.AddProcess(2,pn-2,flavs);
    }
    s_newLibraries = true;
    Copy(rpa->gen.Variable("SHERPA_SHARE_PATH")+"/makelibs",
         rpa->gen.Variable("SHERPA_CPP_PATH")+"/makelibs");
    Copy(rpa->gen.Variable("SHERPA_SHARE_PATH")+"/makelibs-gosam",
         rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/GoSam/makelibs");

    return proc_id;
  }
  else {
    if (lhfile.CheckParameterStatus()!=1) {
      THROW(fatal_error,"Bad parameter.");
    }
    int pstatus=lhfile.CheckProcess(2,pn-2,flavs);
    std::string pstr("");
    switch (pstatus) {
    case -2:
    case 0:
      for (size_t i(0);i<flavs.size();++i) pstr+=ToString((long int)flavs[i])+" ";
      THROW(fatal_error,"Process "+pstr+"not found in contract file.");
    case -1:
      THROW(fatal_error,"Process not found in contract file.");
    default:
      if (pstatus!=1) msg_Info()<<"Found "<<pstatus<<" subprocesses. Cannot "
                                <<"handle this yet, only first ID will be used!"<<endl;
      proc_id=lhfile.GetID(2,pn-2,flavs,0);
    }
    return proc_id;
  }
}

int GoSam_Interface::RegisterProcess(const PHASIC::External_ME_Args& args,bool active) {
  if (!active) return -1;
  size_t pn = args.m_inflavs.size()+args.m_outflavs.size();
  Flavour_Vector flavs = args.m_inflavs;
  for(Flavour fl: args.m_outflavs) flavs.push_back(fl);
  int proc_id = -1;
  GoSam_Communicator lhfile(s_contractName);
  if (!s_contract) {
    if(lhfile.CheckProcess(args.m_inflavs.size(),args.m_outflavs.size(),flavs)==-1) {

      lhfile.AddParameter("AmplitudeType         "+s_type);
      lhfile.AddParameter("CouplingPower QCD        "+ToString(args.m_orders[0]));
      lhfile.AddParameter("CouplingPower QED        "+ToString(args.m_orders[1]));
      lhfile.AddProcess(args.m_inflavs.size(),args.m_outflavs.size(),flavs);
    }
    s_newLibraries = true;
    Copy(rpa->gen.Variable("SHERPA_SHARE_PATH")+"/makelibs",
         rpa->gen.Variable("SHERPA_CPP_PATH")+"/makelibs");
    Copy(rpa->gen.Variable("SHERPA_SHARE_PATH")+"/makelibs-gosam",
         rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/GoSam/makelibs");
    return proc_id;
  }
  else {
    if (lhfile.CheckParameterStatus()!=1) {
      THROW(fatal_error,"Bad parameter.");
    }
    int pstatus=lhfile.CheckProcess(args.m_inflavs.size(),args.m_outflavs.size(),flavs);
    std::string pstr("");
    switch (pstatus) {
    case -2:
    case 0:
      for (size_t i(0);i<flavs.size();++i) pstr+=ToString((long int)flavs[i])+" ";
      THROW(fatal_error,"Process "+pstr+"not found in contract file.");
    case -1:
      THROW(fatal_error,"Process not found in contract file.");
    default:
      if (pstatus!=1) msg_Info()<<"Found "<<pstatus<<" subprocesses. Cannot "
                                <<"handle this yet, only first ID will be used!"<<endl;
      proc_id=lhfile.GetID(args.m_inflavs.size(),args.m_outflavs.size(),flavs,0);
    }
    return proc_id;
  }
}

int  GoSam_Interface::PerformTests() { return 1; }
void GoSam_Interface::PrepareTerminate() {}


DECLARE_GETTER(GoSam_Interface,"GoSam",ME_Generator_Base,ME_Generator_Key);

ME_Generator_Base *ATOOLS::Getter<ME_Generator_Base,ME_Generator_Key,
                                  GoSam_Interface>::
operator()(const ME_Generator_Key &key) const
{
  return new GoSam_Interface();
}

void ATOOLS::Getter<ME_Generator_Base,ME_Generator_Key,GoSam_Interface>::
PrintInfo(ostream &str,const size_t width) const
{
  str<<"Interface to the GoSam loop ME generator";
}
