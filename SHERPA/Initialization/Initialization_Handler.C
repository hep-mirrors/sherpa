#include "SHERPA/Initialization/Initialization_Handler.H"

#include "SHERPA/PerturbativePhysics/Hard_Decay_Handler.H"
#include "SHERPA/PerturbativePhysics/Shower_Handler.H"
#include "SHERPA/SoftPhysics/Beam_Remnant_Handler.H"
#include "SHERPA/SoftPhysics/Colour_Reconnection_Handler.H"
#include "SHERPA/SoftPhysics/Hadron_Decay_Handler.H"
#include "SHERPA/SoftPhysics/Hadron_Init.H"
#include "SHERPA/SoftPhysics/Soft_Collision_Handler.H"
#include "SHERPA/PerturbativePhysics/MI_Handler.H"
#include "SHERPA/SoftPhysics/Soft_Photon_Handler.H"
#include "SHERPA/Tools/Event_Reader_Base.H"
#include "SHERPA/Main/Filter.H"
#include "PHASIC++/Scales/Core_Scale_Setter.H"
#include "MODEL/Main/Model_Base.H"
#include "METOOLS/Currents/C_Spinor.H"
#include "PDF/Main/Structure_Function.H"
#include "PDF/Main/Intact.H"
#include "PDF/Main/PDF_Base.H"
#include "REMNANTS/Tools/Remnants_Parameters.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include "ATOOLS/Math/Scaling.H"
#include "ATOOLS/Phys/Spinor.H"
#include "ATOOLS/Phys/Variations.H"
#include "ATOOLS/Phys/Fragmentation_Base.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Math/Variable.H"
#include "ATOOLS/Org/Data_Writer.H"
#include "SHERPA/Single_Events/Hadron_Decays.H"
#include "ATOOLS/Org/Library_Loader.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Selectors/Selector.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/Channels/Channel_Generator.H"
#include "PDF/Main/NLOMC_Base.H"
#include "PDF/Main/Shower_Base.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace SHERPA;
using namespace MODEL;
using namespace BEAM;
using namespace PDF;
using namespace REMNANTS;
using namespace ATOOLS;
using namespace std;

typedef void (*PDF_Init_Function)();
typedef void (*PDF_Exit_Function)();

Initialization_Handler::Initialization_Handler() :
  m_mode(eventtype::StandardPerturbative),
  m_savestatus(false), p_model(NULL), p_beamspectra(NULL),
  p_mehandler(NULL), p_harddecays(NULL),
  p_beamremnants(NULL), p_reconnections(NULL),
  p_fragmentation(NULL), p_hdhandler(NULL),
  p_softphotons(NULL), p_evtreader(NULL),
  p_variations(NULL), p_filter(NULL)
{
  RegisterDefaults();
  Settings& s = Settings::GetMainSettings();

  // configure runtime parameters
  if (s["SAVE_STATUS"].Get<std::string>() != "") {
    std::string savestatus(s["SAVE_STATUS"].Get<std::string>());
    if (savestatus[savestatus.size() - 1] != '/') savestatus += "/";
    rpa->gen.SetVariable("SHERPA_STATUS_PATH",
                         rpa->gen.Variable("SHERPA_RUN_PATH") + "/" + savestatus);
    m_savestatus=true;
  } else {
    rpa->gen.SetVariable("SHERPA_STATUS_PATH", "");
  }

  m_evtform = s["EVENT_INPUT"].Get<std::string>();
  if (m_evtform != "") {
    m_mode=eventtype::EventReader;
    msg_Out()<<"Sherpa will read in events as "<<m_evtform<<endl;
  }

  ATOOLS::s_loader->SetCheck(s["CHECK_LIBLOCK"].Get<int>());

  rpa->Init();
  CheckVersion();
  LoadLibraries();
  ShowParameterSyntax();
  ran->InitExternal();

  rpa->gen.SetSoftSC(s["HADRON_DECAYS"]["Spin_Correlations"].Get<int>());
  rpa->gen.SetHardSC(s["HARD_DECAYS"]["Spin_Correlations"].Get<int>());

  exh->AddTerminatorObject(this);
}

void Initialization_Handler::RegisterDefaults()
{
  Settings& s = Settings::GetMainSettings();
  s["BEAM_REMNANTS"].SetDefault(true);
  s["EVENT_GENERATION_MODE"].SetDefault("PartiallyUnweighted");
  s["EVENT_TYPE"].SetDefault("StandardPerturbative");
  s["SOFT_COLLISIONS"].UseNoneReplacements().SetDefault("None");
  s["BEAM_RESCATTERING"].UseNoneReplacements().SetDefault("None");
  s["MI_HANDLER"].UseNoneReplacements().SetDefault("Amisic");
  s["EVT_FILE_PATH"].SetDefault(".");
  s["ANALYSIS_OUTPUT"].SetDefault("Analysis/");
  s["RESULT_DIRECTORY"].SetDefault("Results");
  s["CHECK_LIBLOCK"].SetDefault(0);
  s["OUTPUT_PRECISION"].SetDefault(12);
  s["FILE_SIZE"].SetDefault(std::numeric_limits<size_t>::max());
  s["WRITE_REFERENCES_FILE"].SetDefault(true);

  s["MODEL"].SetDefault("SM");
  s["FRAGMENTATION"].SetDefault("Ahadic").UseNoneReplacements();
  s["HARD_DECAYS"]["Enabled"].SetDefault(false);
  s["N_COLOR"].SetDefault(3.0);

  std::string frag{ s["FRAGMENTATION"].Get<std::string>() };
  s["HADRON_DECAYS"]["Model"].SetDefault("HADRONS++").UseNoneReplacements();
  s["HADRON_DECAYS"]["Max_Proper_Lifetime"].SetDefault(10.0);

  s["HADRON_DECAYS"]["Spin_Correlations"].SetDefault(0);
  auto hdenabled = s["HARD_DECAYS"]["Enabled"].Get<bool>();
  s["HARD_DECAYS"]["Spin_Correlations"].SetDefault(hdenabled);

  s["EVENT_INPUT"].SetDefault("");
  s["STATUS_PATH"].SetDefault("");
  s["PATH"].SetDefault("");
  s["SAVE_STATUS"].SetDefault("");

  s.DeclareVectorSettingsWithEmptyDefault({
      "EVENT_OUTPUT",
      "ANALYSIS",
      "SHERPA_LDADD",
      "SHERPA_VERSION",
      "PDF_LIBRARY",
      "PDF_SET",
      "PDF_SET_VERSIONS",
      "MPI_PDF_LIBRARY",
      "MPI_PDF_SET",
      "MPI_PDF_SET_VERSIONS",
      "BBR_PDF_LIBRARY",
      "BBR_PDF_SET",
      "BBR_PDF_SET_VERSIONS",
      "VARIATIONS",
      "SCALE_VARIATIONS",
      "PDF_VARIATIONS",
      "QCUT_VARIATIONS",
      "BUNCHES",
      "MASSIVE_PS",
      "MASSLESS_PS"
      });
  s["MC@NLO"].DeclareVectorSettingsWithEmptyDefault({
      "DISALLOW_FLAVOUR",
      });
  s.DeclareMatrixSettingsWithEmptyDefault({
      "ASSOCIATED_CONTRIBUTIONS_VARIATIONS"
      });
  s["SHOWER"].DeclareMatrixSettingsWithEmptyDefault({
      "ENHANCE",
      });
  s["EVENT_OUTPUT"].UseNoneReplacements();
  s["VARIATIONS"].UseNoneReplacements();
  s["SCALE_VARIATIONS"].UseNoneReplacements();
  s["PDF_VARIATIONS"].UseNoneReplacements();
  s["QCUT_VARIATIONS"].UseNoneReplacements().SetSynonyms({"CKKW_VARIATIONS"});
  s["PDF_LIBRARY"].UseNoneReplacements();
  s["MPI_PDF_LIBRARY"].UseNoneReplacements();
  s["BBR_PDF_LIBRARY"].UseNoneReplacements();
  s["ANALYSIS"].UseNoneReplacements();

  s["SHOW_ME_GENERATORS"].SetDefault(0);
  s["SHOW_PS_GENERATORS"].SetDefault(0);
  s["SHOW_NLOMC_GENERATORS"].SetDefault(0);
  s["SHOW_SHOWER_GENERATORS"].SetDefault(0);
  s["SHOW_KFACTOR_SYNTAX"].SetDefault(0);
  s["SHOW_SCALE_SYNTAX"].SetDefault(0);
  s["SHOW_SELECTOR_SYNTAX"].SetDefault(0);
  s["SHOW_MODEL_SYNTAX"].SetDefault(0);
  s["SHOW_FILTER_SYNTAX"].SetDefault(0);
  s["SHOW_ANALYSIS_SYNTAX"].SetDefault(0);
  s["SHOW_VARIABLE_SYNTAX"].SetDefault(0);
  s["SHOW_PDF_SETS"].SetDefault(0);

  s["ISR_E_ORDER"].SetDefault(1);
  s["ISR_E_SCHEME"].SetDefault(2);

  s["KFACTOR"].SetDefault("None").UseNoneReplacements();
  s["SCALES"].SetDefault("METS{MU_F2}{MU_R2}{MU_Q2}");
  s["SCALE_FACTOR"].SetDefault(1.0);
  s["FACTORIZATION_SCALE_FACTOR"].SetDefault(1.0);
  s["RENORMALIZATION_SCALE_FACTOR"].SetDefault(1.0);
  s["RESUMMATION_SCALE_FACTOR"].SetDefault(1.0);
  s["USR_WGT_MODE"].SetDefault(true);

  Scoped_Settings metssettings{ Settings::GetMainSettings()["METS"] };
  metssettings["CLUSTER_MODE"].SetDefault(0);

  s["NNLOqT_FOMODE"].SetDefault(0);

  // m_mtmode != 0 to reweight the whole cross section by full mt-dependent LO
  // gg->H cross section and add mt-dependence higher order correcions of the
  // Wilson coefficient for the ggH coupling
  s["HNNLO_MTOP_MODE"].SetDefault(0);
  // m_kfmode = [001]_2 to enable factorized matching of the Wilson coefficient for ggH coupling;
  //          = [010]_2 to enable individual matching of the Wilson coefficient for ggH coupling;
  //          = [100]_2 to remove delta(pT) part of NNLO K factor for a separate LO parton shower.
  s["HNNLO_KF_MODE"].SetDefault(0);

  // shower settings (shower classes are rarely singletons, so we either
  // register settings here or we prevent SetDefault... to called more than once
  // otherwise
  s["SHOWER_GENERATOR"].SetDefault("CSS").UseNoneReplacements();
  std::string showergen{ s["SHOWER_GENERATOR"].Get<std::string>() };
  s["JET_CRITERION"].SetDefault(showergen);
  s["NLOMC_GENERATOR"].SetDefault(showergen);
  auto pss = s["SHOWER"], nlopss = s["MC@NLO"];
  pss["EVOLUTION_SCHEME"].SetDefault(30+30*100);
  pss["KFACTOR_SCHEME"].SetDefault(1);
  pss["SCALE_SCHEME"].SetDefault(14);
  pss["SCALE_VARIATION_SCHEME"].SetDefault(1);
  // TODO: Should this be set to 3.0 for the new Dire default? See the manual
  // Sherpa section on master for details
  pss["FS_PT2MIN"].SetDefault(1.0);
  pss["IS_PT2MIN"].SetDefault(2.0);
  pss["PT2MIN_GSPLIT_FACTOR"].SetDefault(1.0);
  pss["FS_AS_FAC"].SetDefault(1.0);
  pss["IS_AS_FAC"].SetDefault(0.25);
  pss["PDF_FAC"].SetDefault(1.0);
  pss["SCALE_FACTOR"].SetDefault(1.);
  pss["MASS_THRESHOLD"].SetDefault(0.0);
  pss["FORCED_IS_QUARK_SPLITTING"].SetDefault(true);
  pss["FORCED_SPLITTING_GLUON_SCALING"].SetDefault(3./2.);
  s["VIRTUAL_EVALUATION_FRACTION"].SetDefault(1.0);
  pss["RECO_CHECK"].SetDefault(0);
  pss["MAXEM"].SetDefault(std::numeric_limits<size_t>::max());
  pss["REWEIGHT"].SetDefault(showergen != "None");
  s["OUTPUT_ME_ONLY_VARIATIONS"].SetDefault(showergen != "None");
  pss["MAX_REWEIGHT_FACTOR"].SetDefault(1e3);
  nlopss["REWEIGHT_EM"].SetDefault(1);
  pss["REWEIGHT_SCALE_CUTOFF"].SetDefault(5.0);
  pss["KIN_SCHEME"].SetDefault(1);
  nlopss["KIN_SCHEME"].SetDefault(1);
  pss["OEF"].SetDefault(3.0);
  pss["KMODE"].SetDefault(2);
  pss["RESPECT_Q2"].SetDefault(false);
  pss["CKFMODE"].SetDefault(1);
  pss["PDFCHECK"].SetDefault(1);
  pss["QCD_MODE"].SetDefault(1);
  pss["EW_MODE"].SetDefault(false);
  pss["RECO_DECAYS"].SetDefault(0);
  pss["MAXPART"].SetDefault(std::numeric_limits<int>::max());
  pss["PDF_MIN"].SetDefault(1.0e-4);
  pss["PDF_MIN_X"].SetDefault(1.0e-2);
  pss["WEIGHT_CHECK"].SetDefault(false);
  pss["CMODE"].SetDefault(1);
  pss["NCOL"].SetDefault(3);
  pss["RECALC_FACTOR"].SetDefault(4.0);
  pss["TC_ENHANCE"].SetDefault(1.0);
  pss["COUPLING_SCHEME"].SetDefault(1);
  pss["ME_CORRECTION"].SetDefault(0);
  pss["KERNEL_TYPE"].SetDefault(15);
  nlopss["RECALC_FACTOR"].SetDefault(2.0);
  nlopss["PSMODE"].SetDefault(0);
  nlopss["WEIGHT_CHECK"].SetDefault(0);
  nlopss["MAXEM"].SetDefault(1);
  pss["MI_KFACTOR_SCHEME"].SetDefault(0);
  pss["MI_IS_PT2MIN"].SetDefault(4.0);
  pss["MI_FS_PT2MIN"].SetDefault(1.0);
  pss["MI_PT2MIN_GSPLIT_FACTOR"].SetDefault(1.0);
  pss["MI_IS_AS_FAC"].SetDefault(0.66);
  pss["MI_FS_AS_FAC"].SetDefault(0.66);
  pss["MI_KIN_SCHEME"].SetDefault(1);

  s["COMIX_DEFAULT_GAUGE"].SetDefault(1);

  s["DIPOLES"]["SCHEME"].SetDefault(subscheme::CSS);
  s["DIPOLES"]["KAPPA"].SetDefault(2.0/3.0);

  s["COUPLINGS"].SetDefault("Alpha_QCD 1");

  s["EXTRAXS_CSS_APPROX_ME"].SetDefault(false);

  s["RESPECT_MASSIVE_FLAG"].SetDefault(false);

  s["BBR_MODE"].SetDefault(0);
  s["BBR_VETO"].SetDefault(false);
}

Initialization_Handler::~Initialization_Handler()
{
  if (m_savestatus) {
    msg_Error()<<METHOD<<"(): Status saved to '"
	       <<rpa->gen.Variable("SHERPA_STATUS_PATH")<<"'."<<std::endl;
    MakeDir(rpa->gen.Variable("SHERPA_STATUS_PATH"),493);
    exh->PrepareTerminate();
  }
  if (p_evtreader)     { delete p_evtreader;     p_evtreader     = NULL; }
  if (p_mehandler)     { delete p_mehandler;     p_mehandler     = NULL; }
  if (p_reconnections) { delete p_reconnections; p_reconnections = NULL; }
  if (p_fragmentation) { delete p_fragmentation; p_fragmentation = NULL; }
  if (p_beamremnants)  { delete p_beamremnants;  p_beamremnants  = NULL; }
  if (p_harddecays)    { delete p_harddecays;    p_harddecays    = NULL; }
  if (p_hdhandler)     { delete p_hdhandler;     p_hdhandler     = NULL; }
  if (p_softphotons)   { delete p_softphotons;   p_softphotons   = NULL; }
  if (p_beamspectra)   { delete p_beamspectra;   p_beamspectra   = NULL; }
  if (p_model)         { delete p_model;         p_model         = NULL; }
  if (p_variations)    { delete p_variations;    p_variations    = NULL; }
  if (p_filter)        { delete p_filter;        p_filter        = NULL; }
  if (p_yfshandler)    { delete p_yfshandler;    p_yfshandler    = NULL; }
  while (m_analyses.size()>0) {
    delete m_analyses.back();
    m_analyses.pop_back();
  }
  while (m_outputs.size()>0) {
    delete m_outputs.back();
    m_outputs.pop_back();
  }
  while (m_isrhandlers.size()>0) {
    delete m_isrhandlers.begin()->second;
    m_isrhandlers.erase(m_isrhandlers.begin());
  }
  while (m_mihandlers.size()>0) {
    delete m_mihandlers.begin()->second;
    m_mihandlers.erase(m_mihandlers.begin());
  }
  while (m_remnanthandlers.size()>0) {
    REMNANTS::Remnant_Handler * rh = m_remnanthandlers.begin()->second;
    for (REMNANTS::Remnant_Handler_Map::iterator rit=m_remnanthandlers.begin();
	 rit!=m_remnanthandlers.end();) {
      if (rit->second==rh) m_remnanthandlers.erase(rit++);
      else rit++;
    }
    delete rh;
  }
  while (m_schandlers.size()>0) {
    delete m_schandlers.begin()->second;
    m_schandlers.erase(m_schandlers.begin());
  }
  while (m_showerhandlers.size()>0) {
    delete m_showerhandlers.begin()->second;
    m_showerhandlers.erase(m_showerhandlers.begin());
  }
  PHASIC::Phase_Space_Handler::DeleteInfo();
  exh->RemoveTerminatorObject(this);
  for (set<string>::iterator pdflib=m_pdflibs.begin(); pdflib!=m_pdflibs.end();
       ++pdflib) {
    if (*pdflib=="None") continue;
    void *exit(s_loader->GetLibraryFunction(*pdflib,"ExitPDFLib"));
    if (exit==NULL)
      PRINT_INFO("Error: Cannot unload PDF library "+*pdflib);
    else ((PDF_Exit_Function)exit)();
  }
}

void Initialization_Handler::CheckVersion()
{
  std::vector<std::string> versioninfo{
    Settings::GetMainSettings()["SHERPA_VERSION"].GetVector<std::string>() };
  if (versioninfo.empty()) return;
  std::string currentversion(ToString(SHERPA_VERSION)+"."
                                      +ToString(SHERPA_SUBVERSION));
  if (versioninfo.size()==1 && versioninfo[0]!=currentversion) {
    THROW(normal_exit,"Run card request Sherpa "+versioninfo[0]
                      +". This is Sherpa "+currentversion);
  }
  else if (versioninfo.size()==2) {
    if (versioninfo[0]==currentversion || versioninfo[1]==currentversion) return;
    size_t min1(versioninfo[0].find(".",0)),
           min2(versioninfo[0].find(".",min1+1)),
           max1(versioninfo[1].find(".",0)),
           max2(versioninfo[1].find(".",max1+1));
    size_t minmajvers(ToType<size_t>(versioninfo[0].substr(0,min1))),
           minminvers(ToType<size_t>(versioninfo[0].substr(min1+1,min2))),
           minbugvers(ToType<size_t>(versioninfo[0].substr(min2+1))),
           maxmajvers(ToType<size_t>(versioninfo[1].substr(0,max1))),
           maxminvers(ToType<size_t>(versioninfo[1].substr(max1+1,max2))),
           maxbugvers(ToType<size_t>(versioninfo[1].substr(max2+1))),
           curmajvers(ToType<size_t>(currentversion.substr(0,max1))),
           curminvers(ToType<size_t>(currentversion.substr(max1+1,max2))),
           curbugvers(ToType<size_t>(currentversion.substr(max2+1)));
    if (!(CompareVersions(minmajvers,minminvers,minbugvers,
                          curmajvers,curminvers,curbugvers)
          *CompareVersions(curmajvers,curminvers,curbugvers,
                           maxmajvers,maxminvers,maxbugvers)))
      THROW(normal_exit,"Run card request Sherpa "+versioninfo[0]
                        +"-"+versioninfo[1]
                        +". This is Sherpa "+currentversion);
  }
  else THROW(not_implemented,"SHERPA_VERSION information not recognised.");
}

bool Initialization_Handler::CompareVersions
(const size_t& a1,const size_t& b1,const size_t& c1,
 const size_t& a2,const size_t& b2,const size_t& c2)
{
  if (a1<a2) return true;
  if (a1==a2) {
    if (b1<b2) return true;
    if (b1==b2) {
      if (c1<=c2) return true;
    }
  }
  return false;
}

void Initialization_Handler::LoadLibraries()
{
  std::vector<std::string> ldadd =
    Settings::GetMainSettings()["SHERPA_LDADD"].GetVector<std::string>();
  for (size_t i(0);i<ldadd.size();++i) {
    if (!s_loader->LoadLibrary(ldadd[i])) {
      THROW(fatal_error,"Cannot load extra library.");
    }
    else msg_Info()<<METHOD<<"(): Library lib"<<ldadd[i]<<".so loaded.\n";
  }
}

void Initialization_Handler::ShowParameterSyntax()
{
  Settings& s = Settings::GetMainSettings();
  int helpi(s["SHOW_ME_GENERATORS"].Get<int>());
  if (helpi>0) {
    msg->SetLevel(2);
    PHASIC::ME_Generator_Base::ShowSyntax(helpi);
    THROW(normal_exit,"Syntax shown.");
  }
  helpi = s["SHOW_PS_GENERATORS"].Get<int>();
  if (helpi>0) {
    msg->SetLevel(2);
    PHASIC::Channel_Generator::ShowSyntax(helpi);
    THROW(normal_exit,"Syntax shown.");
  }
  helpi = s["SHOW_NLOMC_GENERATORS"].Get<int>();
  if (helpi>0) {
    msg->SetLevel(2);
    PDF::NLOMC_Base::ShowSyntax(helpi);
    THROW(normal_exit,"Syntax shown.");
  }
  helpi = s["SHOW_SHOWER_GENERATORS"].Get<int>();
  if (helpi>0) {
    msg->SetLevel(2);
    PDF::Shower_Base::ShowSyntax(helpi);
    THROW(normal_exit,"Syntax shown.");
  }
  helpi = s["SHOW_KFACTOR_SYNTAX"].Get<int>();
  if (helpi>0) {
    msg->SetLevel(2);
    PHASIC::KFactor_Setter_Base::ShowSyntax(helpi);
    THROW(normal_exit,"Syntax shown.");
  }
  helpi = s["SHOW_SCALE_SYNTAX"].Get<int>();
  if (helpi>0) {
    msg->SetLevel(2);
    if (helpi&1) PHASIC::Scale_Setter_Base::ShowSyntax(helpi);
    if (helpi&2) PHASIC::Core_Scale_Setter::ShowSyntax(helpi);
    THROW(normal_exit,"Syntax shown.");
  }
  helpi = s["SHOW_SELECTOR_SYNTAX"].Get<int>();
  if (helpi>0) {
    msg->SetLevel(2);
    PHASIC::Selector_Base::ShowSyntax(helpi);
    THROW(normal_exit,"Syntax shown.");
  }
  helpi = s["SHOW_MODEL_SYNTAX"].Get<int>();
  if (helpi>0) {
    msg->SetLevel(2);
    MODEL::Model_Base::ShowSyntax(helpi);
    THROW(normal_exit,"Syntax shown.");
  }
  helpi = s["SHOW_FILTER_SYNTAX"].Get<int>();
  if (helpi>0) {
    msg->SetLevel(2);
    Filter::ShowSyntax(helpi);
    THROW(normal_exit,"Syntax shown.");
  }
  helpi = s["SHOW_ANALYSIS_SYNTAX"].Get<int>();
  if (helpi>0) {
    msg->SetLevel(2);
    InitializeTheAnalyses();
    for (Analysis_Vector::iterator it=m_analyses.begin(); it!=m_analyses.end(); ++it)
      (*it)->ShowSyntax(helpi);
    THROW(normal_exit,"Syntax shown.");
  }
  helpi = s["SHOW_VARIABLE_SYNTAX"].Get<int>();
  if (helpi>0) {
    msg->SetLevel(2);
    ATOOLS::Variable_Base<double>::ShowVariables(helpi);
    THROW(normal_exit,"Syntax shown.");
  }
}

std::string StripSectionTags(const std::string &name)
{
  if (name.find('|')!=std::string::npos)
    return name.substr(0,name.find('|'));
  return name;
}

void Initialization_Handler::PrepareTerminate()
{
  Settings& s = Settings::GetMainSettings();
  std::string path(rpa->gen.Variable("SHERPA_STATUS_PATH")+"/");
  if (path=="/") return;
  Copy(s.GetPath(), path + s.GetPath());
  Data_Writer writer;
  writer.SetOutputFile(path+"cmd");
  writer.SetVectorType(vtc::vertical);
  std::vector<std::string> lines = {
    "SHERPA_RUN_PATH = "+rpa->gen.Variable("SHERPA_RUN_PATH"),
    "SHERPA_CPP_PATH = "+rpa->gen.Variable("SHERPA_CPP_PATH"),
    "SHERPA_LIB_PATH = "+rpa->gen.Variable("SHERPA_LIB_PATH"),
  };
  writer.VectorToFile(lines);
}

bool Initialization_Handler::InitializeTheFramework(int nr)
{
  Settings& s = Settings::GetMainSettings();
  bool okay = true;
  const int defgauge{ s["COMIX_DEFAULT_GAUGE"].Get<int>() };
  Spinor<double>::SetDefaultGauge(defgauge);
  Spinor<long double>::SetDefaultGauge(defgauge);
  SetGlobalVariables();
  std::string stag(rpa->gen.Variable("RNG_SEED"));
  while (stag.find(' ')!=std::string::npos) stag.replace(stag.find(' '),1,"-");
  s.AddTag("RNG_SEED", stag);
  Hadron_Init().Init();
  okay = okay && InitializeTheModel();

  if (m_mode==eventtype::StandardPerturbative) {
    std::string eventtype{ s["EVENT_TYPE"].Get<std::string>() };
    if (eventtype=="StandardPerturbative")
      m_mode=eventtype::StandardPerturbative;
    else if (eventtype=="MinimumBias") {
      m_mode=eventtype::MinimumBias;
      if (s["SOFT_COLLISIONS"].Get<string>()==string("Amisic"))
        s["MI_HANDLER"].OverrideScalar<std::string>("Amisic");
      else if (s["SOFT_COLLISIONS"].Get<string>()==string("Shrimps"))
        s["MI_HANDLER"].OverrideScalar<std::string>("None");
      s["ME_GENERATORS"].OverrideScalar<std::string>("None");
    }
    else if (eventtype=="HadronDecay") {
      m_mode=eventtype::HadronDecay;
      s["MI_HANDLER"].OverrideScalar<std::string>("None");
      s["ME_GENERATORS"].OverrideScalar<std::string>("None");
    }
    else {
      THROW(not_implemented,"Unknown event type '"+eventtype+"'");
    }
    msg_Out()<<METHOD<<": mode = "<<m_mode<<", type = "<<eventtype<<": "
      //<<"ME = "<<s["ME_GENERATORS"].Get<string>()<<", "
	     <<"MI = "<<s["SOFT_COLLISIONS"].Get<string>()<<"\n";
  }
  okay = okay && InitializeTheBeams();
  okay = okay && InitializeThePDFs();
  if (!p_model->ModelInit(m_isrhandlers))
    THROW(critical_error,"Model cannot be initialized");
  p_model->InitializeInteractionModel();
  okay = okay && InitializeTheYFS();
  // need to initalize yfs before remnants
  okay = okay && InitializeTheRemnants();
  if (!CheckBeamISRConsistency()) return 0.;
  if (m_mode==eventtype::EventReader) {
    std::string infile;
    size_t bpos(m_evtform.find('[')), epos(m_evtform.rfind(']'));
    if (bpos!=std::string::npos && epos!=std::string::npos) {
      infile=m_evtform.substr(bpos+1,epos-bpos-1);
      m_evtform=m_evtform.substr(0,bpos);
    }
    std::string libname(m_evtform);
    if (libname.find('_')) libname=libname.substr(0,libname.find('_'));
    if (!s_loader->LoadLibrary("Sherpa"+libname+"Input"))
      THROW(missing_module,"Cannot load input library Sherpa"+libname+"Input.");
    p_evtreader = Event_Reader_Base::Getter_Function::GetObject
      (m_evtform,Input_Arguments(s.GetPath(), infile,
				 p_model, m_isrhandlers[isr::hard_process], p_yfshandler));
    if (p_evtreader==NULL) THROW(fatal_error,"Event reader not found");
    msg_Events()<<"SHERPA will read in the events."<<std::endl
  		<<"   The full framework is not needed."<<std::endl;
    InitializeTheAnalyses();
    InitializeTheHardDecays();
    InitializeTheBeamRemnants();
    InitializeTheIO();
    InitializeTheReweighting(Variations_Mode::all);
    return true;
  }
  PHASIC::Phase_Space_Handler::GetInfo();
  okay = okay && InitializeTheShowers();
  okay = okay && InitializeTheHardDecays();
  okay = okay && InitializeTheMatrixElements();
  okay = okay && InitializeTheBeamRemnants();
  if (rpa->gen.NumberOfEvents()>0) {
    okay = okay && InitializeTheUnderlyingEvents();
    okay = okay && InitializeTheSoftCollisions();
    okay = okay && InitializeTheColourReconnections();
    okay = okay && InitializeTheFragmentation();
    okay = okay && InitializeTheHadronDecays();
    okay = okay && InitializeTheSoftPhotons();
    okay = okay && InitializeTheIO();
    okay = okay && InitializeTheFilter();
    okay = okay && InitializeTheReweighting(Variations_Mode::all);
    okay = okay && InitializeTheAnalyses();
  } else {
    okay = okay && InitializeTheReweighting(Variations_Mode::nominal_only);
  }
  return okay;
}

bool Initialization_Handler::CheckBeamISRConsistency()
{
  if (p_model->Name()==std::string("ADD")) {
    double ms = p_model->ScalarConstant("M_s");
    if (ms<rpa->gen.Ecms()) {
      msg_Error()<<"WARNING in Initialization_Handler::CheckBeamISRConsistency :"<<std::endl
	       <<"   You might be using the ADD model beyond its valid range ! "<<endl;
    }
  }

  double smin=0;
  double smax=sqr(rpa->gen.Ecms());
  string name=p_model->Name();
  if (name==std::string("ADD")) {
    double mcut2 = sqr(p_model->ScalarConstant("M_cut"));
    // if ISR & beam -> apply mcut on ISR only
    // if beam only  -> apply mcut on Beam
    smax = Min(smax,mcut2);
    for (size_t i=1;i<3;++i) {
      isr::id id=(isr::id)i;
      if (m_isrhandlers[id]->On()) {
	m_isrhandlers[id]->SetFixedSprimeMax(smax);
	m_isrhandlers[id]->SetFixedSprimeMin(smin);
      }
      else if (p_beamspectra->On()) {
	p_beamspectra->SetSprimeMax(smax);
      }
    }
  }

  if (!(p_beamspectra->CheckConsistency(m_bunch_particles))) {
    msg_Error()<<"Error in Initialization of the Sherpa framework : "<<endl
	       <<"    Detected a mismatch of flavours from beams to bunches : "<<endl
	       <<"    "<<p_beamspectra->GetBeam(0)<<" -> "
	       <<m_isrhandlers[isr::hard_process]->Flav(0)<<" and "
	       <<p_beamspectra->GetBeam(1)<<" -> "
	       <<m_isrhandlers[isr::hard_process]->Flav(1)<<endl;
    return 0;
  }

  return 1;
}

bool Initialization_Handler::InitializeTheIO()
{
  Settings& s = Settings::GetMainSettings();
  auto outputs = s["EVENT_OUTPUT"].GetVector<std::string>();
  std::string outpath=s["EVT_FILE_PATH"].Get<std::string>();
  for (size_t i=0; i<outputs.size(); ++i) {
    if (outputs[i]=="None") continue;
    std::string outfile;
    size_t bpos(outputs[i].find('[')), epos(outputs[i].rfind(']'));
    if (bpos!=std::string::npos && epos!=std::string::npos) {
      outfile=outputs[i].substr(bpos+1,epos-bpos-1);
      outputs[i]=outputs[i].substr(0,bpos);
    }
    std::string libname(outputs[i]);
    if (libname.find('_')) libname=libname.substr(0,libname.find('_'));
    Output_Base* out=Output_Base::Getter_Function::GetObject
      (outputs[i], Output_Arguments(outpath, outfile));
    if (out==NULL) {
      if (!s_loader->LoadLibrary("Sherpa"+libname+"Output"))
	THROW(missing_module,"Cannot load output library Sherpa"+libname+"Output.");
      out=Output_Base::Getter_Function::GetObject
	(outputs[i], Output_Arguments(outpath, outfile));
    }
    if (out==NULL) THROW(fatal_error,"Cannot initialize "+outputs[i]+" output");
    m_outputs.push_back(out);
  }

  return true;
}

bool Initialization_Handler::InitializeTheModel()
{
  Settings& s = Settings::GetMainSettings();
  if (p_model) delete p_model;
  std::string name(s["MODEL"].Get<std::string>());
  p_model=Model_Base::Model_Getter_Function::
    GetObject(name, Model_Arguments(true));
  if (p_model==NULL) {
    if (!s_loader->LoadLibrary("Sherpa"+name))
      THROW(missing_module,"Cannot load model library Sherpa"+name+".");
    p_model=Model_Base::Model_Getter_Function::
      GetObject(name, Model_Arguments(true));
  }
  if (p_model==NULL) THROW(not_implemented,"Model not implemented");
  MODEL::s_model=p_model;
  return 1;
}


bool Initialization_Handler::InitializeTheBeams()
{
  if (p_beamspectra) { delete p_beamspectra; p_beamspectra = NULL; }
  p_beamspectra = new Beam_Spectra_Handler();
  return 1;
}

bool Initialization_Handler::InitializeTheYFS(){
  p_yfshandler = new YFS::YFS_Handler();
  if(p_yfshandler->Mode()!=YFS::yfsmode::off) {
    msg_Info()<<"Initialized YFS for Soft Photon Resummation"<<std::endl;
    for (const auto &pdf: m_pdflibs) {
      if(pdf!="None") THROW(fatal_error,"Cannot use PDFs with initial state YFS. Disable the PDF (PDF_LIBRARY: None) or YFS (YFS: MODE: OFF)");
    }
    for (size_t beam=0;beam<2;++beam) {
      p_yfshandler->SetInFlav(m_bunch_particles[beam]);
      p_yfshandler->SetBeam(p_beamspectra);
    }
  }
  return true;
}

bool Initialization_Handler::InitializeThePDFs()
{
  Settings& settings = Settings::GetMainSettings();
  // Load all necessary PDF libraries and check that they are available
  LoadPDFLibraries(settings);
  // Define bunch flavours
  DefineBunchFlavours(settings);
  // Initialisation of PDF sets
  for (size_t i = 1; i < 4; ++i) InitISRHandler((isr::id)(i), settings);
  msg_Info() << "Initializing PDFs ..." << endl;
  for (const auto& it : m_isrhandlers) {
    size_t pid      = it.first;
    PDF::isr::id pc = it.second->Id();
    if (pid==1) {
      msg_Info()<<"  Hard scattering:    ";

    }
    if (pid==2) {
      msg_Info()<<"  MPI:                ";

    }
    if (pid == 3) {
      msg_Info() << "  Beam re-scattering: ";
    }
    msg_Info() << (it.second->PDF(0) ? it.second->PDF(0)->Set() : "None")
               << " + "
               << (it.second->PDF(1) ? it.second->PDF(1)->Set() : "None")
               << "\n";
  }
  return true;
}

void Initialization_Handler::LoadPDFLibraries(Settings& settings) {
  std::vector<std::string> pdflibs = settings["PDF_LIBRARY"].GetVector<std::string>();
  std::vector<std::string> mpilibs = settings["MPI_PDF_LIBRARY"].GetVector<std::string>();
  std::vector<std::string> bbrlibs = settings["BBR_PDF_LIBRARY"].GetVector<std::string>();
  m_defsets[PDF::isr::hard_process]    = std::array<std::string, 2>();
  m_defsets[PDF::isr::hard_subprocess] = std::array<std::string, 2>();
  m_defsets[PDF::isr::bunch_rescatter] = std::array<std::string, 2>();
  msg_Out()<<METHOD<<": "
	   <<"bunch(0) = "<<p_beamspectra->GetBeam(0)->Bunch(0)<<" & "
	   <<"bunch(1) = "<<p_beamspectra->GetBeam(1)->Bunch(0)<<".\n";
  for (size_t beam=0;beam<2;++beam) {
    /////////////////////////////////////////////////////////
    // define bunch particle-dependent PDF libraries and sets here
    /////////////////////////////////////////////////////////
    std::string deflib("None"), defset;
    if (p_beamspectra->GetBeam(beam)->Bunch(0).Kfcode()==kf_p_plus) {
      deflib = PDF::pdfdefs->DefaultPDFLibrary(kf_p_plus);
      defset = PDF::pdfdefs->DefaultPDFSet(kf_p_plus);
    }
    else if (p_beamspectra->GetBeam(beam)->Bunch(0).Kfcode()==kf_e ||
	     p_beamspectra->GetBeam(beam)->Bunch(0).Kfcode()==kf_mu) {
      deflib = PDF::pdfdefs->DefaultPDFLibrary(kf_e);
      defset = PDF::pdfdefs->DefaultPDFSet(kf_e);
    }
    else if (p_beamspectra->GetBeam(beam)->Bunch(0).IsPhoton()) {
      deflib = PDF::pdfdefs->DefaultPDFLibrary(kf_photon);
      defset = PDF::pdfdefs->DefaultPDFSet(kf_photon);
    } else if (p_beamspectra->GetBeam(beam)->Bunch(0).Kfcode() == kf_pomeron) {
      deflib = PDF::pdfdefs->DefaultPDFLibrary(kf_pomeron);
      defset = PDF::pdfdefs->DefaultPDFSet(kf_pomeron);
    } else if (p_beamspectra->GetBeam(beam)->Bunch(0).Kfcode() == kf_reggeon) {
      deflib = PDF::pdfdefs->DefaultPDFLibrary(kf_reggeon);
      defset = PDF::pdfdefs->DefaultPDFSet(kf_reggeon);
    }
    // fix PDFs and default sets for the hard_process here
    if (pdflibs.empty()) m_pdflibs.insert(deflib);
    else m_pdflibs.insert(pdflibs[Min(beam,pdflibs.size()-1)]);
    m_defsets[PDF::isr::hard_process][beam] = defset;
    // fix PDFs and default sets for the MPI's / hard_subprocesses here
    // we may have to define defaults here.
    if (!mpilibs.empty()) {
      std::string libname = mpilibs[Min(beam,mpilibs.size()-1)];
      if (m_pdflibs.find(libname)==m_pdflibs.end()) m_pdflibs.insert(libname);
    } else if (!p_beamspectra->GetBeam(beam)->Bunch(0).IsLepton())
      m_pdflibs.insert(deflib);
    m_defsets[PDF::isr::hard_subprocess][beam] = defset;
    // fix PDFs and default sets for the beam rescattering here
    // EPA is the only configuration at the moment where we allow
    // additional scattering/interactions of the incoming beams.
    // we default to beam-beam rescatterings
    if (m_mode == eventtype::StandardPerturbative &&
        settings["BEAM_RESCATTERING"].Get<string>() == string("Amisic") &&
        p_beamspectra->On()) {
      kf_code kfbeam = p_beamspectra->GetBeam(beam)->Beam().Kfcode();
      if (!bbrlibs.empty()) {
        std::string libname = bbrlibs[Min(beam, bbrlibs.size() - 1)];
        if (m_pdflibs.find(libname) == m_pdflibs.end())
          m_pdflibs.insert(libname);
      } else
        m_pdflibs.insert(PDF::pdfdefs->DefaultPDFLibrary(kfbeam));
      defset = PDF::pdfdefs->DefaultPDFSet(kfbeam);
    }
    m_defsets[PDF::isr::bunch_rescatter][beam] = defset;
  }
  msg_Out()<<"   * (mode = "<<m_mode<<"): "
	   <<"beam 0 = "<<m_defsets[PDF::isr::hard_process][0]<<" + "
	   <<m_defsets[PDF::isr::hard_subprocess][0]<<", "
	   <<"beam 0 = "<<m_defsets[PDF::isr::hard_process][1]<<" + "
	   <<m_defsets[PDF::isr::hard_subprocess][1]<<"\n";
  // add LHAPDF if necessary and load the relevant libraries
  if (Variations::NeedsLHAPDF6Interface()) {
    m_pdflibs.insert("LHAPDFSherpa");
  }
  for (set<string>::iterator pdflib=m_pdflibs.begin(); pdflib!=m_pdflibs.end();++pdflib) {
    if (*pdflib=="None") continue;
    if (*pdflib=="LHAPDFSherpa") {
      #ifdef USING__LHAPDF
        s_loader->AddPath(std::string(LHAPDF_PATH)+"/lib");
        s_loader->LoadLibrary("LHAPDF");
      #else
        THROW(fatal_error, "LHAPDF has not been enabled during configuration.")
      #endif
    }
    void *init(s_loader->GetLibraryFunction(*pdflib,std::string("InitPDFLib")));
    if (init==NULL) THROW(fatal_error,"Cannot load PDF library "+*pdflib);
    ((PDF_Init_Function)init)();
  }

  // PDF set listing output
  int helpi{ settings["SHOW_PDF_SETS"].Get<int>() };
  if (helpi>0) {
    msg->SetLevel(2);
    PDF::PDF_Base::ShowSyntax(helpi);
    THROW(normal_exit,"Syntax shown.");
  }
}

void Initialization_Handler::InitISRHandler(const PDF::isr::id & pid,Settings& settings) {
  if (m_isrhandlers.find(pid)!=m_isrhandlers.end()) delete m_isrhandlers[pid];
  bool needs_resc  = settings["BEAM_RESCATTERING"].Get<string>()!=string("None");
  /////////////////////////////////////////////////////////////
  // make sure rescatter ISR bases are only initialised if necessary
  /////////////////////////////////////////////////////////////
  if (pid == PDF::isr::bunch_rescatter && !needs_resc) return;
  std::string tag  = ( pid==PDF::isr::hard_process ? string("PDF_SET") :
		       pid==PDF::isr::hard_subprocess ? string("MPI_PDF_SET") :
		       string("BBR_PDF_SET") );
  std::string vtag = tag + string("_VERSIONS");
  /////////////////////////////////////////////////////////////
  // read sets and versions for relevant part of event generation
  // (hard process, MPI, beam rescattering)
  /////////////////////////////////////////////////////////////
  std::vector<std::string> sets = settings[tag].GetVector<std::string>();
  std::vector<int> versions     = settings[vtag].GetVector<int>();
  if (sets.size() > 2) {
    THROW(fatal_error, "You can not specify more than two PDF sets.");
  }
  if (versions.size() > 2) {
    THROW(fatal_error, "You can not specify more than two PDF set versions.");
  }
  std::array<ISR_Base*, 2> isrbases = {};
  for (size_t beam = 0; beam < 2; beam++) {
    isrbases[beam] = nullptr;
    // fix actual set and version for beam number and part of
    // event generation here - assume central (i.e. version 0 by default)
    std::string set =
            ((sets.empty() || sets[Min(beam, sets.size() - 1)] == "Default")
                     ? m_defsets[pid][beam]
                     : sets[Min(beam, sets.size() - 1)]);
    int version =
            (versions.empty() ? 0 : versions[Min(beam, versions.size() - 1)]);
    int order = -1, scheme = -1;
    // special treatment of electron PDF
    if ((pid == PDF::isr::hard_process || pid == PDF::isr::hard_subprocess) &&
        set == "PDFe") {
      order  = settings["ISR_E_ORDER"].Get<int>();
      scheme = settings["ISR_E_SCHEME"].Get<int>();
    }
    // Initialise the actual PDF and make sure it arrives at the right place.
    Flavour flav;
    m_bunchtags[pid] = {0, 0};
    if (pid != PDF::isr::bunch_rescatter) flav = m_bunch_particles[beam];
    else {
      // use a mode to select which beam/bunch combination is used for
      // rescattering defaults to 0. The modes are
      // Mode | Beam1 - Beam2
      // -----|--------------
      // 0    | beam  - beam
      // 1    | bunch - beam
      // 2    | beam  - bunch
      // 3    | bunch - bunch (disabled, because that would be normal MPI)
      int bbrmode = settings["BBR_MODE"].Get<int>();
      if (bbrmode < 0 || bbrmode > 2)
        THROW(fatal_error, "BBR_Mode: Invalid mode, use MPI settings instead. ")
      flav                   = bbrmode & (beam + 1) ? m_bunch_particles[beam]
                                                    : p_beamspectra->GetBeam(beam)->Beam();
      m_bunchtags[pid][beam] = bbrmode & (beam + 1) ? 0 : 1;
    }
    PDF_Arguments args = PDF_Arguments(flav, beam, set, version, order, scheme);
    if (pid != PDF::isr::bunch_rescatter) {
      msg_Out()<<METHOD<<"("<<pid<<"): before getter for set = "<<set<<".\n";
      PDF_Base* pdfbase = PDF_Base::PDF_Getter_Function::GetObject(set, args);
      if (m_bunch_particles[beam].IsHadron() && pdfbase == nullptr)
        THROW(critical_error,
              "PDF '" + set + "' does not exist in any of the loaded" +
                      " libraries for " + ToString(m_bunch_particles[beam]) +
                      " bunch.");
      if (pid == PDF::isr::hard_process) rpa->gen.SetPDF(beam, pdfbase);
      msg_Out()<<METHOD<<"[beam = "<<beam<<", "<<pid<<"]: "<<m_bunch_particles[beam]
	       <<" --> "<<pdfbase<<"\n"
	       <<"   from flav = "<<flav<<", set = "<<set<<", version = "<<version<<", "
	       <<"order = "<<order<<", scheme = "<<scheme<<"\n";
      if (pdfbase == nullptr) {
        isrbases[beam] = new Intact(flav);
        needs_resc     = false;
      } else {
        pdfbase->SetBounds();
        isrbases[beam] = new Structure_Function(pdfbase, flav);
      }
      ATOOLS::rpa->gen.SetBunch(m_bunch_particles[beam], beam);
    } else if (pid == PDF::isr::bunch_rescatter && needs_resc) {
      PDF_Base* pdfbase = PDF_Base::PDF_Getter_Function::GetObject(set, args);
      if (pdfbase == nullptr)
        THROW(critical_error,
              "PDF '" + set +
                      "' for rescattering does not exist in any of the loaded" +
                      " libraries for " + ToString(flav) + " bunch.");
      pdfbase->SetBounds();
      isrbases[beam] = new Structure_Function(pdfbase, flav);
    }
  }
  if ((pid == PDF::isr::bunch_rescatter && needs_resc) ||
      pid != isr::bunch_rescatter) {
    msg_Out()<<METHOD<<"["<<pid<<"]: "
	     <<m_bunch_particles[0]<<" ("<<isrbases[0]<<") & "
	     <<m_bunch_particles[1]<<" ("<<isrbases[1]<<").\n";
    ISR_Handler* isr = new ISR_Handler(isrbases, pid);
    for (size_t beam = 0; beam < 2; beam++)
      isr->SetBeam(p_beamspectra->GetBeam(beam), beam);
    isr->Init();
    if (!(p_beamspectra->CheckConsistency(m_bunch_particles))) {
      msg_Error() << "Error in Environment::InitializeThePDFs()" << endl
                  << "   Inconsistent ISR & Beam:" << endl
                  << "   Abort program." << endl;
      Abort();
    }
    isr->Output();
    m_isrhandlers[pid] = isr;
  }
}

void Initialization_Handler::DefineBunchFlavours(Settings& settings) {
  std::vector<int> bunches{ settings["BUNCHES"].GetVector<int>() };
  if (bunches.size() > 2) {
    THROW(fatal_error, "You can not specify more than two bunches.");
  }
  for (size_t beam=0;beam<2;beam++) {
    if (bunches.empty()) m_bunch_particles[beam] = p_beamspectra->GetBeam(beam)->Bunch(0);
    else {
      int flav = bunches[Min(beam,bunches.size()-1)];
      m_bunch_particles[beam] = Flavour((kf_code)abs(flav));
      if (flav<0) m_bunch_particles[beam] = m_bunch_particles[beam].Bar();
    }
  }
}

bool Initialization_Handler::InitializeTheRemnants() {
  ///////////////////////////////////////////////////////////
  // define two sets of remnants, if necessary (i.e. if we have bunch rescattering).
  // we will have to make sure that the remnants for hard process and hard subprocess -
  // the MPI related to the hard process - are the same.
  // I have the feeling we will have to communicate the mode to the Remnant_Handler in question
  ///////////////////////////////////////////////////////////
  REMNANTS::Remnants_Parameters();
  m_remnanthandlers[isr::hard_process] =
    new Remnant_Handler(m_isrhandlers[isr::hard_process],p_yfshandler,p_beamspectra,
			m_bunchtags[isr::hard_process]);
  m_remnanthandlers[isr::hard_subprocess] = m_remnanthandlers[isr::hard_process];
  // Set up the correct combination of remnant handlers for the bunch rescattering, depending on the
  // mode
  if (m_isrhandlers.find(isr::bunch_rescatter)!=m_isrhandlers.end()) {
    int bbrmode = Settings::GetMainSettings()["BBR_MODE"].Get<int>();
    if (bbrmode == 0)
      m_remnanthandlers[isr::bunch_rescatter] =
      new Remnant_Handler(m_isrhandlers[isr::bunch_rescatter],p_yfshandler,p_beamspectra,
                                  m_bunchtags[isr::bunch_rescatter]);
    else {
      std::array<std::shared_ptr<REMNANTS::Remnant_Base>, 2> remnants = {
              nullptr, nullptr};
      // For mixed beam-bunch rescattering, the remnant for the bunch has already been created for
      // the hard process; we need to get it from there and hand it over
      remnants[bbrmode - 1] = m_remnanthandlers[PDF::isr::hard_process]->GetRemnant(bbrmode - 1);
      m_remnanthandlers[isr::bunch_rescatter] =
              new Remnant_Handler(remnants, m_isrhandlers[isr::bunch_rescatter],
                                  p_yfshandler, p_beamspectra,
                                  m_bunchtags[isr::bunch_rescatter]);
    }
  }
  msg_Info()<<"Initializing remnants ...\n"
	    <<"  Hard process: "
	    <<m_remnanthandlers[isr::hard_process]->GetRemnant(0)->GetBeam()->Bunch(0)<<" ("
	    <<m_remnanthandlers[isr::hard_process]->GetRemnant(0)->Type()<<") + "
	    <<m_remnanthandlers[isr::hard_process]->GetRemnant(1)->GetBeam()->Bunch(0)<<" ("
	    <<m_remnanthandlers[isr::hard_process]->GetRemnant(1)->Type()<<")\n";
  if (m_remnanthandlers.find(isr::bunch_rescatter) != m_remnanthandlers.end())
    msg_Info()<<"  Rescattering: "
	      <<m_remnanthandlers[isr::bunch_rescatter]->GetRemnant(0)->GetBeam()->Bunch(m_bunchtags[isr::bunch_rescatter][0])<<" ("
	      <<m_remnanthandlers[isr::bunch_rescatter]->GetRemnant(0)->Type()<<") + "
	      <<m_remnanthandlers[isr::bunch_rescatter]->GetRemnant(1)->GetBeam()->Bunch(m_bunchtags[isr::bunch_rescatter][1])<<" ("
	      <<m_remnanthandlers[isr::bunch_rescatter]->GetRemnant(1)->Type()<<")\n";
  return true;
}

bool Initialization_Handler::InitializeTheHardDecays()
{
  if (!Settings::GetMainSettings()["HARD_DECAYS"]["Enabled"].Get<bool>())
    return true;
  if (p_harddecays) {
    delete p_harddecays;
    p_harddecays = NULL;
  }
  p_harddecays = new Hard_Decay_Handler();
  return true;
}

bool Initialization_Handler::InitializeTheMatrixElements()
{
  msg_Info()<<"Initializing matrix elements for the hard processes ...\n";
#ifdef USING__EWSud
  // in case that KFACTOR=EWSud is used we need to be ready when the ME handler
  // sets up the KFactor setters
  if (!s_loader->LoadLibrary("SherpaEWSud"))
    THROW(missing_module,"Cannot load EWSud library.");
#endif
  if (p_mehandler) delete p_mehandler;
  p_mehandler = new Matrix_Element_Handler(p_model);
  p_mehandler->SetShowerHandler(m_showerhandlers[isr::hard_process]);
  p_mehandler->SetRemnantHandler(m_remnanthandlers[isr::hard_process]);
  auto ret = p_mehandler->InitializeProcesses(p_beamspectra,
                                              m_isrhandlers[isr::hard_process],
                                              p_yfshandler);
  return ret==1;
}

bool Initialization_Handler::InitializeTheShowers()
{
  ///////////////////////////////////////////////////////////
  // define up to three shower handlers ...
  ///////////////////////////////////////////////////////////
  msg_Info()<<"Initializing showers ...\n";
  std::vector<isr::id> isrtypes = {isr::hard_process, isr::hard_subprocess};
  if (m_isrhandlers.find(isr::bunch_rescatter)!=m_isrhandlers.end())
    isrtypes.push_back(isr::bunch_rescatter);
  for (isr::id id : isrtypes) {
    as->SetActiveAs(id);
    if (m_showerhandlers.find(id) != m_showerhandlers.end())
      delete m_showerhandlers[id];
    m_showerhandlers[id] = new Shower_Handler(p_model, m_isrhandlers[id], id-1);
    m_showerhandlers[id]->SetRemnants(m_remnanthandlers[id]);
    m_isrhandlers[id]->SetRemnants(m_remnanthandlers[id]->GetRemnants());
  }
  as->SetActiveAs(isr::hard_process);
  return true;
}


bool Initialization_Handler::InitializeTheUnderlyingEvents()
{
  ///////////////////////////////////////////////////////////
  // define up to three multiple interaction handlers ...
  ///////////////////////////////////////////////////////////
  if (m_isrhandlers[isr::hard_process]->Mode() != PDF::isrmode::hadron_hadron)
    Settings::GetMainSettings()["MI_HANDLER"].OverrideScalar<std::string>(
            "None");
  std::vector<isr::id> isrtypes = {isr::hard_subprocess};
  if (m_isrhandlers.find(isr::bunch_rescatter) != m_isrhandlers.end())
    isrtypes.push_back(isr::bunch_rescatter);
  for (isr::id id : isrtypes) {
    as->SetActiveAs(isr::hard_subprocess);
    msg_Out()<<METHOD<<"["<<isrtypes[i]<<"]:\n";
    MI_Handler * mih = new MI_Handler(p_model,m_isrhandlers[id], p_yfshandler, m_remnanthandlers[id]);
    mih->SetShowerHandler(m_showerhandlers[id]);
    as->SetActiveAs(isr::hard_process);// really needed?
    m_mihandlers[id] = mih;
  }
  msg_Info()<<"Underlying event/multiple interactions initialized\n";
  for (isr::id id : isrtypes) {
    MI_Handler * mih = m_mihandlers[id];
    msg_Info() << "    MI[" << id << "]: on = " << mih->On() << " "
               << "(type = " << mih->Type() << ", " << mih->Name() << ")\n";
  }
  return true;
}

bool Initialization_Handler::InitializeTheSoftCollisions()
{
  ///////////////////////////////////////////////////////////
  // define up to two soft collision handlers -
  // they will have to differ in how they fill blobs.
  // modify the beam remnants to take care of the rescatter
  ///////////////////////////////////////////////////////////
  std::vector<isr::id> isrtypes;
  isrtypes.push_back(isr::hard_subprocess);
  if (m_isrhandlers.find(isr::bunch_rescatter)!=m_isrhandlers.end())
    isrtypes.push_back(isr::bunch_rescatter);
  for (isr::id id : isrtypes) {
    if (m_schandlers.find(id) != m_schandlers.end()) delete m_schandlers[id];
    MI_Handler* mih = m_mihandlers[id];
    m_schandlers[id] =
            (mih->On()
                     ? new Soft_Collision_Handler(mih->Amisic(), mih->Shrimps(),
                                                  id == isr::bunch_rescatter)
                     : nullptr);
    if (id == isr::bunch_rescatter) {
      p_beamremnants->AddBunchRescattering(
              m_remnanthandlers[isr::bunch_rescatter],
              m_schandlers[isr::bunch_rescatter]);
    }
  }
  bool did_print_header {false};
  for (isr::id id : isrtypes) {
    if (m_schandlers[id]!=NULL) {
      if (!did_print_header) {
        msg_Info()<<"Soft-collision handlers:\n";
        did_print_header = true;
      }
      msg_Info()<<"  Type["<<id<<"]: "
		<<m_schandlers[id]->Soft_CollisionModel()<<"\n";
    }
  }
  return true;
}

bool Initialization_Handler::InitializeTheBeamRemnants()
{
  msg_Info()<<"Initializing the beam remnants ...\n";
  if (p_beamremnants)  delete p_beamremnants;
  p_beamremnants = new Beam_Remnant_Handler(p_beamspectra,
					    m_remnanthandlers[isr::hard_process],
					    m_schandlers[isr::hard_subprocess]);
  return 1;
}

bool Initialization_Handler::InitializeTheColourReconnections()
{
  if (p_reconnections) { delete p_reconnections; p_reconnections = NULL; }
  p_reconnections = new Colour_Reconnection_Handler();
  p_reconnections->Output();
  return 1;
}

bool Initialization_Handler::InitializeTheFragmentation()
{
  if (p_fragmentation) { delete p_fragmentation; p_fragmentation = NULL; }
  as->SetActiveAs(isr::hard_subprocess);
  Settings& s = Settings::GetMainSettings();
  string fragmentationmodel = s["FRAGMENTATION"].Get<std::string>();
  if (fragmentationmodel!="None") {
    if (!s["BEAM_REMNANTS"].Get<bool>())
      msg_Error()<<METHOD<<om::red<<": Fragmentation called without beam remnants, "<<
        "hadronization might not be possible due to missing colour partners "<<
        "in the beam!\nFragmentation might stall, please consider aborting manually.\n"<<om::reset;
    Hadron_Init().Init();
    ATOOLS::OutputHadrons(msg->Tracking());
  }
  p_fragmentation = Fragmentation_Getter::GetObject
    (fragmentationmodel,
     Fragmentation_Getter_Parameters(m_showerhandlers[isr::hard_process]->ShowerGenerator()));
  if (p_fragmentation==NULL)
    THROW(fatal_error, "  Fragmentation model '"+fragmentationmodel+"' not found.");
  as->SetActiveAs(isr::hard_process);
  msg_Info()<<"Initialized fragmentation\n";
  return 1;
}

bool Initialization_Handler::InitializeTheHadronDecays()
{
  Settings& s = Settings::GetMainSettings();
  if (s["FRAGMENTATION"].Get<std::string>() == "None")
    return true;
  std::string decmodel{ s["HADRON_DECAYS"]["Model"].Get<std::string>() };
  msg_Tracking()<<"Decaymodel = "<<decmodel<<std::endl;
  if (decmodel=="None") return true;
  else if (decmodel==std::string("HADRONS++")) {
    as->SetActiveAs(isr::hard_subprocess);
    Hadron_Decay_Handler* hd=new Hadron_Decay_Handler();
    as->SetActiveAs(isr::hard_process);
    p_hdhandler=hd;
  }
  else {
    THROW(fatal_error,"Hadron decay model '"+decmodel+"' not implemented.");
  }
  msg_Info()<<"Initialized hadron decays (model = "
            <<decmodel<<")\n";
  return true;
}

bool Initialization_Handler::InitializeTheSoftPhotons()
{
  if (p_softphotons) { delete p_softphotons; p_softphotons = NULL; }
  p_softphotons = new Soft_Photon_Handler(p_mehandler);
  if (p_harddecays) p_harddecays->SetSoftPhotonHandler(p_softphotons);
  if (p_hdhandler)  p_hdhandler->SetSoftPhotonHandler(p_softphotons);
  msg_Info()<<"Initialized soft photons"<<endl;
  return true;
}

bool Initialization_Handler::InitializeTheAnalyses()
{
  Settings& s = Settings::GetMainSettings();
  std::string outpath=s["ANALYSIS_OUTPUT"].Get<std::string>();
  const Analysis_Arguments args{ Analysis_Arguments(outpath) };
  std::vector<std::string> analyses=s["ANALYSIS"].GetVector<std::string>();
  for (size_t i=0; i<analyses.size(); ++i) {
    if (analyses[i]=="1") analyses[i]="Internal";
    if (analyses[i]=="None") continue;
    if (analyses[i]=="Internal")
      if (!s_loader->LoadLibrary("SherpaAnalysis"))
        THROW(missing_module,"Cannot load Analysis library (-DSHERPA_ENABLE_ANALYSIS=ON).");
    if (analyses[i]=="Rivet" || analyses[i]=="RivetME" || analyses[i]=="RivetShower") {
      bool hepmc_loaded {false};
#ifdef USING__HEPMC3
      hepmc_loaded |= (s_loader->LoadLibrary("SherpaHepMC3Output") != nullptr);
#endif
      if (!hepmc_loaded) {
        THROW(missing_module,
              "Cannot load HepMC library (-DSHERPA_ENABLE_HEPMC3=ON).");
      }
      if (!s_loader->LoadLibrary("SherpaRivetAnalysis"))
        THROW(missing_module,"Cannot load RivetAnalysis library (-DSHERPA_ENABLE_RIVET=ON).");
    }
    Analysis_Interface* ana =
      Analysis_Interface::Analysis_Getter_Function::GetObject(analyses[i], args);
    if (ana==NULL) {
      if (!s_loader->LoadLibrary("Sherpa"+analyses[i]+"Analysis"))
	THROW(missing_module,"Cannot load Analysis library '"+analyses[i]+"'.");
      ana = Analysis_Interface::Analysis_Getter_Function::GetObject(
          analyses[i], args);
      if (ana==NULL) THROW(fatal_error,"Cannot initialize Analysis "+analyses[i]);
    }
    m_analyses.push_back(ana);
  }
  return true;
}

bool Initialization_Handler::InitializeTheReweighting(Variations_Mode mode)
{
  if (p_variations) {
    delete p_variations;
  }
  if (mode != Variations_Mode::nominal_only)
    Variations::CheckConsistencyWithBeamSpectra(p_beamspectra);
  p_variations = new Variations(mode);
  s_variations = p_variations;
  if (p_mehandler)
    p_mehandler->InitializeTheReweighting(mode);
  if (mode != Variations_Mode::nominal_only)
    msg_Info()<<"Initialized on-the-fly reweighting"<<endl;
  return true;
}

bool Initialization_Handler::InitializeTheFilter()
{
  if (p_filter)
    delete p_filter;
  p_filter = new Filter();
  if (!p_filter->Init()) { delete p_filter; p_filter = NULL; }
  return true;
}

bool Initialization_Handler::CalculateTheHardProcesses()
{
  if (m_mode!=eventtype::StandardPerturbative) return true;

  msg_Events()<<"===================================================================\n"
              <<"Start calculating the hard cross sections. This may take some time.\n";
  as->SetActiveAs(isr::hard_process);
  p_variations->DisableVariations();
  int ok = p_mehandler->CalculateTotalXSecs();
  p_variations->EnableVariations();
  if (ok) {
    msg_Events()<<"Calculating the hard cross sections has been successful.\n"
		<<"====================================================================\n";
  }
  else {
    msg_Events()<<"Calculating the hard cross sections failed. Check this carefully.\n"
		<<"=======================================================================\n";
  }
  return ok;
}

void Initialization_Handler::SetGlobalVariables()
{
  Settings& s = Settings::GetMainSettings();
  double sf(s["SCALE_FACTOR"].Get<double>());
  double fsf(sf*s["FACTORIZATION_SCALE_FACTOR"].Get<double>());
  double rsf(sf*s["RENORMALIZATION_SCALE_FACTOR"].Get<double>());
  double qsf(sf*s["RESUMMATION_SCALE_FACTOR"].Get<double>());
  rpa->gen.SetVariable("FACTORIZATION_SCALE_FACTOR", ToString(fsf));
  rpa->gen.SetVariable("RENORMALIZATION_SCALE_FACTOR", ToString(rsf));
  rpa->gen.SetVariable("RESUMMATION_SCALE_FACTOR", ToString(qsf));
  msg_Debugging()<<METHOD<<"(): Set scale factors {\n"
		 <<"  fac scale: "<<rpa->gen.Variable("FACTORIZATION_SCALE_FACTOR")<<"\n"
		 <<"  ren scale: "<<rpa->gen.Variable("RENORMALIZATION_SCALE_FACTOR")<<"\n"
		 <<"  res scale: "<<rpa->gen.Variable("RESUMMATION_SCALE_FACTOR")<<"\n}\n";

  // TODO: remove from rpa?
  double virtfrac = s["VIRTUAL_EVALUATION_FRACTION"].Get<double>();
  rpa->gen.SetVariable("VIRTUAL_EVALUATION_FRACTION", ToString(virtfrac));

  std::string evtm{ s["EVENT_GENERATION_MODE"].Get<std::string>() };
  int eventmode{ 0 };
  if (evtm=="Unweighted" || evtm=="U") eventmode=1;
  else if (evtm=="PartiallyUnweighted" || evtm=="P") eventmode=2;
  rpa->gen.SetVariable("EVENT_GENERATION_MODE",ToString(eventmode));
}
