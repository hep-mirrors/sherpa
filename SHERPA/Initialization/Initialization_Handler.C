#include "SHERPA/Initialization/Initialization_Handler.H"

#include "SHERPA/PerturbativePhysics/Hard_Decay_Handler.H"
#include "SHERPA/PerturbativePhysics/Shower_Handler.H"
#include "SHERPA/SoftPhysics/Beam_Remnant_Handler.H"
#include "SHERPA/SoftPhysics/Fragmentation_Handler.H"
#include "SHERPA/SoftPhysics/Hadron_Decay_Handler.H"
#include "SHERPA/SoftPhysics/Lund_Decay_Handler.H"
#include "SHERPA/SoftPhysics/Soft_Collision_Handler.H"
#include "SHERPA/PerturbativePhysics/MI_Handler.H"
#include "SHERPA/SoftPhysics/Soft_Photon_Handler.H"
#include "SHERPA/LundTools/Lund_Interface.H"
#include "SHERPA/Tools/Event_Reader_Base.H"
#include "SHERPA/Main/Filter.H"
#include "PHASIC++/Scales/Core_Scale_Setter.H"
#include "MODEL/Main/Model_Base.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "METOOLS/Currents/C_Spinor.H"
#include "PDF/Main/Structure_Function.H"
#include "PDF/Main/Intact.H"
#include "PDF/Main/PDF_Base.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include "ATOOLS/Math/Scaling.H"
#include "ATOOLS/Phys/Spinor.H"
#include "ATOOLS/Phys/Variations.H"
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

#include <sys/stat.h>
#include <time.h>

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
  p_mehandler(NULL), p_harddecays(NULL), p_beamremnants(NULL),
  p_fragmentation(NULL), p_softcollisions(NULL), p_hdhandler(NULL), 
  p_mihandler(NULL), p_softphotons(NULL), p_evtreader(NULL),
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

  rpa->gen.SetSoftSC(s["SOFT_SPIN_CORRELATIONS"].Get<int>());
  rpa->gen.SetHardSC(s["HARD_SPIN_CORRELATIONS"].Get<int>());

  exh->AddTerminatorObject(this);
}

void Initialization_Handler::RegisterDefaults()
{
  Settings& s = Settings::GetMainSettings();
  s["EVENT_GENERATION_MODE"].SetDefault("PartiallyUnweighted");
  s["EVENT_TYPE"].SetDefault("StandardPerturbative");
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
  s["DECAYMODEL"]
    .SetDefault((frag == "Lund") ? "Lund" : "Hadrons")
    .UseNoneReplacements();

  s["SOFT_SPIN_CORRELATIONS"].SetDefault(0);
  auto hdenabled = s["HARD_DECAYS"]["Enabled"].Get<bool>();
  s["HARD_SPIN_CORRELATIONS"].SetDefault(hdenabled);

  s["EVENT_INPUT"].SetDefault("");
  s["STATUS_PATH"].SetDefault("");
  s["PATH"].SetDefault("");
  s["SAVE_STATUS"].SetDefault("");
  s["MI_PDF_LIBRARY"].SetDefault("");

  s.DeclareVectorSettingsWithEmptyDefault({
      "EVENT_OUTPUT",
      "ANALYSIS",
      "SHERPA_LDADD",
      "SHERPA_VERSION",
      "PDF_LIBRARY",
      "PDF_SET",
      "MPI_PDF_SET",
      "VARIATIONS",
      "BUNCHES",
      "PDF_SET_VERSIONS",
      "MPI_PDF_SET_VERSIONS",
      "NLO_CSS_DISALLOW_FLAVOUR",
      "MASSIVE_PS",
      "MASSLESS_PS"
      });
  s.DeclareMatrixSettingsWithEmptyDefault({ "CSS_ENHANCE" });
  s["EVENT_OUTPUT"].UseNoneReplacements();
  s["VARIATIONS"].UseNoneReplacements();
  s["PDF_LIBRARY"].UseNoneReplacements();
  s["ANALYSIS"].UseNoneReplacements();

  s["SHOW_ME_GENERATORS"].SetDefault(0);
  s["SHOW_PS_GENERATORS"].SetDefault(0);
  s["SHOW_NLOMC_GENERATORS"].SetDefault(0);
  s["SHOW_SHOWER_GENERATORS"].SetDefault(0);
  s["SHOW_SCALE_SYNTAX"].SetDefault(0);
  s["SHOW_SELECTOR_SYNTAX"].SetDefault(0);
  s["SHOW_MODEL_SYNTAX"].SetDefault(0);
  s["SHOW_FILTER_SYNTAX"].SetDefault(0);
  s["SHOW_ANALYSIS_SYNTAX"].SetDefault(0);
  s["SHOW_VARIABLE_SYNTAX"].SetDefault(0);
  s["SHOW_PDF_SETS"].SetDefault(0);

  s["ISR_SMIN"].SetDefault(1e-10);
  s["ISR_SMAX"].SetDefault(1.0);
  s["ISR_E_ORDER"].SetDefault(1);
  s["ISR_E_SCHEME"].SetDefault(2);

  s["KFACTOR"].SetDefault("None").UseNoneReplacements();
  s["SCALES"].SetDefault("METS{MU_F2}{MU_R2}{MU_Q2}");
  s["SCALE_FACTOR"].SetDefault(1.0);
  s["FACTORIZATION_SCALE_FACTOR"].SetDefault(1.0);
  s["RENORMALIZATION_SCALE_FACTOR"].SetDefault(1.0);
  s["USR_WGT_MODE"].SetDefault(true);

  s["OVERRIDE_PDF_INFO"].SetDefault(false);

  s["NLO_SUBTRACTION_SCHEME"].SetDefault(0);

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
  s["SHOWER_GENERATOR"].SetDefault("Dire").UseNoneReplacements();
  std::string showergen{ s["SHOWER_GENERATOR"].Get<std::string>() };
  s["JET_CRITERION"].SetDefault(showergen);
  s["NLOMC_GENERATOR"].SetDefault(showergen);
  s["CSS_EVOLUTION_SCHEME"].SetDefault(1);
  s["CSS_KFACTOR_SCHEME"].SetDefault(9);
  s["CSS_SCALE_SCHEME"].SetDefault(0);
  s["CSS_SCALE_VARIATION_SCHEME"].SetDefault(1);
  // TODO: Should this be set to 3.0 for the new Dire default? See the manual
  // Sherpa section on master for details
  s["CSS_FS_PT2MIN"].SetDefault(2.0);
  s["CSS_IS_PT2MIN"].SetDefault(2.0);
  s["CSS_FS_AS_FAC"].SetDefault(1.0);
  s["CSS_IS_AS_FAC"].SetDefault(1.0);
  s["CSS_SCALE_FACTOR"].SetDefault(1.);
  s["CSS_MASS_THRESHOLD"].SetDefault(0.0);
  s["VIRTUAL_EVALUATION_FRACTION"].SetDefault(1.0);
  s["CSS_RECO_CHECK"].SetDefault(0);
  s["CSS_MAXEM"].SetDefault(std::numeric_limits<size_t>::max());
  s["CSS_REWEIGHT_ALPHAS"].SetDefault(1);
  s["CSS_REWEIGHT_PDFS"].SetDefault(1);
  s["REWEIGHT_MAXEM"].SetDefault(std::numeric_limits<unsigned int>::max());
  s["REWEIGHT_MCATNLO_EM"].SetDefault(1);
  s["CSS_REWEIGHT_SCALE_CUTOFF"].SetDefault(5.0);
  s["CSS_KIN_SCHEME"].SetDefault(1);
  s["NLO_CSS_KIN_SCHEME"].SetDefault(1);
  s["CSS_OEF"].SetDefault(3.0);
  s["CSS_KMODE"].SetDefault(2);
  s["CSS_RESPECT_Q2"].SetDefault(false);
  s["CSS_CKFMODE"].SetDefault(1);
  s["CSS_PDFCHECK"].SetDefault(1);
  s["CSS_QCD_MODE"].SetDefault(1);
  s["CSS_EW_MODE"].SetDefault(false);
  s["CSS_USE_BBW"].SetDefault(1);
  s["CSS_RECO_DECAYS"].SetDefault(0);
  s["CSS_MAXPART"].SetDefault(std::numeric_limits<int>::max());
  s["CSS_PDF_MIN"].SetDefault(1.0e-4);
  s["CSS_PDF_MIN_X"].SetDefault(1.0e-2);
  s["CSS_WEIGHT_CHECK"].SetDefault(false);
  s["CSS_CMODE"].SetDefault(1);
  s["CSS_NCOL"].SetDefault(3);
  s["CSS_RECALC_FACTOR"].SetDefault(4.0);
  s["CSS_TC_ENHANCE"].SetDefault(1.0);
  s["CSS_COUPLING_SCHEME"].SetDefault(1);
  s["CSS_ME_CORRECTION"].SetDefault(0);
  s["CSS_KERNEL_TYPE"].SetDefault(15);
  s["NLO_CSS_RECALC_FACTOR"].SetDefault(2.0);
  s["NLO_CSS_PSMODE"].SetDefault(0);
  s["NLO_CSS_WEIGHT_CHECK"].SetDefault(0);
  s["NLO_CSS_MAXEM"].SetDefault(1);
  s["MI_CSS_KFACTOR_SCHEME"].SetDefault(0);
  s["MI_CSS_IS_PT2MIN"].SetDefault(4.0);
  s["MI_CSS_FS_PT2MIN"].SetDefault(1.0);
  s["MI_CSS_IS_AS_FAC"].SetDefault(0.66);
  s["MI_CSS_FS_AS_FAC"].SetDefault(0.66);
  s["MI_CSS_KIN_SCHEME"].SetDefault(1);

  s["COMIX_DEFAULT_GAUGE"].SetDefault(1);

  s["DIPOLES"]["KAPPA"].SetDefault(2.0/3.0);

  s["COUPLINGS"].SetDefault("Alpha_QCD 1");

  s["EXTRAXS_CSS_APPROX_ME"].SetDefault(false);

  s["RESPECT_MASSIVE_FLAG"].SetDefault(false);
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
  if (p_fragmentation) { delete p_fragmentation; p_fragmentation = NULL; }
  if (p_beamremnants)  { delete p_beamremnants;  p_beamremnants  = NULL; }
  if (p_harddecays)    { delete p_harddecays;    p_harddecays    = NULL; }
  if (p_hdhandler)     { delete p_hdhandler;     p_hdhandler     = NULL; }
  if (p_softphotons)   { delete p_softphotons;   p_softphotons   = NULL; } 
  if (p_softcollisions){ delete p_softcollisions;p_softcollisions= NULL; } 
  if (p_mihandler)     { delete p_mihandler;     p_mihandler     = NULL; }
  if (p_remnants)      { delete p_remnants;      p_remnants      = NULL; }
  if (p_beamspectra)   { delete p_beamspectra;   p_beamspectra   = NULL; }
  if (p_model)         { delete p_model;         p_model         = NULL; }
  if (p_variations)    { delete p_variations;    p_variations    = NULL; }
  if (p_filter)        { delete p_filter;        p_filter        = NULL; }
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
  okay = okay && InitializeTheModel();

  if (m_mode==eventtype::StandardPerturbative) {
    std::string eventtype{ s["EVENT_TYPE"].Get<std::string>() };
    if (eventtype=="StandardPerturbative")
      m_mode=eventtype::StandardPerturbative;
    else if (eventtype=="MinimumBias") {
      m_mode=eventtype::MinimumBias;
      s["MI_HANDLER"].OverrideScalar<std::string>("None");
    }
    else if (eventtype=="HadronDecay") {
      m_mode=eventtype::HadronDecay;
      s["MI_HANDLER"].OverrideScalar<std::string>("None");
    }
    else {
      THROW(not_implemented,"Unknown event type '"+eventtype+"'");
    }
  }
  okay = okay && InitializeTheBeams();
  okay = okay && InitializeThePDFs();
  okay = okay && InitializeTheRemnants();
  if (!p_model->ModelInit(m_isrhandlers))
    THROW(critical_error,"Model cannot be initialized");
  okay = okay && p_beamspectra->Init();
  p_model->InitializeInteractionModel();
  okay = okay && InitializeTheAnalyses();
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
      THROW(missing_module,"Cannot load output library Sherpa"+libname+"Input.");
    p_evtreader = Event_Reader_Base::Getter_Function::GetObject
      (m_evtform,Input_Arguments(s.GetPath(), infile,
				 p_model, m_isrhandlers[isr::hard_process]));
    if (p_evtreader==NULL) THROW(fatal_error,"Event reader not found");
    msg_Events()<<"SHERPA will read in the events."<<std::endl
  		<<"   The full framework is not needed."<<std::endl;
    InitializeTheHardDecays();
    InitializeTheBeamRemnants();
    InitializeTheIO();
    InitializeTheReweighting();
    return true;
  }
  PHASIC::Phase_Space_Handler::GetInfo();
  if (rpa->gen.NumberOfEvents()>0) {
  }
  okay = okay && InitializeTheShowers();
  okay = okay && InitializeTheMatrixElements();
  okay = okay && InitializeTheBeamRemnants();
  okay = okay && InitializeTheHardDecays();
  //  only if events:
  if (rpa->gen.NumberOfEvents()>0) {
    okay = okay && InitializeTheFragmentation();
    okay = okay && InitializeTheSoftCollisions();
    okay = okay && InitializeTheHadronDecays();
    okay = okay && InitializeTheUnderlyingEvents();
    okay = okay && InitializeTheSoftPhotons();
    okay = okay && InitializeTheIO();
    okay = okay && InitializeTheReweighting();
    okay = okay && InitializeTheFilter();
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
  smin = Max(smin,p_beamspectra->SprimeMin());
  smax = Min(smax,p_beamspectra->SprimeMax());
  if (m_isrhandlers[isr::hard_process]->On()) {
    smin = Max(smin,m_isrhandlers[isr::hard_process]->SprimeMin());
    smax = Min(smax,m_isrhandlers[isr::hard_process]->SprimeMax());
  }
  if (p_beamspectra->On()) {
    p_beamspectra->SetSprimeMin(smin);
  }
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
  p_beamspectra->Output();
  return 1;
}

bool Initialization_Handler::InitializeThePDFs()
{
  Settings& s = Settings::GetMainSettings();
  // load PDF libraries
  std::string defset[2];
  for (int beam(0);beam<=1;++beam) {
    std::string deflib("None");
    if (p_beamspectra->GetBeam(beam)->Bunch().Kfcode()==kf_p_plus) {
      deflib="NNPDFSherpa";
      defset[beam]="NNPDF31_nnlo_as_0118_mc";
    }
    else if (p_beamspectra->GetBeam(beam)->Bunch().Kfcode()==kf_e) {
      deflib="PDFESherpa";
      defset[beam]="PDFe";
    }
    else if (p_beamspectra->GetBeam(beam)->Bunch().IsPhoton()) {
      deflib="GRVSherpa";
      defset[beam]="GRV";
    }
    std::vector<std::string> pdflibs{
      s["PDF_LIBRARY"].GetVector<std::string>() };
    if (pdflibs.size()==0) m_pdflibs.insert(deflib);
    for (size_t i(0);i<pdflibs.size();++i) m_pdflibs.insert(pdflibs[i]);
    std::string mpilib{ s["MI_PDF_LIBRARY"].Get<std::string>() };
    if (mpilib != "") {
      m_pdflibs.insert(mpilib);
    }
  }
  if (Variations::NeedsLHAPDF6Interface()) {
    m_pdflibs.insert("LHAPDFSherpa");
  }
  for (set<string>::iterator pdflib=m_pdflibs.begin(); pdflib!=m_pdflibs.end();
       ++pdflib) {
    if (*pdflib=="None") continue;
    if (*pdflib=="LHAPDFSherpa") {
#ifdef USING__LHAPDF
      s_loader->AddPath(std::string(LHAPDF_PATH)+"/lib");
      s_loader->LoadLibrary("LHAPDF");
#else
      THROW(fatal_error, "Sherpa not compiled with LHAPDF support.");
#endif
    }
    void *init(s_loader->GetLibraryFunction(*pdflib,"InitPDFLib"));
    if (init==NULL) THROW(fatal_error,"Cannot load PDF library "+*pdflib);
    ((PDF_Init_Function)init)();
  }

  // PDF set listing output
  int helpi{ s["SHOW_PDF_SETS"].Get<int>() };
  if (helpi>0) {
    msg->SetLevel(2);
    PDF::PDF_Base::ShowSyntax(helpi);
    THROW(normal_exit,"Syntax shown.");
  }

  // Initialisation of PDF sets
  for (size_t i=0;i<2;++i) {
    isr::id id=(isr::id)(i+1);
    if (m_isrhandlers.find(id)!=m_isrhandlers.end()) 
      delete m_isrhandlers[id]; 
    PDF_Base * pdfbase;
    ISR_Base ** isrbases = new ISR_Base*[2];
    double m_bunch_splimits[2];
    for (int j=0;j<2;++j) {
      std::string indextag("_" + ToString(j + 1));

      // read bunch flavour
      int defaultflav(p_beamspectra->GetBeam(j)->Bunch());
      int flav{ 0 };
      std::vector<int> bunches{ s["BUNCHES"].GetVector<int>() };
      if (bunches.size() > 2) {
        THROW(fatal_error, "You can not specify more than two bunches.");
      } else if (bunches.size() == 0) {
        flav = defaultflav;
      } else {
        flav = bunches[Min(j, (int)bunches.size() - 1)];
      }
      m_bunch_particles[j] = Flavour((kf_code)abs(flav));
      if (flav<0) m_bunch_particles[j] = m_bunch_particles[j].Bar();

      // read PDF set
      std::string set;
      std::vector<std::string> sets{ s["PDF_SET"].GetVector<std::string>() };
      if (sets.size() > 2) {
        THROW(fatal_error, "You can not specify more than two PDF sets.");
      } else if (sets.size() == 0) {
        set = defset[j];
      } else {
        set = sets[Min(j, (int)sets.size() - 1)];
        if (set == "Default")
          set = defset[j];
      }

      // read PDF set versions
      int version(0);
      std::vector<int> versions{ s["PDF_SET_VERSIONS"].GetVector<int>() };
      if (versions.size() > 2) {
        THROW(fatal_error, "You can not specify more than two PDF set versions.");
      } else if (versions.size() == 0) {
        version = 0;
      } else {
        version = versions[Min(j, (int)versions.size() - 1)];
      }

      // override PDF set and version for hard subprocesses (multiple
      // interactions)
      bool specializedformi(false);
      if (id==isr::hard_subprocess) {

        std::vector<std::string> mpisets{
          s["MPI_PDF_SET"].GetVector<std::string>() };
        if (mpisets.size() > 2) {
          THROW(fatal_error, "You can not specify more than two MPI PDF sets.");
        } else if (mpisets.size() == 0) {
        } else {
          set = mpisets[Min(j, (int)mpisets.size() - 1)];
          // If using a special MI PDF set, then do not use normal PDF set
          // version
          version = 0;
          specializedformi = true;
        }

        // read PDF set versions
        std::vector<int> mpiversions{
          s["MPI_PDF_SET_VERSIONS"].GetVector<int>() };
        if (mpiversions.size() > 2) {
          THROW(fatal_error, "You can not specify more than two MPI PDF set versions.");
        } else if (mpiversions.size() == 0) {
        } else {
          version = mpiversions[Min(j, (int)mpiversions.size() - 1)];
          specializedformi = true;
        }

      }

      // read further configuration
      int order{ -1 };
      int scheme{ -1 };
      if (set == "PDFe") {
        order = s["ISR_E_ORDER"].Get<int>();
        scheme = s["ISR_E_SCHEME"].Get<int>();
      }

      // Load PDF
      pdfbase = PDF_Base::PDF_Getter_Function::GetObject(
          set,
          PDF_Arguments(m_bunch_particles[j], j, set, version, order, scheme));
      if (i==0) rpa->gen.SetPDF(j,pdfbase);
      if (m_bunch_particles[j].IsHadron() && pdfbase==NULL)
	THROW(critical_error,"PDF '"+set+"' does not exist in any of the loaded"
              +" libraries for "+ToString(m_bunch_particles[j])+" bunch.");
      if (pdfbase && i==0) {
	msg_Info()<<"PDF set '"<<set<<"' loaded for beam "<<j+1<<" ("
		  <<m_bunch_particles[j]<<")."<<std::endl;
      } else if (pdfbase && specializedformi) {
	msg_Info()<<"PDF set '"<<set<<"' loaded for beam "<<j+1<<" ("
		  <<m_bunch_particles[j]<<") for multiple interactions."<<std::endl;
      }
      if (pdfbase==NULL) isrbases[j] = new Intact(m_bunch_particles[j]);     
      else {
	pdfbase->SetBounds();
	isrbases[j] = new Structure_Function(pdfbase,m_bunch_particles[j]);
      }
      ATOOLS::rpa->gen.SetBunch(m_bunch_particles[j],j);
    }
    m_bunch_splimits[0] = s["ISR_SMIN"].Get<double>();
    m_bunch_splimits[1] = s["ISR_SMAX"].Get<double>();
    m_isrhandlers[id] = new ISR_Handler(isrbases);
    m_isrhandlers[id]->SetBeam(p_beamspectra->GetBeam(0),0);
    m_isrhandlers[id]->SetBeam(p_beamspectra->GetBeam(1),1);
    m_isrhandlers[id]->Init(m_bunch_splimits);
    if (!(p_beamspectra->CheckConsistency(m_bunch_particles))) {
      msg_Error()<<"Error in Environment::InitializeThePDFs()"<<endl
		 <<"   Inconsistent ISR & Beam:"<<endl
		 <<"   Abort program."<<endl;
      Abort();
    }
  }
  msg_Info() << "Initialized the ISR." << endl;
  return 1;
}

bool Initialization_Handler::InitializeTheRemnants() {
  isr::id id=isr::hard_process;
  p_remnants = new Remnant_Handler(m_isrhandlers[id],p_beamspectra);
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
  if (p_mehandler) delete p_mehandler;
  p_mehandler = new Matrix_Element_Handler(p_model);
  p_mehandler->SetShowerHandler(m_showerhandlers[isr::hard_process]);
  p_mehandler->SetRemnantHandler(p_remnants);
  auto ret = p_mehandler->InitializeProcesses(p_beamspectra,
                                              m_isrhandlers[isr::hard_process]);
  msg_Info()<<"Initialized the Matrix_Element_Handler for the hard processes."
            <<endl;
  return ret==1;
}

bool Initialization_Handler::InitializeTheUnderlyingEvents()
{
  as->SetActiveAs(isr::hard_subprocess);
  p_mihandler = new MI_Handler(p_model,
			       m_isrhandlers[isr::hard_subprocess]);
  p_mihandler->SetShowerHandler(m_showerhandlers[isr::hard_subprocess]);
  p_mihandler->SetRemnantHandler(p_remnants);
  as->SetActiveAs(isr::hard_process);
  if (p_mihandler->Type()!=0)
    msg_Info()<<"Initialized the Multiple_Interactions_Handler (MI_Handler)."<<endl;
  return true;
}

bool Initialization_Handler::InitializeTheShowers()
{
  std::vector<isr::id> isrtypes;
  isrtypes.push_back(isr::hard_process);
  isrtypes.push_back(isr::hard_subprocess);
  for (size_t i=0; i<isrtypes.size(); ++i) {
    as->SetActiveAs(isrtypes[i]);
    Shower_Handler_Map::iterator it=m_showerhandlers.find(isrtypes[i]);
    if (it!=m_showerhandlers.end()) delete it->second;
    m_showerhandlers[isrtypes[i]] =
      new Shower_Handler(p_model, m_isrhandlers[isrtypes[i]], i);
    m_showerhandlers[isrtypes[i]]->SetRemnants(p_remnants);
  }
  as->SetActiveAs(isr::hard_process);
  msg_Info()<<"Initialized the Shower_Handler."<<endl;
  return 1;
}


bool Initialization_Handler::InitializeTheSoftCollisions() 
{
  if (p_softcollisions) { delete p_softcollisions; p_softcollisions = NULL; }
  p_softcollisions = new Soft_Collision_Handler(p_beamspectra,
                                                m_isrhandlers[isr::hard_process]);
  msg_Info()<<"Initialized the Soft_Collision_Handler."<<endl;
  return 1;
}

bool Initialization_Handler::InitializeTheBeamRemnants() 
{
  if (p_beamremnants)  delete p_beamremnants;
  p_beamremnants = 
    new Beam_Remnant_Handler(p_beamspectra,
			     p_remnants,
			     p_softcollisions);
  msg_Info()<<"Initialized the Beam_Remnant_Handler."<<endl;
  return 1;
}

bool Initialization_Handler::InitializeTheFragmentation() 
{
  if (p_fragmentation) { delete p_fragmentation; p_fragmentation = NULL; }
  as->SetActiveAs(isr::hard_subprocess);
  const auto shower = m_showerhandlers[isr::hard_process]->ShowerGenerator();
  p_fragmentation = new Fragmentation_Handler(shower);
  as->SetActiveAs(isr::hard_process);
  msg_Info()<<"Initialized the Fragmentation_Handler."<<endl;
  return 1;
}

bool Initialization_Handler::InitializeTheHadronDecays() 
{
  Settings& s = Settings::GetMainSettings();
  if (s["FRAGMENTATION"].Get<std::string>() == "None")
    return true;
  std::string decmodel{ s["DECAYMODEL"].Get<std::string>() };
  msg_Tracking()<<"Decaymodel = "<<decmodel<<std::endl;
  if (decmodel=="None") return true;
  else if (decmodel==std::string("Hadrons")) {
    as->SetActiveAs(isr::hard_subprocess);
    Hadron_Decay_Handler* hd=new Hadron_Decay_Handler();
    as->SetActiveAs(isr::hard_process);
    p_hdhandler=hd;
  }
  else if ((decmodel==string("Lund")) ) {
#ifdef USING__PYTHIA
    as->SetActiveAs(isr::hard_subprocess);
    Lund_Interface * lund(NULL);
    if (p_fragmentation->GetLundInterface()==NULL) {
      lund = new Lund_Interface();
    }
    else lund = p_fragmentation->GetLundInterface();
    Lund_Decay_Handler* hd=new Lund_Decay_Handler(lund);
    as->SetActiveAs(isr::hard_process);
    p_hdhandler=hd;
#else
    THROW(fatal_error, string("Pythia not enabled during compilation. ")+
          "Use the configure option --enable-pythia to enable it.");
#endif
  }
  else {
    THROW(fatal_error,"Hadron decay model '"+decmodel+"' not implemented.");
  }
  msg_Info()<<"Initialized the Hadron_Decay_Handler, Decay model = "
            <<decmodel<<endl;
  return true;
}

bool Initialization_Handler::InitializeTheSoftPhotons()
{
  if (p_softphotons) { delete p_softphotons; p_softphotons = NULL; }
  p_softphotons = new Soft_Photon_Handler(p_mehandler);
  if (p_harddecays) p_harddecays->SetSoftPhotonHandler(p_softphotons);
  if (p_hdhandler)  p_hdhandler->SetSoftPhotonHandler(p_softphotons);
  msg_Info()<<"Initialized the Soft_Photon_Handler."<<endl;
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
        THROW(missing_module,"Cannot load Analysis library (--enable-analysis).");
    if (analyses[i]=="Rivet" || analyses[i]=="RivetME" || analyses[i]=="RivetShower") {
      if (!s_loader->LoadLibrary("SherpaHepMCOutput")) 
        THROW(missing_module,"Cannot load HepMC library (--enable-hepmc2).");
      if (!s_loader->LoadLibrary("SherpaRivetAnalysis")) 
        THROW(missing_module,"Cannot load RivetAnalysis library (--enable-rivet).");
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

bool Initialization_Handler::InitializeTheReweighting()
{
  if (p_variations) {
    delete p_variations;
  }
  Variations::CheckConsistencyWithBeamSpectra(p_beamspectra);
  p_variations = new Variations();
  msg_Info()<<"Initialized the Reweighting."<<endl;
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
  int ok = p_mehandler->CalculateTotalXSecs();
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
  rpa->gen.SetVariable("FACTORIZATION_SCALE_FACTOR", ToString(fsf));
  rpa->gen.SetVariable("RENORMALIZATION_SCALE_FACTOR", ToString(rsf));
  msg_Debugging()<<METHOD<<"(): Set scale factors {\n"
		 <<"  fac scale: "<<rpa->gen.Variable("FACTORIZATION_SCALE_FACTOR")<<"\n"
		 <<"  ren scale: "<<rpa->gen.Variable("RENORMALIZATION_SCALE_FACTOR")<<"\n}\n";

  // TODO: remove from rpa?
  double virtfrac = s["VIRTUAL_EVALUATION_FRACTION"].Get<double>();
  rpa->gen.SetVariable("VIRTUAL_EVALUATION_FRACTION", ToString(virtfrac));

  std::string evtm{ s["EVENT_GENERATION_MODE"].Get<std::string>() };
  int eventmode{ 0 };
  if (evtm=="Unweighted" || evtm=="U") eventmode=1;
  else if (evtm=="PartiallyUnweighted" || evtm=="P") eventmode=2;
  rpa->gen.SetVariable("EVENT_GENERATION_MODE",ToString(eventmode));
}
