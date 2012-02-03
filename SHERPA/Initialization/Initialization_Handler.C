#include "SHERPA/Initialization/Initialization_Handler.H"

#include "SHERPA/PerturbativePhysics/Hard_Decay_Handler.H"
#include "SHERPA/PerturbativePhysics/Shower_Handler.H"
#include "SHERPA/SoftPhysics/Beam_Remnant_Handler.H"
#include "SHERPA/SoftPhysics/Fragmentation_Handler.H"
#include "SHERPA/PerturbativePhysics/MI_Handler.H"
#include "SHERPA/SoftPhysics/Soft_Photon_Handler.H"
#include "SHERPA/LundTools/Lund_Interface.H"
#include "SHERPA/Tools/Event_Reader_Base.H"
#include "SHERPA/Tools/Input_Output_Handler.H"
#include "MODEL/Main/Model_Base.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "PDF/Main/Structure_Function.H"
#include "PDF/Main/Intact.H"
#include "PDF/Main/PDF_Base.H"
#ifdef USING__Amisic
#include "AMISIC++/Main/MI_Base.H"
#endif
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include "ATOOLS/Math/Scaling.H"
#include "ATOOLS/Phys/Spinor.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Math/Variable.H"
#include "ATOOLS/Org/Data_Writer.H"
#include "SHERPA/Single_Events/Hadron_Decays.H"
#include "ATOOLS/Org/Library_Loader.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Selectors/Selector.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PDF/Main/NLOMC_Base.H"
#include "PDF/Main/Shower_Base.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Math/Random.H"

#include "SHERPA/Tools/Event_Reader.H"
#include "SHERPA/Tools/RootNtuple_Reader.H"

#include <sys/stat.h>
#include <time.h>

using namespace SHERPA;
using namespace MODEL;
using namespace BEAM;
using namespace PDF;
using namespace ATOOLS;
using namespace std;

typedef void (*PDF_Init_Function)();
typedef void (*PDF_Exit_Function)();

Initialization_Handler::Initialization_Handler(int argc,char * argv[]) : 
  m_mode(0), m_savestatus(false), p_model(NULL), p_beamspectra(NULL), 
  p_mehandler(NULL), p_harddecays(NULL), p_beamremnants(NULL),
  p_fragmentation(NULL), p_hdhandler(NULL), p_mihandler(NULL), p_softphotons(NULL),
  p_iohandler(NULL), p_evtreader(NULL)
{
  m_path=std::string("./");
  m_file=std::string("Run.dat");

  std::vector<std::string> names(4);
  names[0]="Decaydata";
  names[1]="Run.dat";
  My_In_File::SetNoComplains(names);

  ExtractCommandLineParameters(argc, argv);

  if (m_mode==9999) {
    ShowParameterSyntax();
    p_evtreader   = new Event_Reader(m_path,m_evtfile);
    p_dataread    = new Data_Reader(" ",";","!","=");
    p_dataread->AddComment("#");
    p_dataread->AddWordSeparator("\t");
    p_dataread->SetInputPath(m_path);
    p_dataread->SetInputFile(m_file);
    m_analysisdat = p_dataread->GetValue<string>("ANALYSIS_DATA_FILE",string("Analysis.dat"));
    rpa->Init(m_path,m_file,argc,argv);
    LoadLibraries();
    InitializeTheIO();
    return;
  }  

  SetFileNames();

  ATOOLS::s_loader->SetCheck(p_dataread->GetValue<int>("CHECK_LIBLOCK",0));

  rpa->Init(m_path,m_file,argc,argv);
  LoadLibraries();
  ShowParameterSyntax();
  ran->InitExternal(m_path,m_file);

  rpa->gen.SetSoftSC(p_dataread->GetValue<int>("SOFT_SPIN_CORRELATIONS",0));
  rpa->gen.SetHardSC(p_dataread->GetValue<int>("HARD_SPIN_CORRELATIONS",0));
  exh->AddTerminatorObject(this);
}

void Initialization_Handler::SetFileNames()
{
  p_dataread    = new Data_Reader(" ",";","!","=");
  p_dataread->AddComment("#");
  p_dataread->AddWordSeparator("\t");
  p_dataread->SetInputPath(m_path);
  p_dataread->SetInputFile(m_file);
  m_modeldat         = p_dataread->GetValue<string>("MODEL_DATA_FILE",string("Model.dat"));
  m_beamdat          = p_dataread->GetValue<string>("BEAM_DATA_FILE",string("Beam.dat"));
  m_isrdat[0]        = p_dataread->GetValue<string>("ISR_DATA_FILE",string("ISR.dat"));
  m_isrdat[1]        = p_dataread->GetValue<string>("MI_ISR_DATA_FILE",m_isrdat[0]);
  m_medat            = p_dataread->GetValue<string>("ME_DATA_FILE",string("ME.dat"));
  m_midat            = p_dataread->GetValue<string>("MI_DATA_FILE",string("MI.dat"));
  m_showerdat        = p_dataread->GetValue<string>("SHOWER_DATA_FILE",string("Shower.dat"));
  m_beamremnantdat   = p_dataread->GetValue<string>("BEAMREMNANT_DATA_FILE",string("Beam.dat"));
  m_fragmentationdat = p_dataread->GetValue<string>("FRAGMENTATION_DATA_FILE",string("Fragmentation.dat"));
  m_hadrondecaysdat  = p_dataread->GetValue<string>("FRAGMENTATION_DATA_FILE",string("Fragmentation.dat"));
  m_softphotonsdat   = p_dataread->GetValue<string>("SOFT_PHOTON_DATA_FILE",string("Fragmentation.dat"));
  m_analysisdat      = p_dataread->GetValue<string>("ANALYSIS_DATA_FILE",string("Analysis.dat"));
  std::string integrationdat=p_dataread->GetValue<string>
    ("INTEGRATION_DATA_FILE","Integration.dat");
  m_processesdat=p_dataread->GetValue<string>
    ("PROCESSFILE",string("Processes.dat"));
  m_selectordat=p_dataread->
    GetValue<string>("SELECTORFILE",string("Selector.dat"));

  std::string fname(m_file);
  if (fname.find("|")!=std::string::npos) 
    fname=fname.substr(0,fname.find("|"));
  Read_Write_Base cf(1,0," ",";","!","=");
  cf.SetAddCommandLine(false);
  cf.SetInputPath(m_path);
  cf.SetInputFile(fname+"|(beam){|}(beam)");
  if (cf.RereadInFile()) m_beamremnantdat=m_beamdat=fname+"|(beam){|}(beam)";
  cf.ClearFileBegin(); cf.ClearFileEnd();
  cf.SetInputFile(fname+"|(isr){|}(isr)");
  if (cf.RereadInFile()) m_isrdat[0]=m_isrdat[1]=fname+"|(isr){|}(isr)";
  cf.ClearFileBegin(); cf.ClearFileEnd();
  cf.SetInputFile(fname+"|(model){|}(model)");
  if (cf.RereadInFile()) m_modeldat=fname+"|(model){|}(model)";
  cf.ClearFileBegin(); cf.ClearFileEnd();
  cf.SetInputFile(fname+"|(me){|}(me)");
  if (cf.RereadInFile()) m_medat=fname+"|(me){|}(me)";
  cf.ClearFileBegin(); cf.ClearFileEnd();
  cf.SetInputFile(fname+"|(processes){|}(processes)");
  if (cf.RereadInFile()) m_processesdat=fname+"|(processes){|}(processes)";
  cf.ClearFileBegin(); cf.ClearFileEnd();
  cf.SetInputFile(fname+"|(selector){|}(selector)");
  if (cf.RereadInFile()) m_selectordat=fname+"|(selector){|}(selector)";
  cf.ClearFileBegin(); cf.ClearFileEnd();
  cf.SetInputFile(fname+"|(integration){|}(integration)");
  if (cf.RereadInFile()) integrationdat=fname+"|(integration){|}(integration)";
  cf.ClearFileBegin(); cf.ClearFileEnd();
  cf.SetInputFile(fname+"|(mi){|}(mi)");
  if (cf.RereadInFile()) m_midat=fname+"|(mi){|}(mi)";
  cf.ClearFileBegin(); cf.ClearFileEnd();
  cf.SetInputFile(fname+"|(shower){|}(shower)");
  if (cf.RereadInFile()) m_showerdat=fname+"|(shower){|}(shower)";
  cf.ClearFileBegin(); cf.ClearFileEnd();
  cf.SetInputFile(fname+"|(fragmentation){|}(fragmentation)");
  if (cf.RereadInFile()) m_fragmentationdat=m_hadrondecaysdat=
    m_softphotonsdat=fname+"|(fragmentation){|}(fragmentation)";
  cf.ClearFileBegin(); cf.ClearFileEnd();
  cf.SetInputFile(fname+"|(analysis){|}(analysis)");
  if (cf.RereadInFile()) m_analysisdat=fname+"|(analysis){|}(analysis)";

  rpa->gen.SetVariable("MODEL_DATA_FILE",m_modeldat);
  rpa->gen.SetVariable("ME_DATA_FILE",m_medat);
  rpa->gen.SetVariable("MODEL_DATA_FILE",m_modeldat);
  rpa->gen.SetVariable("SHOWER_DATA_FILE",m_showerdat);
  rpa->gen.SetVariable("INTEGRATION_DATA_FILE",integrationdat);
  rpa->gen.SetVariable("FRAGMENTATION_DATA_FILE",m_fragmentationdat);
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
  if (p_iohandler)     { delete p_iohandler;     p_iohandler     = NULL; }
  if (p_mehandler)     { delete p_mehandler;     p_mehandler     = NULL; }
  if (p_fragmentation) { delete p_fragmentation; p_fragmentation = NULL; }
  if (p_beamremnants)  { delete p_beamremnants;  p_beamremnants  = NULL; }
  if (p_harddecays)    { delete p_harddecays;    p_harddecays    = NULL; }
  if (p_hdhandler)     { delete p_hdhandler;     p_hdhandler     = NULL; }
  if (p_softphotons)   { delete p_softphotons;   p_softphotons   = NULL; } 
  if (p_mihandler)     { delete p_mihandler;     p_mihandler     = NULL; }
  if (p_beamspectra)   { delete p_beamspectra;   p_beamspectra   = NULL; }
  if (p_model)         { delete p_model;         p_model         = NULL; }
  if (p_dataread)      { delete p_dataread;      p_dataread      = NULL; }
  while (m_analyses.size()>0) {
    delete m_analyses.back();
    m_analyses.pop_back();
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
    if (exit==NULL) THROW(fatal_error,"Cannot unload PDF library "+*pdflib);
    ((PDF_Exit_Function)exit)();
  }
}

void Initialization_Handler::LoadLibraries() const
{
  Data_Reader read(" ",";","!","=");
  read.AddComment("#");
  read.SetInputFile(m_path+m_file);
  std::vector<std::string> ldadd;
  if (!read.VectorFromFile(ldadd,"SHERPA_LDADD")) return;
  for (size_t i(0);i<ldadd.size();++i) 
    if (!s_loader->LoadLibrary(ldadd[i])) 
      THROW(fatal_error,"Cannot load extra library.");
}

void Initialization_Handler::ShowParameterSyntax()
{
  Data_Reader read(" ",";","!","=");
  int helpi(0);
  if (!read.ReadFromFile(helpi,"SHOW_ME_GENERATORS")) helpi=0;
  if (helpi>0) {
    msg->SetLevel(2);
    PHASIC::ME_Generator_Base::ShowSyntax(helpi);
    THROW(normal_exit,"Syntax shown.");
  }
  if (!read.ReadFromFile(helpi,"SHOW_NLOMC_GENERATORS")) helpi=0;
  if (helpi>0) {
    msg->SetLevel(2);
    PDF::NLOMC_Base::ShowSyntax(helpi);
    THROW(normal_exit,"Syntax shown.");
  }
  if (!read.ReadFromFile(helpi,"SHOW_SHOWER_GENERATORS")) helpi=0;
  if (helpi>0) {
    msg->SetLevel(2);
    PDF::Shower_Base::ShowSyntax(helpi);
    THROW(normal_exit,"Syntax shown.");
  }
  if (!read.ReadFromFile(helpi,"SHOW_SCALE_SYNTAX")) helpi=0;
  if (helpi>0) {
    msg->SetLevel(2);
    PHASIC::Scale_Setter_Base::ShowSyntax(helpi);
    THROW(normal_exit,"Syntax shown.");
  }
  if (!read.ReadFromFile(helpi,"SHOW_SELECTOR_SYNTAX")) helpi=0;
  if (helpi>0) {
    msg->SetLevel(2);
    PHASIC::Selector_Base::ShowSyntax(helpi);
    THROW(normal_exit,"Syntax shown.");
  }
  if (!read.ReadFromFile(helpi,"SHOW_MODEL_SYNTAX")) helpi=0;
  if (helpi>0) {
    msg->SetLevel(2);
    MODEL::Model_Base::ShowSyntax(helpi);
    THROW(normal_exit,"Syntax shown.");
  }
  if (!read.ReadFromFile(helpi,"SHOW_ANALYSIS_SYNTAX")) helpi=0;
  if (helpi>0) {
    msg->SetLevel(2);
    InitializeTheAnalyses();
    for (Analysis_Vector::iterator it=m_analyses.begin(); it!=m_analyses.end(); ++it)
      (*it)->ShowSyntax(helpi);
    THROW(normal_exit,"Syntax shown.");
  }
  if (!read.ReadFromFile(helpi,"SHOW_VARIABLE_SYNTAX")) helpi=0;
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
  std::string path(rpa->gen.Variable("SHERPA_STATUS_PATH")+"/");
  if (path=="/") return;
  Copy(m_path+StripSectionTags(m_file),path+StripSectionTags(m_file));
  Copy(m_path+StripSectionTags(m_modeldat),path+StripSectionTags(m_modeldat));
  Copy(m_path+StripSectionTags(m_beamdat),path+StripSectionTags(m_beamdat));
  Copy(m_path+StripSectionTags(m_isrdat[0]),path+StripSectionTags(m_isrdat[0]));
  Copy(m_path+StripSectionTags(m_isrdat[1]),path+StripSectionTags(m_isrdat[1]));
  Copy(m_path+StripSectionTags(m_medat),path+StripSectionTags(m_medat));
  Copy(m_path+StripSectionTags(m_midat),path+StripSectionTags(m_midat));
  Copy(m_path+StripSectionTags(m_showerdat),path+StripSectionTags(m_showerdat));
  Copy(m_path+StripSectionTags(m_beamremnantdat),path+StripSectionTags(m_beamremnantdat));
  Copy(m_path+StripSectionTags(m_fragmentationdat),path+StripSectionTags(m_fragmentationdat));
  Copy(m_path+StripSectionTags(m_hadrondecaysdat),path+StripSectionTags(m_hadrondecaysdat));
  Copy(m_path+StripSectionTags(m_analysisdat),path+StripSectionTags(m_analysisdat));
  Copy(m_path+StripSectionTags(m_selectordat),
	   path+StripSectionTags(m_selectordat));
  Copy(m_path+StripSectionTags(m_processesdat),
	   path+StripSectionTags(m_processesdat));
  Copy(m_path+StripSectionTags(rpa->gen.Variable("INTEGRATION_DATA_FILE")),
	   path+StripSectionTags(rpa->gen.Variable("INTEGRATION_DATA_FILE")));
  Data_Writer writer;
  writer.SetOutputFile(path+"cmd");
  writer.SetVectorType(vtc::vertical);
  writer.AddCommandLine("SHERPA_RUN_PATH = "+
			rpa->gen.Variable("SHERPA_RUN_PATH"));
  writer.AddCommandLine("SHERPA_CPP_PATH = "+
			rpa->gen.Variable("SHERPA_CPP_PATH"));
  writer.AddCommandLine("SHERPA_LIB_PATH = "+
			rpa->gen.Variable("SHERPA_LIB_PATH"));
  writer.VectorToFile(writer.CommandLine());
}

bool Initialization_Handler::InitializeTheFramework(int nr)
{
  bool okay = true;
  Spinor<double>::SetDefaultGauge(1);
  Spinor<long double>::SetDefaultGauge(1);
  SetGlobalVariables();
  okay = okay && InitializeTheModel();  
  
  if (m_mode>=9000 && !(m_mode==9990 || m_mode==9991)) {
    msg_Events()<<"SHERPA will read in the events."<<std::endl
		<<"   The full framework is not needed."<<std::endl;
    InitializeTheAnalyses();
    return true;
  }
  okay = okay && InitializeTheBeams();
  okay = okay && InitializeThePDFs();
  if (!p_model->ModelInit(m_isrhandlers))
    THROW(critical_error,"Model cannot be initialized");
  okay = okay && p_beamspectra->Init();
  p_model->InitializeInteractionModel();
  okay = okay && InitializeTheAnalyses();
  if (!CheckBeamISRConsistency()) return 0.;
  if (m_mode==9990 || m_mode==9991) {
    p_evtreader = new RootNtuple_Reader
      (m_path,m_evtfile,m_mode-9990,p_model,
       m_isrhandlers[isr::hard_process]);
    msg_Events()<<"SHERPA will read in the events."<<std::endl
		<<"   The full framework is not needed."<<std::endl;
    InitializeTheIO();
    InitializeTheAnalyses();
    return true;
  }
  PHASIC::Phase_Space_Handler::GetInfo();
  okay = okay && InitializeTheBeamRemnants();
  okay = okay && InitializeTheHardDecays();
  okay = okay && InitializeTheShowers();
  okay = okay && InitializeTheFragmentation();
  okay = okay && InitializeTheMatrixElements();
  //  only if events:
  if (rpa->gen.NumberOfEvents()>0) {
    okay = okay && InitializeTheSoftPhotons();
    okay = okay && InitializeTheHadronDecays();
    okay = okay && InitializeTheUnderlyingEvents();
    okay = okay && InitializeTheIO();
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
  p_iohandler = new Input_Output_Handler(p_dataread);
  p_iohandler->SetMEHandler(p_mehandler);
  return true;
}

bool Initialization_Handler::InitializeTheModel()
{
  if (p_model) delete p_model;
  //determine and set scale for coupling initialization
  Data_Reader beamer(" ",";","!","=");
  beamer.AddComment("#");
  beamer.AddWordSeparator("\t");
  beamer.SetInputFile(m_path+m_beamdat);
  double beam1 = beamer.GetValue<double>("BEAM_ENERGY_1",0.0);
  double beam2 = beamer.GetValue<double>("BEAM_ENERGY_2",0.0);
  rpa->gen.SetCplScale(4.*beam1*beam2);
  Data_Reader read(" ",";","!","=");
  read.AddWordSeparator("\t");
  read.SetInputPath(m_path);
  read.SetInputFile(m_modeldat);
  std::string name;
  if (!read.ReadFromFile(name,"MODEL")) name="SM";
  p_model=Model_Base::Model_Getter_Function::
    GetObject(name,Model_Arguments(m_path,m_modeldat,true));
  if (p_model==NULL) THROW(not_implemented,"Model not implemented");
  MODEL::s_model=p_model;
  return 1;
}


bool Initialization_Handler::InitializeTheBeams()
{
  if (p_beamspectra) { delete p_beamspectra; p_beamspectra = NULL; }
  Data_Reader dataread(" ",";","!","=");
  dataread.AddComment("#");
  dataread.AddWordSeparator("\t");
  dataread.SetInputPath(m_path);
  dataread.SetInputFile(m_beamdat);
  p_beamspectra        = new Beam_Spectra_Handler(&dataread);
  msg_Info()<<"Initialized the beams "<<p_beamspectra->Type()<<endl;
  return 1;
}

bool Initialization_Handler::InitializeThePDFs()
{
  Data_Reader dataread(" ",";","!","=");
  dataread.AddComment("#");
  dataread.AddWordSeparator("\t");
  dataread.SetInputPath(m_path);
  dataread.SetInputFile(m_isrdat[0]);

  // load PDF libraries
  std::string defset[2];
  for (int beam(0);beam<=1;++beam) {
    std::string deflib("None");
    if (p_beamspectra->GetBeam(beam)->Bunch().Kfcode()==kf_p_plus) {
      deflib="CTEQ6Sherpa";
      defset[beam]="cteq6.6m";
    }
    else if (p_beamspectra->GetBeam(beam)->Bunch().Kfcode()==kf_e) {
      deflib="PDFESherpa";
      defset[beam]="PDFe";
    }
    else if (p_beamspectra->GetBeam(beam)->Bunch().IsPhoton()) {
      deflib="GRVSherpa";
      defset[beam]="GRV";
    }
    m_pdflibs.insert(dataread.GetValue<std::string>("PDF_LIBRARY",deflib));
    std::string mpilib, beamlib;
    if (dataread.ReadFromFile(mpilib,"PDF_LIBRARY_MPI"))
      m_pdflibs.insert(mpilib);
    if (dataread.ReadFromFile(beamlib,"PDF_LIBRARY_"+ToString(beam+1)))
      m_pdflibs.insert(beamlib);
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
  int helpi(0);
  if (!dataread.ReadFromFile(helpi,"SHOW_PDF_SETS")) helpi=0;
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
    dataread.SetInputFile(m_isrdat[i]);
    PDF_Base * pdfbase;
    ISR_Base ** isrbases = new ISR_Base*[2];
    double m_bunch_splimits[2];
    for (int j=0;j<2;++j) {
      int defaultflav(p_beamspectra->GetBeam(j)->Bunch());
      int flav = dataread.GetValue<int>("BUNCH_"+ToString(j+1),defaultflav);
      m_bunch_particles[j] = Flavour((kf_code)abs(flav));
      if (flav<0) m_bunch_particles[j] = m_bunch_particles[j].Bar();
      std::string set = dataread.GetValue<std::string>("PDF_SET",defset[j]);
      std::string specialset;
      if (dataread.ReadFromFile(specialset,"PDF_SET_"+ToString(j+1)))
	set=specialset;
      if (id==isr::hard_subprocess) {
        std::string mpiset;
        if (dataread.ReadFromFile(mpiset,"PDF_SET_MPI")) {
          set=mpiset;
        }
      }
      pdfbase = PDF_Base::PDF_Getter_Function::GetObject
	(set,PDF_Arguments(m_bunch_particles[j],&dataread, j));
      if (m_bunch_particles[j].IsHadron() && pdfbase==NULL)
	THROW(critical_error,"PDF '"+set+"' does not exist in any of the loaded"
              +" libraries for "+ToString(m_bunch_particles[j])+" bunch.");
      if (pdfbase && i==0) {
	msg_Info()<<"PDF set '"<<set<<"' loaded for beam "<<j+1<<" ("
		  <<m_bunch_particles[j]<<")."<<std::endl;
      }
      if (pdfbase==NULL) isrbases[j] = new Intact(m_bunch_particles[j]);     
      else {
	pdfbase->SetBounds();
	isrbases[j] = new Structure_Function(pdfbase,m_bunch_particles[j]);
      }
      ATOOLS::rpa->gen.SetBunch(m_bunch_particles[j],j);
    }
    m_bunch_splimits[0] = dataread.GetValue<double>("ISR_SMIN",1e-10);
    m_bunch_splimits[1] = dataread.GetValue<double>("ISR_SMAX",1.);
    m_isrhandlers[id] = new ISR_Handler(isrbases);
    m_isrhandlers[id]->SetBeam(p_beamspectra->GetBeam(0),0);
    m_isrhandlers[id]->SetBeam(p_beamspectra->GetBeam(1),1);
    m_isrhandlers[id]->Init(m_bunch_splimits);
    if (i==0)
      msg_Info()<<"Initialized the ISR: "<<m_isrhandlers[id]->Type()<<endl;
    if (!(p_beamspectra->CheckConsistency(m_bunch_particles))) {
      msg_Error()<<"Error in Environment::InitializeThePDFs()"<<endl
		 <<"   Inconsistent ISR & Beam:"<<endl
		 <<"   Abort program."<<endl;
      abort();
    }
  }
  return 1;
}

bool Initialization_Handler::InitializeTheHardDecays()
{
  Data_Reader dr(" ",";","!","=");
  dr.AddWordSeparator("\t");
  dr.SetInputPath(m_path);
  dr.SetInputFile(m_medat);
  std::string decays=dr.GetValue<string>("HARD_DECAYS",string("Off"));
  if (decays=="Off") return true;

  if (p_harddecays)    { delete p_harddecays;    p_harddecays    = NULL; }
  p_harddecays = new Hard_Decay_Handler(m_path,m_medat);
  return 1;
}

bool Initialization_Handler::InitializeTheMatrixElements()
{
  if (p_mehandler) delete p_mehandler;
  p_mehandler = new Matrix_Element_Handler(m_path,m_medat,m_processesdat,m_selectordat);
  p_mehandler->SetShowerHandler(m_showerhandlers[isr::hard_process]);
  p_mehandler->InitializeProcesses(p_model,p_beamspectra,m_isrhandlers[isr::hard_process]);
  msg_Info()<<"Initialized the Matrix_Element_Handler for the hard processes."
            <<endl;
  return 1;
}

bool Initialization_Handler::InitializeTheUnderlyingEvents()
{
  as->SetActiveAs(isr::hard_subprocess);
  p_mihandler = new MI_Handler(m_path,m_midat,p_model,p_beamspectra,
			       m_isrhandlers[isr::hard_subprocess]);
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
    m_showerhandlers[isrtypes[i]]=new Shower_Handler
        (m_path, m_showerdat, p_model, m_isrhandlers[isrtypes[i]]);
  }
  as->SetActiveAs(isr::hard_process);
  msg_Info()<<"Initialized the Shower_Handler."<<endl;
  return 1;
}


bool Initialization_Handler::InitializeTheBeamRemnants() 
{
  if (p_beamremnants)  delete p_beamremnants;
  p_beamremnants = 
    new Beam_Remnant_Handler(m_path,m_beamremnantdat,
			     m_isrhandlers[isr::hard_process],p_beamspectra);
  p_beamremnants->SetScale(4.0);
  msg_Info()<<"Initialized the Beam_Remnant_Handler."<<endl;
  return 1;
}

bool Initialization_Handler::InitializeTheFragmentation() 
{
  if (p_fragmentation) { delete p_fragmentation; p_fragmentation = NULL; }
  as->SetActiveAs(isr::hard_subprocess);
  p_fragmentation = new Fragmentation_Handler(m_path,m_fragmentationdat);
  as->SetActiveAs(isr::hard_process);
  msg_Info()<<"Initialized the Fragmentation_Handler."<<endl;
  return 1;
}

bool Initialization_Handler::InitializeTheHadronDecays() 
{
  Data_Reader dr(" ",";","!","=");
  dr.AddComment("#");
  dr.AddWordSeparator("\t");
  dr.SetInputPath(m_path);
  dr.SetInputFile(m_hadrondecaysdat);
  std::string frag=dr.GetValue<string>("FRAGMENTATION",string("Ahadic"));
  if (frag=="Off") return true;

  string decmodel = dr.GetValue<string>("DECAYMODEL",string("Hadrons"));
  msg_Tracking()<<"Decaymodel = "<<decmodel<<std::endl;
  if (decmodel=="Off") return true;
  else if (decmodel==std::string("Hadrons")) {
    as->SetActiveAs(isr::hard_subprocess);
    Hadron_Decay_Handler* hd=new Hadron_Decay_Handler(m_path,m_hadrondecaysdat);
    hd->SetSoftPhotonHandler(p_softphotons);
    as->SetActiveAs(isr::hard_process);
    p_hdhandler=hd;
  }
  else if ((decmodel==string("Lund")) ) {
#ifdef USING__PYTHIA
    THROW(not_implemented, "Lund hadron decays not implemented.");
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
  p_softphotons = new Soft_Photon_Handler(m_path,m_softphotonsdat);
  msg_Info()<<"Initialized the Soft_Photon_Handler."<<endl;
  return true;
}

bool Initialization_Handler::InitializeTheAnalyses()
{
  std::string outpath=p_dataread->GetValue<std::string>("ANALYSIS_OUTPUT","Analysis/");
  std::string analysis=p_dataread->GetValue<std::string>("ANALYSIS","0");
  std::vector<std::string> analyses;
  Data_Reader readline(",",";","#","");
  readline.SetString(analysis);
  readline.VectorFromString(analyses);
  for (size_t i=0; i<analyses.size(); ++i) {
    if (analyses[i]=="0") continue;
    if (analyses[i]=="1") analyses[i]="Internal";
    if (analyses[i]=="Internal")
      if (!s_loader->LoadLibrary("SherpaAnalysis")) 
        THROW(missing_module,"Cannot load Analysis library (--enable-analysis).");
    if (analyses[i]=="Rivet" || analyses[i]=="RivetME" || analyses[i]=="RivetShower")
      if (!s_loader->LoadLibrary("SherpaRivetAnalysis")) 
        THROW(missing_module,"Cannot load RivetAnalysis library (--enable-rivet).");
    Analysis_Interface* ana=Analysis_Interface::Analysis_Getter_Function::GetObject
                            (analyses[i],Analysis_Arguments(m_path,m_analysisdat,outpath));
    if (ana==NULL) THROW(fatal_error,"Cannot initialize Analysis "+analyses[i]);
    m_analyses.push_back(ana);
  }
  return true;
}

bool Initialization_Handler::CalculateTheHardProcesses()
{
  if (m_mode>9000) return true;
  msg_Events()<<"=========================================================================="<<std::endl
              <<"Start calculating the hard cross sections. This may take some time.       "<<std::endl;
  ATOOLS::Data_Reader read(" ",";","!","=");
  ATOOLS::msg->SetLevel(read.GetValue<int>("INT_OUTPUT",ATOOLS::msg->Level()));
  as->SetActiveAs(isr::hard_process);
  int ok = p_mehandler->CalculateTotalXSecs();
  if (ok) {
    msg_Events()<<"Calculating the hard cross sections has been successful.                  "<<std::endl
	     <<"=========================================================================="<<std::endl;
  }
  else {
    msg_Events()<<"Calculating the hard cross sections failed. Check this carefully.         "<<std::endl
	     <<"=========================================================================="<<std::endl;
  }
  return ok;
}

void Initialization_Handler::SetGlobalVariables() 
{
  Data_Reader dr(" ",";","!","=");
  dr.AddComment("#");
  dr.AddWordSeparator("\t");
  dr.SetInputPath(m_path);
  dr.SetInputFile(m_medat);
  double sf(dr.GetValue<double>("SCALE_FACTOR",1.));
  rpa->gen.SetVariable("FACTORIZATION_SCALE_FACTOR",
		      ToString(sf*dr.GetValue<double>("FACTORIZATION_SCALE_FACTOR",1.0)));
  rpa->gen.SetVariable("RENORMALIZATION_SCALE_FACTOR",
		      ToString(sf*dr.GetValue<double>("RENORMALIZATION_SCALE_FACTOR",1.0)));
  msg_Debugging()<<METHOD<<"(): Set scale factors {\n"
		 <<"  fac scale: "<<rpa->gen.Variable("FACTORIZATION_SCALE_FACTOR")<<"\n"
		 <<"  ren scale: "<<rpa->gen.Variable("RENORMALIZATION_SCALE_FACTOR")<<"\n}\n";
  int cmode=dr.GetValue<int>("METS_CLUSTER_MODE",16);
  rpa->gen.SetVariable("METS_CLUSTER_MODE",ToString(cmode));
  if (cmode!=0) msg_Info()<<METHOD<<"(): Set cluster mode "<<cmode<<".\n";
}

bool Initialization_Handler::ExtractValArg
(std::vector<std::string> &args,std::vector<std::string>::iterator &it,
 const std::string &arg,const std::string &tag,const std::string &def) const
{
  if (it->find(arg)!=0) return false;
  std::string val=it->substr(2);
  if (def!="") {
    if (val!="") *it="-"+val;
    else it=args.erase(it);
    it=args.insert(it,tag+"="+def);
    return true;
  }
  if (val=="") {
    if (it+1!=args.end()) {
      it=args.erase(it);
      val=*it;
    }
  }
  if (val=="") {
    msg_Error()<<METHOD<<"(): No argument to '"
	       <<arg<<"'. Abort."<<std::endl;
    exit(1);
  }
  *it=tag+"="+val;
  return true;
}

int Initialization_Handler::ExtractCommandLineParameters(int argc,char * argv[])
{
  std::string datpath;
  std::vector<std::string> helpsv(argc-1);
  for (int i(0);i<argc-1;++i) helpsv[i]=argv[i+1];
  for (std::vector<std::string>::iterator oit(helpsv.begin());
       oit!=helpsv.end();) {
    string par = *oit;
    string key,value;
    size_t equal=par.find("=");
    if (equal!=std::string::npos) {
      value = par.substr(equal+1);
      key   = par = par.substr(0,equal);
      if (key=="PATH") {
	if (value[value.length()-1]=='/') value.erase(value.length()-1,1);
	m_path=value;
        oit=helpsv.erase(oit);
      }
      else if (key=="RUNDATA") {
	m_file=value;
        oit=helpsv.erase(oit);
      }
      else if (key=="STATUS_PATH") {
	if (value[value.length()-1]!='/') value+=std::string("/");
	datpath=value;
        oit=helpsv.erase(oit);
      }
      else if (key=="SAVE_STATUS") {
	if (value[value.length()-1]!='/') value+=std::string("/");
	rpa->gen.SetVariable
	  ("SHERPA_STATUS_PATH",rpa->gen.Variable("SHERPA_RUN_PATH")+"/"+value);
	m_savestatus=true;
        oit=helpsv.erase(oit);
      }
      else if (key=="EVTDATA") {
	m_mode       = 9999;
	m_evtfile    = value;
	msg_Out()<<" Sherpa will read in events from : "<<value<<endl;
        oit=helpsv.erase(oit);
      }
      else if (key=="EVTDATA_ROOTNTUPLE") {
	m_mode       = 9990;
	m_evtfile    = value;
	msg_Out()<<" Sherpa will read in root ntuple events from : "<<value<<endl;
        Read_Write_Base::AddCommandLine("ROOTNTUPLE_OUTPUT=");
        oit=helpsv.erase(oit);
      }
      else if (key=="EVTDATA_ROOTNTUPLE_RC") {
	m_mode       = 9991;
	m_evtfile    = value;
	msg_Out()<<" Sherpa will read in root ntuple events from : "<<value<<endl;
        Read_Write_Base::AddCommandLine("ROOTNTUPLE_OUTPUT=");
        oit=helpsv.erase(oit);
      }
      else {
	++oit;
      }
    }
    else if (ExtractValArg(helpsv,oit,"-f","RUNDATA"));
    else if (ExtractValArg(helpsv,oit,"-p","PATH"));
    else if (ExtractValArg(helpsv,oit,"-e","EVENTS"));
    else if (ExtractValArg(helpsv,oit,"-r","RESULT_DIRECTORY"));
    else if (ExtractValArg(helpsv,oit,"-L","SHERPA_CPP_PATH"));
    else if (ExtractValArg(helpsv,oit,"-R","RANDOM_SEED"));
    else if (ExtractValArg(helpsv,oit,"-m","ME_SIGNAL_GENERATOR"));
    else if (ExtractValArg(helpsv,oit,"-w","EVENT_GENERATION_MODE"));
    else if (ExtractValArg(helpsv,oit,"-s","SHOWER_GENERATOR"));
    else if (ExtractValArg(helpsv,oit,"-F","FRAGMENTATION"));
    else if (ExtractValArg(helpsv,oit,"-D","DECAYMODEL"));
    else if (ExtractValArg(helpsv,oit,"-a","ANALYSIS"));
    else if (ExtractValArg(helpsv,oit,"-A","ANALYSIS_OUTPUT"));
    else if (ExtractValArg(helpsv,oit,"-b","BATCH_MODE","0"));
    else if (ExtractValArg(helpsv,oit,"-O","OUTPUT"));
    else if (ExtractValArg(helpsv,oit,"-o","EVT_OUTPUT"));
    else if (ExtractValArg(helpsv,oit,"-l","LOG_FILE"));
    else if (ExtractValArg(helpsv,oit,"-j","PG_THREADS"));
    else if (ExtractValArg(helpsv,oit,"-g","GENERATE_RESULT_DIRECTORY","0"));
    else if (par=="--version" || par=="-v"){
      msg_Out()<<"Sherpa Version "<<SHERPA_VERSION<<"."<<SHERPA_SUBVERSION<<endl;
      exit(0);
    }
    else {
      if (par!="-h" && par!="--help")
	msg_Out()<<"Unrecognized option '"<<par<<"'.\n"<<endl;
      msg_Out()<<"Usage:\n"<<endl;
      msg_Out()<<"  Sherpa [options] [<parameter>=<value>] [<tag>:=<value>]\n"<<endl;
      msg_Out()<<"Options:\t-f <file>         read input from file <file>"<<endl;
      msg_Out()<<"\t\t-p <path>         read input from path <path>"<<endl;
      msg_Out()<<"\t\t-e <events>       set number of events <events>"<<endl;
      msg_Out()<<"\t\t-r <results>      set result directory <results>"<<endl;
      msg_Out()<<"\t\t-m <generator>    set me generator <generator>"<<endl;
      msg_Out()<<"\t\t-w <mode>         set event generation mode <mode>"<<endl;
      msg_Out()<<"\t\t-s <generator>    set ps generator <generator>"<<endl;
      msg_Out()<<"\t\t-F <module>       set fragmentation module <module>"<<endl;
      msg_Out()<<"\t\t-D <module>       set decay module <module>"<<endl;
      msg_Out()<<"\t\t-a <analysis>     set analysis handler <analysis>"<<endl;
      msg_Out()<<"\t\t-A <path>         set analysis output path <path>"<<endl;
      msg_Out()<<"\t\t-O <level>        set general output level <level>"<<endl;
      msg_Out()<<"\t\t-o <level>        set output level for event generation"<<endl;
      msg_Out()<<"\t\t-l <logfile>      set log file name <logfile>"<<endl;
      msg_Out()<<"\t\t-j <threads>      set number of threads <threads>"<<endl;
      msg_Out()<<"\t\t-g                do not create result directory"<<endl;
      msg_Out()<<"\t\t-b                run in non-batch mode"<<endl;
      msg_Out()<<"\t\t-v,--version      print the version number"<<endl;
      msg_Out()<<"\t\t-h,--help         print this help message\n"<<endl;
      exit(0);
    }
  }
  if (m_file.find("|")==std::string::npos) {
    Read_Write_Base cf(1,0," ",";","!","=");
    cf.SetInputPath(m_path);
    cf.SetInputFile(m_file+"|(run){|}(run)");
    if (cf.OpenInFile()) m_file+="|(run){|}(run)";
  }

  if (datpath!="") m_path=datpath;
  std::vector<std::string> searchpaths;
  searchpaths.push_back(rpa->gen.Variable("SHERPA_RUN_PATH")+"/"+m_path);
  My_Out_File::SetSearchPaths(searchpaths);
  searchpaths.push_back(rpa->gen.Variable("SHERPA_DAT_PATH")+"/"+m_path);
  searchpaths.push_back(rpa->gen.Variable("SHERPA_DAT_PATH"));
  searchpaths.push_back(rpa->gen.Variable("SHERPA_SHARE_PATH")+"/"+m_path);
  searchpaths.push_back(rpa->gen.Variable("SHERPA_SHARE_PATH"));
  My_In_File::SetSearchPaths(searchpaths);
  rpa->gen.SetVariable("PATH_PIECE",m_path);
  m_path="";

  std::vector<std::string> helpsv2;
  // Add parameters from possible global.dat to command line
  Data_Reader dr(" ",";","!","=");
  dr.AddWordSeparator("\t");
  dr.AddComment("#");
  dr.SetInputPath(rpa->gen.Variable("HOME")+"/.sherpa/");
  dr.SetInputFile("global.dat");
  std::vector<std::vector<std::string> > helpsvv;
  if (dr.MatrixFromFile(helpsvv,"")) {
    msg_Out()<<METHOD<<"(): Reading parameters from '"
	     <<rpa->gen.Variable("HOME")<<"/.sherpa/global.dat'."<<std::endl;
    helpsv2.resize(helpsvv.size());
    for (size_t i(0);i<helpsvv.size();++i) {
      helpsv2[i]=helpsvv[i][0];
      for (size_t j(1);j<helpsvv[i].size();++j) helpsv2[i]+=" "+helpsvv[i][j];
    }
  }
  // Add parameters from Run.dat to command line
  // (this makes it possible to overwrite particle properties in Run.dat)
  dr.SetInputPath(m_path);
  dr.SetInputFile(m_file);
  dr.RereadInFile();
  if (dr.MatrixFromFile(helpsvv,"")) {
    size_t oldsize(helpsv2.size());
    helpsv2.resize(oldsize+helpsvv.size());
    for (size_t i(0);i<helpsvv.size();++i) {
      helpsv2[oldsize+i]=helpsvv[i][0];
      for (size_t j(1);j<helpsvv[i].size();++j)
	helpsv2[oldsize+i]+=" "+helpsvv[i][j];
    }
  }
  helpsv2.insert(helpsv2.end(),helpsv.begin(),helpsv.end());
  for (size_t i(0);i<helpsv2.size();++i) {
    string par = helpsv2[i];
    string key,value;
    size_t equal=Min(par.find("="),par.find(" "));
    if (equal!=std::string::npos) {
      value = par.substr(equal+1);
      key   = par = par.substr(0,equal);
      if (key[key.length()-1]==':') {
        key.erase(key.length()-1,1);
        Read_Write_Base::AddGlobalTag(key,value);
      }
      else {
        Read_Write_Base::AddCommandLine(key+" = "+value+"; ");
      }
    }
    else {
      Read_Write_Base::AddCommandLine(par+";");
    }
  }
  rpa->gen.SetVariable("RUN_DATA_FILE",m_file);

  return m_mode;
}


