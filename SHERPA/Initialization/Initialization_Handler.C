#include <time.h>
#include "SHERPA/Initialization/Initialization_Handler.H"

#include "SHERPA/Tools/Input_Output_Handler.H"
#include "MODEL/Main/Model_Base.H"
#include "PDF/Main/Structure_Function.H"
#include "PDF/Main/Intact.H"
#include "PDF/Main/PDF_Base.H"
#ifdef USING__Amisic
#include "AMISIC++/Main/MI_Base.H"
#endif
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Scaling.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Math/Variable.H"
#ifdef USING__PYTHIA
#include "SHERPA/LundTools/Lund_Interface.H"
#endif
#include "ATOOLS/Org/Data_Writer.H"
#include "SHERPA/Single_Events/Hadron_Decays.H"
#include "ATOOLS/Org/Library_Loader.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Selectors/Selector.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PDF/Main/Shower_Base.H"
#include "ATOOLS/Org/MyStrStream.H"

#include "ATOOLS/Phys/Spin_Correlation_Tensor.H"

#ifdef USING__Hadrons
#include "HADRONS++/Main/Hadrons.H"
#endif

#include <sys/stat.h>

using namespace SHERPA;
using namespace MODEL;
using namespace BEAM;
using namespace PDF;
using namespace ATOOLS;
using namespace std;

typedef void (*PDF_Init_Function)(const std::string &);
typedef void (*PDF_Exit_Function)();

#define PTP long unsigned int

Initialization_Handler::Initialization_Handler(int argc,char * argv[]) : 
  m_mode(0), m_savestatus(false), p_model(NULL), p_beamspectra(NULL), 
  p_harddecays(NULL), p_showerhandler(NULL), p_beamremnants(NULL), 
  p_fragmentation(NULL), p_mihandler(NULL), p_softphotons(NULL),
  p_iohandler(NULL), p_pythia(NULL), p_evtreader(NULL)
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
    p_dataread->AddWordSeparator("\t");
    p_dataread->SetInputPath(m_path);
    p_dataread->SetInputFile(m_file);
    m_analysisdat = p_dataread->GetValue<string>("ANALYSIS_DATA_FILE",string("Analysis.dat"));
    rpa.Init(m_path,m_file,argc,argv);
    return;
  }  

  SetFileNames();

  rpa.Init(m_path,m_file,argc,argv);
  LoadLibraries();
  ShowParameterSyntax();
  ran.InitExternal(m_path,m_file);

  m_spincorrelations = bool(p_dataread->GetValue<int>("SPIN_CORRELATIONS",0));
  rpa.gen.SetSpinCorrelation(m_spincorrelations);
  exh->AddTerminatorObject(this);
}

void Initialization_Handler::SetFileNames()
{
  p_dataread    = new Data_Reader(" ",";","!","=");
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
  if (cf.RereadInFile()) m_fragmentationdat=
    m_hadrondecaysdat=fname+"|(fragmentation){|}(fragmentation)";
  cf.ClearFileBegin(); cf.ClearFileEnd();
  cf.SetInputFile(fname+"|(analysis){|}(analysis)");
  if (cf.RereadInFile()) m_analysisdat=fname+"|(analysis){|}(analysis)";

  rpa.gen.SetVariable("MODEL_DATA_FILE",m_modeldat);
  rpa.gen.SetVariable("ME_DATA_FILE",m_medat);
  rpa.gen.SetVariable("MODEL_DATA_FILE",m_modeldat);
  rpa.gen.SetVariable("SHOWER_DATA_FILE",m_showerdat);
  rpa.gen.SetVariable("INTEGRATION_DATA_FILE",integrationdat);
}


Initialization_Handler::~Initialization_Handler()
{
  if (m_savestatus) {
    msg_Error()<<METHOD<<"(): Status saved to '"
	       <<rpa.gen.Variable("SHERPA_STATUS_PATH")<<"'."<<std::endl;
    MakeDir(rpa.gen.Variable("SHERPA_STATUS_PATH"),493);
    exh->PrepareTerminate();
  }
  if (p_evtreader)     { delete p_evtreader;     p_evtreader     = NULL; }
  if (p_iohandler)     { delete p_iohandler;     p_iohandler     = NULL; }
  if (p_fragmentation) { delete p_fragmentation; p_fragmentation = NULL; }
  if (p_beamremnants)  { delete p_beamremnants;  p_beamremnants  = NULL; }
  if (p_showerhandler) { delete p_showerhandler; p_showerhandler = NULL; }
  if (p_harddecays)    { delete p_harddecays;    p_harddecays    = NULL; }
  if (p_softphotons)   { delete p_softphotons;   p_softphotons   = NULL; } 
  if (p_mihandler)     { delete p_mihandler;     p_mihandler     = NULL; }
  if (p_beamspectra)   { delete p_beamspectra;   p_beamspectra   = NULL; }
  if (p_model)         { delete p_model;         p_model         = NULL; }
  if (p_pythia)        { delete p_pythia;        p_pythia        = NULL; }
  if (p_dataread)      { delete p_dataread;      p_dataread      = NULL; }
  while (m_analyses.size()>0) {
    delete m_analyses.begin()->second;
    m_analyses.erase(m_analyses.begin());
  }
  std::set<Matrix_Element_Handler*> deletedme;
  while (m_mehandlers.size()>0) {
    if (deletedme.find(m_mehandlers.begin()->second)==deletedme.end()) {
      deletedme.insert(m_mehandlers.begin()->second);
      delete m_mehandlers.begin()->second;
    }
    m_mehandlers.erase(m_mehandlers.begin());
  }
  std::set<Hadron_Decay_Handler*> deletedhd;
  while (m_hdhandlers.size()>0) {
    if (deletedhd.find(m_hdhandlers.begin()->second)==deletedhd.end()) {
      deletedhd.insert(m_hdhandlers.begin()->second);
      delete m_hdhandlers.begin()->second;
    }
    m_hdhandlers.erase(m_hdhandlers.begin());
  }
  while (m_isrhandlers.size()>0) {
    delete m_isrhandlers.begin()->second;
    m_isrhandlers.erase(m_isrhandlers.begin());
  }
  PHASIC::Phase_Space_Handler::DeleteInfo();
  exh->RemoveTerminatorObject(this);
  void *exit(s_loader->GetLibraryFunction(m_pdflib,"ExitPDFLib"));
  if (exit==NULL) THROW(fatal_error,"Cannot unload PDF library.");
  ((PDF_Exit_Function)(PTP)exit)();
}

void Initialization_Handler::LoadLibraries() const
{
  Data_Reader read(" ",";","!","=");
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
  if (!read.ReadFromFile(helpi,"SHOW_SHOWER_GENERATORS")) helpi=0;
  if (helpi>0) {
    msg->SetLevel(2);
    PDF::Shower_Base::ShowSyntax(helpi);
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
    for (Analysis_Map::iterator it=m_analyses.begin(); it!=m_analyses.end(); ++it)
      it->second->ShowSyntax(helpi);
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
  std::string path(rpa.gen.Variable("SHERPA_STATUS_PATH")+"/");
  if (path=="/") return;
  CopyFile(m_path+StripSectionTags(m_file),path+StripSectionTags(m_file));
  CopyFile(m_path+StripSectionTags(m_modeldat),path+StripSectionTags(m_modeldat));
  CopyFile(m_path+StripSectionTags(m_beamdat),path+StripSectionTags(m_beamdat));
  CopyFile(m_path+StripSectionTags(m_isrdat[0]),path+StripSectionTags(m_isrdat[0]));
  CopyFile(m_path+StripSectionTags(m_isrdat[1]),path+StripSectionTags(m_isrdat[1]));
  CopyFile(m_path+StripSectionTags(m_medat),path+StripSectionTags(m_medat));
  CopyFile(m_path+StripSectionTags(m_midat),path+StripSectionTags(m_midat));
  CopyFile(m_path+StripSectionTags(m_showerdat),path+StripSectionTags(m_showerdat));
  CopyFile(m_path+StripSectionTags(m_beamremnantdat),path+StripSectionTags(m_beamremnantdat));
  CopyFile(m_path+StripSectionTags(m_fragmentationdat),path+StripSectionTags(m_fragmentationdat));
  CopyFile(m_path+StripSectionTags(m_hadrondecaysdat),path+StripSectionTags(m_hadrondecaysdat));
  CopyFile(m_path+StripSectionTags(m_analysisdat),path+StripSectionTags(m_analysisdat));
  CopyFile(m_path+StripSectionTags(m_selectordat),
	   path+StripSectionTags(m_selectordat));
  CopyFile(m_path+StripSectionTags(m_processesdat),
	   path+StripSectionTags(m_processesdat));
  CopyFile(m_path+StripSectionTags(rpa.gen.Variable("INTEGRATION_DATA_FILE")),
	   path+StripSectionTags(rpa.gen.Variable("INTEGRATION_DATA_FILE")));
  Data_Writer writer;
  writer.SetOutputFile(path+"cmd");
  writer.SetVectorType(vtc::vertical);
  writer.AddCommandLine("SHERPA_RUN_PATH = "+
			rpa.gen.Variable("SHERPA_RUN_PATH"));
  writer.AddCommandLine("SHERPA_CPP_PATH = "+
			rpa.gen.Variable("SHERPA_CPP_PATH"));
  writer.AddCommandLine("SHERPA_LIB_PATH = "+
			rpa.gen.Variable("SHERPA_LIB_PATH"));
  writer.VectorToFile(writer.CommandLine());
}

bool Initialization_Handler::InitializeTheFramework(int nr)
{
  bool okay = true;
  SetScaleFactors();
  okay = okay && InitializeTheModel();  
  
  if (m_mode==9999) {
    msg_Events()<<"SHERPA will read in the events."<<std::endl
		<<"   The full framework is not needed."<<std::endl;
    InitializeTheAnalyses();
    return true;
  }
  okay = okay && InitializeTheBeams();
  okay = okay && InitializeThePDFs();
  okay = okay && InitializeTheAnalyses();
  ATOOLS::Integration_Info *info=PHASIC::Phase_Space_Handler::GetInfo();
  m_isrhandlers[isr::hard_process]->AssignKeys(info);
  if (m_isrhandlers.find(isr::hard_subprocess)!=m_isrhandlers.end()) {
    m_isrhandlers[isr::hard_subprocess]->AssignKeys(info);
  }
  if (!CheckBeamISRConsistency()) return 0.;
  if (m_mode>8999) {
    okay &= InitializeTheExternalMC();
    return true;
  }
  okay = okay && InitializeTheBeamRemnants();
  okay = okay && InitializeTheHardDecays();
  okay = okay && InitializeTheShowers();
  okay = okay && InitializeTheFragmentation();
  okay = okay && InitializeTheMatrixElements();
  //  only if events:
  if (rpa.gen.NumberOfEvents()>0) {
    okay = okay && InitializeTheHadronDecays();
    okay = okay && InitializeTheUnderlyingEvents();
    okay = okay && InitializeTheSoftPhotons();
    okay = okay && InitializeTheIO();
  }
  return okay;
}

bool Initialization_Handler::CheckBeamISRConsistency()
{
  if (p_model->Name()==std::string("ADD")) {
    double ms = p_model->ScalarConstant("M_s");
    if (ms<rpa.gen.Ecms()) {
      msg_Error()<<"WARNING in Initialization_Handler::CheckBeamISRConsistency :"<<std::endl
	       <<"   You are using the ADD model beyond its valid range ! "<<endl;
    }
  }

  double smin=0;
  double smax=sqr(rpa.gen.Ecms());
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
  return true;
}

bool Initialization_Handler::InitializeTheExternalMC()
{
  std::string file;
  switch (m_mode) {
#ifdef USING__PYTHIA
  case 9000: 
    p_pythia  = new Lund_Interface(m_path,m_evtfile,false);
    return true;
#endif
  default: 
    m_mode = 9999;
    msg_Info()<<"Initialization_Handler::InitializeTheExternalMC :"<<std::endl
	      <<"   SHERPA will read in the events, the full framework is not needed."<<std::endl;
  }
  return false;
}

bool Initialization_Handler::InitializeTheModel()
{
  if (p_model) delete p_model;
  //determine and set scale for coupling initialization
  Data_Reader beamer(" ",";","!","=");
  beamer.AddWordSeparator("\t");
  beamer.SetInputFile(m_path+m_beamdat);
  double beam1 = beamer.GetValue<double>("BEAM_ENERGY_1",0.0);
  double beam2 = beamer.GetValue<double>("BEAM_ENERGY_2",0.0);
  rpa.gen.SetCplScale(4.*beam1*beam2);
  Data_Reader read(" ",";","!","=");
  read.AddWordSeparator("\t");
  read.SetInputPath(m_path);
  read.SetInputFile(m_modeldat);
  std::string name;
  if (!read.ReadFromFile(name,"MODEL")) name="SM";
  p_model=Model_Base::Model_Getter_Function::
    GetObject(name,Model_Arguments(m_path,m_modeldat,true));
  if (p_model==NULL) THROW(not_implemented,"Model not implemented");
  p_model->InitializeInteractionModel();
  MODEL::s_model=p_model;
  return 1;
}


bool Initialization_Handler::InitializeTheBeams()
{
  if (p_beamspectra) { delete p_beamspectra; p_beamspectra = NULL; }
  Data_Reader dataread(" ",";","!","=");
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
  dataread.AddWordSeparator("\t");
  dataread.SetInputPath(m_path);
  dataread.SetInputFile(m_isrdat[0]);
#ifdef USING__LHAPDF
  std::string defaultlib("LHAPDFSherpa");
#else
  std::string defaultlib("CTEQ6Sherpa");
#endif  
  if (rpa.gen.Beam1().IsLepton() &&
      rpa.gen.Beam2().IsLepton()) defaultlib="PDFESherpa";
  m_pdflib=dataread.GetValue<std::string>("PDF_LIBRARY", defaultlib);
  void *init(s_loader->GetLibraryFunction(m_pdflib,"InitPDFLib"));
  if (init==NULL) THROW(fatal_error,"Cannot load PDF library.");
  std::string defset, defpath;
  if (m_pdflib=="LHAPDFSherpa") {
    defset="cteq6l.LHpdf";
    defpath="PDFSets";
  }
  else if (m_pdflib=="CTEQ6Sherpa") {
    defset="cteq6l";
    defpath="CTEQ6Grid";
  }
  else if (m_pdflib=="MRST04QEDSherpa") {
    defset="MRST04QED";
    defpath="MRST04Grid";
  }
  std::string grid_path=dataread.GetValue<string>("PDF_GRID_PATH",defpath);
  if (grid_path.length()==0 || grid_path[0]!='/')
    grid_path=rpa.gen.Variable("SHERPA_SHARE_PATH")+"/"+grid_path;
  ((PDF_Init_Function)(PTP)init)(grid_path);
  int helpi(0);
  if (!dataread.ReadFromFile(helpi,"SHOW_PDF_SETS")) helpi=0;
  if (helpi>0) {
    msg->SetLevel(2);
    PDF::PDF_Base::ShowSyntax(helpi);
    THROW(normal_exit,"Syntax shown.");
  }
  for (size_t i=0;i<2;++i) {
    isr::id id=(isr::id)(i+1);
    if (m_isrhandlers.find(id)!=m_isrhandlers.end()) 
      delete m_isrhandlers[id]; 
    dataread.SetInputFile(m_isrdat[i]);
    PDF_Base * pdfbase;
    ISR_Base ** isrbases = new ISR_Base*[2];
    double m_bunch_splimits[2];
    for (int j=0;j<2;++j) {
      int defaultflav(0);
      if (j==0) {
	defaultflav=rpa.gen.Beam1().IsAnti()? 
	  -rpa.gen.Beam1().Kfcode() : rpa.gen.Beam1().Kfcode();
      }
      else if (j==1) {
	defaultflav=rpa.gen.Beam2().IsAnti()? 
	  -rpa.gen.Beam2().Kfcode() : rpa.gen.Beam2().Kfcode();
      }
      int flav = dataread.GetValue<int>("BUNCH_"+ToString(j+1),defaultflav);
      m_bunch_particles[j] = Flavour((kf_code)abs(flav));
      if (flav<0) m_bunch_particles[j] = m_bunch_particles[j].Bar();
      std::string set = dataread.GetValue<std::string>("PDF_SET",defset);
      pdfbase = PDF_Base::PDF_Getter_Function::GetObject
	(set,PDF_Arguments(m_bunch_particles[j],grid_path,&dataread));
      if (m_bunch_particles[j].IsHadron() && pdfbase==NULL)
	THROW(critical_error,"PDF '"+set+"' does not exist in 'lib"+m_pdflib
	      +"' for "+ToString(m_bunch_particles[j])+" bunch.");
      if (i==0) msg_Info()<<"PDF set '"<<set<<"' loaded from 'lib"
			  <<m_pdflib<<"'."<<std::endl;
      if (pdfbase==NULL) isrbases[j] = new Intact(m_bunch_particles[j]);     
      else isrbases[j] = new Structure_Function(pdfbase,m_bunch_particles[j]);
      ATOOLS::rpa.gen.SetBunch(m_bunch_particles[j],j);
    }
    m_bunch_splimits[0] = dataread.GetValue<double>("ISR_SMIN",1e-10);
    m_bunch_splimits[1] = dataread.GetValue<double>("ISR_SMAX",1.);
    double kplimits[2];
    kplimits[0] = dataread.GetValue<double>("ISR_KPMIN",m_bunch_splimits[0]);
    kplimits[1] = dataread.GetValue<double>("ISR_KPMAX",m_bunch_splimits[1]);
    m_isrhandlers[id] = new ISR_Handler(isrbases);
    m_isrhandlers[id]->SetBeam(p_beamspectra->GetBeam(0),0);
    m_isrhandlers[id]->SetBeam(p_beamspectra->GetBeam(1),1);
    m_isrhandlers[id]->Init(m_bunch_splimits,kplimits);
    if (i==0)
      msg_Info()<<"Initialized the ISR: "<<m_isrhandlers[id]->Type()<<endl;
    if (!(p_beamspectra->CheckConsistency(m_bunch_particles))) {
      msg_Error()<<"Error in Environment::InitializeThePDFs()"<<endl
		 <<"   Inconsistent ISR & Beam:"<<endl
		 <<"   Abort program."<<endl;
      abort();
    }
  }
#ifdef USING__PYTHIA
  Lund_Interface::SetISRHandler(m_isrhandlers[isr::hard_process]);
#endif
  return 1;
}

bool Initialization_Handler::InitializeTheHardDecays()
{
  if (p_harddecays)    { delete p_harddecays;    p_harddecays    = NULL; }
  p_harddecays = new Hard_Decay_Handler(m_path,m_medat);
  p_harddecays->InitializeDecayMap();
  return 1;
}

bool Initialization_Handler::InitializeTheMatrixElements()
{
  Matrix_Element_Handler * me = NULL;
  me = new Matrix_Element_Handler(m_path,m_medat,m_processesdat,m_selectordat);
  me->SetShowerHandler(p_showerhandler);
  me->InitializeProcesses(p_model,p_beamspectra,m_isrhandlers[isr::hard_process]);
//   me->SetSpinCorrelations(m_spincorrelations);
  MEHandlersMap::iterator it=m_mehandlers.find("SignalMEs");
  if (it!=m_mehandlers.end()) delete it->second;
  m_mehandlers["SignalMEs"]=me; 
  msg_Info()<<"Initialized the Matrix_Element_Handler for the hard processes."
            <<endl;
  return 1;
}

Matrix_Element_Handler * const Initialization_Handler::GetMatrixElementHandler(std::string _key) { 
  MEHandlerIter pos = m_mehandlers.find(_key);
  if (pos!=m_mehandlers.end()) return pos->second;
  msg_Error()<<"Error in Initialization_Handler::GetMatrixElementHandler("<<_key<<") :"
		     <<"   Key not found. Return Null pointer."<<endl;
  return NULL;
}


bool Initialization_Handler::InitializeTheUnderlyingEvents()
{
  p_mihandler = new MI_Handler(m_path,m_midat,p_model,p_beamspectra,
			       m_isrhandlers[isr::hard_subprocess]);
  if (p_mihandler->Type()!=0)
    msg_Info()<<"Initialized the Multiple_Interactions_Handler (MI_Handler)."<<endl;
  if (p_mihandler->Name()=="None") {
    ISR_Handler_Map::iterator iit=m_isrhandlers.find(isr::hard_subprocess);
    delete iit->second;
    m_isrhandlers.erase(iit);
  }
  return true;
}

bool Initialization_Handler::InitializeTheShowers()
{
  if (p_showerhandler) delete p_showerhandler;
  p_showerhandler = new Shower_Handler(m_path,m_showerdat,p_model,
				       m_isrhandlers[isr::hard_process]);
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
  p_fragmentation = new Fragmentation_Handler(m_path,m_fragmentationdat);
  msg_Info()<<"Initialized the Fragmentation_Handler."<<endl;
  return 1;
}

bool Initialization_Handler::InitializeTheHadronDecays() 
{
  Data_Reader dr(" ",";","!","=");
  dr.AddWordSeparator("\t");
  dr.SetInputPath(m_path);
  dr.SetInputFile(m_hadrondecaysdat);
  std::string frag=dr.GetValue<string>("FRAGMENTATION",string("Ahadic"));
  if (frag=="Off") return true;

  if (m_hdhandlers.size()>0) {
    for (HDHandlersIter hditer=m_hdhandlers.begin();hditer!=m_hdhandlers.end();hditer++) {
      msg_Info()<<"Delete Hadron_Decay_Handler for "<<hditer->first<<endl;
      if (hditer->second!=NULL) { delete hditer->second; hditer->second=NULL; }
    }
    m_hdhandlers.clear();
  }

  double max_propertime = dr.GetValue<double>("MAX_PROPER_LIFETIME",-1.0);
  if( max_propertime > 0.0) {
    for(KFCode_ParticleInfo_Map::const_iterator kfit(s_kftable.begin());
	kfit!=s_kftable.end();++kfit) {
      Flavour flav(kfit->first);
      if (flav.IsOn() && flav.IsHadron() && !flav.IsStable() &&
          0.197e-12>max_propertime*flav.Width() && flav.Kfcode()!=kf_K)
      {
        flav.SetStable(true);
      }
    }
  }
  Hadron_Decays::SetMassSmearing(dr.GetValue<int>("MASS_SMEARING",1));
  
  set<kf_code>* hadrons_cans=NULL;
  Hadron_Decay_Handler * hdhandler = NULL;
  string decmodel = dr.GetValue<string>("DECAYMODEL",string("Hadrons"));
  msg_Tracking()<<"Decaymodel = "<<decmodel<<std::endl;
  if (decmodel=="Off") return true;
  if (decmodel==std::string("Hadrons")) {
    string decaypath       = dr.GetValue<string>("DECAYPATH",string("Decaydata/"));
    string decayfile       = dr.GetValue<string>("DECAYFILE",string("HadronDecays.dat"));
    string decayconstfile  = dr.GetValue<string>("DECAYCONSTFILE",string("HadronConstants.dat"));
    HADRONS::Hadrons* hadrons = new HADRONS::Hadrons(decaypath,decayfile,decayconstfile);
    hadrons->SetSpinCorrelations(m_spincorrelations);
    hdhandler              = new Hadron_Decay_Handler(hadrons);
    hadrons_cans = hdhandler->GetCans();
    m_hdhandlers["Hadrons"] = hdhandler;
  }
  if ((decmodel==string("Lund")) ) {
#ifdef USING__PYTHIA
    Lund_Interface * lund(NULL);
    if (p_fragmentation->GetLundInterface()==NULL) {
      string lfile = dr.GetValue<std::string>("LUND_FILE",std::string("Lund.dat"));
      lund         = new Lund_Interface(m_path,lfile,true);
    }
    else lund      = p_fragmentation->GetLundInterface();
    if(hadrons_cans) {
      for(set<kf_code>::iterator cankf=hadrons_cans->begin();cankf!=hadrons_cans->end();cankf++) {
        lund->SwitchOffDecays((*cankf));
      }
    }
    hdhandler      = new Hadron_Decay_Handler(lund);
    m_hdhandlers["Lund"]   = hdhandler;
#else
    THROW(fatal_error, string("Pythia not enabled during compilation. ")+
          "Use the configure option --enable-pythia to enable it.");
#endif
  }
  if (decmodel!=std::string("Hadrons") && decmodel!=string("Lund")) {
    THROW(fatal_error,"Hadron decay model not implemented.");
  }
  msg_Info()<<"Initialized the Hadron_Decay_Handler, Decay model = "<<decmodel<<endl;
  return true;
}

Hadron_Decay_Handler * const Initialization_Handler::GetHadronDecayHandler(std::string _key) { 
  HDHandlersIter pos = m_hdhandlers.find(_key);
  if (pos!=m_hdhandlers.end()) return pos->second;
  msg_Error()<<"Error in Initialization_Handler::GetHadronDecayHandler("<<_key<<") :"
	     <<"   Key not found. Return Null pointer."<<endl;
  return NULL;
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
    if (analyses[i]=="Rivet")
      if (!s_loader->LoadLibrary("SherpaRivetAnalysis")) 
        THROW(missing_module,"Cannot load RivetAnalysis library (--enable-rivet).");
    Analysis_Interface* ana=Analysis_Interface::Analysis_Getter_Function::GetObject
                            (analyses[i],Analysis_Arguments(m_path,m_analysisdat,outpath));
    if (ana==NULL) THROW(fatal_error,"Cannot initialize Analysis "+analyses[i]);
    m_analyses.insert(make_pair(analyses[i],ana));
  }
  return true;
}

bool Initialization_Handler::CalculateTheHardProcesses()
{
  if (m_mode>8999) {
    switch (m_mode) {
    case 9000:
      msg_Out()<<"SHERPA will generate the events through Pythia."<<std::endl
	       <<"   No cross sections for hard processes to be calculated."<<std::endl;
      return true;
    case 9999:
      msg_Out()<<"SHERPA will read in the events."<<std::endl
	       <<"   No cross sections for hard processes to be calculated."<<std::endl;
      return true;
    }
  }
  Matrix_Element_Handler * me = GetMatrixElementHandler(std::string("SignalMEs"));
  msg_Events()<<"=========================================================================="<<std::endl
              <<"Start calculating the hard cross sections. This may take some time.       "<<std::endl;
  int ok = me->CalculateTotalXSecs();
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

void Initialization_Handler::SetScaleFactors() 
{
  Data_Reader dr(" ",";","!","=");
  dr.AddWordSeparator("\t");
  dr.SetInputPath(m_path);
  dr.SetInputFile(m_medat);
  double sf(dr.GetValue<double>("SCALE_FACTOR",1.));
  rpa.gen.SetVariable("FACTORIZATION_SCALE_FACTOR",
		      ToString(sf*dr.GetValue<double>("FACTORIZATION_SCALE_FACTOR",1.0)));
  rpa.gen.SetVariable("RENORMALIZATION_SCALE_FACTOR",
		      ToString(sf*dr.GetValue<double>("RENORMALIZATION_SCALE_FACTOR",1.0)));
  msg_Debugging()<<METHOD<<"(): Set scale factors {\n"
		 <<"  fac scale: "<<rpa.gen.Variable("FACTORIZATION_SCALE_FACTOR")<<"\n"
		 <<"  ren scale: "<<rpa.gen.Variable("RENORMALIZATION_SCALE_FACTOR")<<"\n}\n";
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
	rpa.gen.SetVariable
	  ("SHERPA_STATUS_PATH",rpa.gen.Variable("SHERPA_RUN_PATH")+"/"+value);
	m_savestatus=true;
        oit=helpsv.erase(oit);
      }
      else if (key=="PYTHIA") {
	m_mode       = 9000;
	m_evtfile    = value;
	msg_Out()<<" Sherpa will produce Pythia events according to "<<value<<endl;
        oit=helpsv.erase(oit);
      }
      else if (key=="EVTDATA") {
	m_mode       = 9999;
	m_evtfile    = value;
	msg_Out()<<" Sherpa will read in events from : "<<value<<endl;
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
    else if (ExtractValArg(helpsv,oit,"-a","ANALYSIS"));
    else if (ExtractValArg(helpsv,oit,"-g","GENERATE_RESULT_DIRECTORY","1"));
    else if (ExtractValArg(helpsv,oit,"-b","BATCH_MODE","0"));
    else if (ExtractValArg(helpsv,oit,"-O","OUTPUT"));
    else if (par=="--version" || par=="-v"){
      msg_Out()<<"Sherpa Version "<<SHERPA_VERSION<<"."<<SHERPA_SUBVERSION<<endl;
      exit(0);
    }
    else {
      if (par!="-h" && par!="--help")
	msg_Out()<<"Unrecognized option '"<<par<<"'.\n"<<endl;
      msg_Out()<<"Usage:\n"<<endl;
      msg_Out()<<"  Sherpa [options] [<option>=<value>] [<tag>:=<value>]\n"<<endl;
      msg_Out()<<"Options:\t-f <file>      read input from file <file>"<<endl;
      msg_Out()<<"\t\t-p <path>      read input from path <path>"<<endl;
      msg_Out()<<"\t\t-e <events>    set number of events <events>"<<endl;
      msg_Out()<<"\t\t-r <results>   set result directory <results>"<<endl;
      msg_Out()<<"\t\t-a <analysis>  set analysis handler <analysis>"<<endl;
      msg_Out()<<"\t\t-O <level>     set output level <level>"<<endl;
      msg_Out()<<"\t\t-g             create result directory automatically"<<endl;
      msg_Out()<<"\t\t-b             run in non-batch mode"<<endl;
      msg_Out()<<"\t\t-v,--version   print the version number"<<endl;
      msg_Out()<<"\t\t-h,--help      print this help message\n"<<endl;
      exit(0);
    }
  }
  if (m_file.find("|")==std::string::npos) {
    Read_Write_Base cf(1,0," ",";","!","=");
    cf.SetInputPath(m_path);
    cf.SetInputFile(m_file+"|(run){|}(run)");
    if (cf.OpenInFile()) m_file+="|(run){|}(run)";
  }

  std::vector<std::string> helpsv2;
  // Add parameters from possible global.dat to command line
  Data_Reader dr(" ",";","!","=");
  dr.AddWordSeparator("\t");
  dr.AddComment("#");
  dr.SetInputPath(rpa.gen.Variable("HOME")+"/.sherpa/");
  dr.SetInputFile("global.dat");
  std::vector<std::vector<std::string> > helpsvv;
  if (dr.MatrixFromFile(helpsvv,"")) {
    msg_Out()<<METHOD<<"(): Reading parameters from '"
	     <<rpa.gen.Variable("HOME")<<"/.sherpa/global.dat'."<<std::endl;
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
  rpa.gen.SetVariable("RUN_DATA_FILE",m_file);

  if (datpath!="") m_path=datpath;
  std::vector<std::string> searchpaths;
  searchpaths.push_back(rpa.gen.Variable("SHERPA_RUN_PATH")+"/"+m_path);
  My_Out_File::SetSearchPaths(searchpaths);
  searchpaths.push_back(rpa.gen.Variable("SHERPA_DAT_PATH")+"/"+m_path);
  searchpaths.push_back(rpa.gen.Variable("SHERPA_DAT_PATH"));
  searchpaths.push_back(rpa.gen.Variable("SHERPA_SHARE_PATH")+"/"+m_path);
  searchpaths.push_back(rpa.gen.Variable("SHERPA_SHARE_PATH"));
  My_In_File::SetSearchPaths(searchpaths);
  rpa.gen.SetVariable("PATH_PIECE",m_path);
  m_path="";
  return m_mode;
}


