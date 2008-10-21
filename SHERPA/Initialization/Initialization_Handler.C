#include <time.h>
#include "Initialization_Handler.H"

#include "Input_Output_Handler.H"
#include "Model_Base.H"
#include "Structure_Function.H"
#include "Intact.H"
#include "PDF_Handler.H"
#include "PDF_Base.H"
#include "Initial_State_Shower.H"
#include "MI_Base.H"
#include "Data_Reader.H"
#include "Message.H"
#include "Scaling.H"
#include "Shell_Tools.H"
#include "Particle_Qualifier.H"
#include "Variable.H"
#include "CXXFLAGS.H"
#ifdef USING__PYTHIA
#include "Lund_Interface.H"
#endif
#include "Data_Writer.H"
#include "Hadron_Decays.H"
#include "Library_Loader.H"

#include "Spin_Correlation_Tensor.H"

#ifdef USING__Hadrons
#include "Hadrons.H"
#endif

#include <sys/stat.h>

using namespace SHERPA;
using namespace MODEL;
using namespace BEAM;
using namespace PDF;
using namespace ATOOLS;
using namespace std;

Initialization_Handler::Initialization_Handler(int argc,char * argv[]) : 
  m_mode(0), m_savestatus(false), p_model(NULL), p_beamspectra(NULL), 
  p_harddecays(NULL), p_showerhandler(NULL), p_beamremnants(NULL), 
  p_fragmentation(NULL), p_mihandler(NULL), p_softphotons(NULL),
  p_iohandler(NULL), p_pythia(NULL), p_evtreader(NULL), 
  p_analysis(NULL)
{
  m_path=std::string("./");
  m_file=std::string("Run.dat");

  std::vector<std::string> names(4);
  names[0]="Decaydata";
  names[1]="Run.dat";
  My_In_File::SetNoComplains(names);

  ExtractCommandLineParameters(argc, argv);
  ShowParameterSyntax();

  if (m_mode==9999) {
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
  ran.InitExternal(m_path,m_file);

  CheckFlagConsistency();
  
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
  m_decaydat         = p_dataread->GetValue<string>("DECAY_DATA_FILE",string("Decays.dat"));
  m_showerdat        = p_dataread->GetValue<string>("SHOWER_DATA_FILE",string("Shower.dat"));
  m_beamremnantdat   = p_dataread->GetValue<string>("BEAMREMNANT_DATA_FILE",string("Beam.dat"));
  m_fragmentationdat = p_dataread->GetValue<string>("FRAGMENTATION_DATA_FILE",string("Fragmentation.dat"));
  m_hadrondecaysdat  = p_dataread->GetValue<string>("FRAGMENTATION_DATA_FILE",string("Fragmentation.dat"));
  m_softphotonsdat   = p_dataread->GetValue<string>("SOFT_PHOTON_DATA_FILE",string("Fragmentation.dat"));
  m_analysisdat      = p_dataread->GetValue<string>("ANALYSIS_DATA_FILE",string("Analysis.dat"));
  std::string integrationdat=p_dataread->GetValue<string>
    ("INTEGRATION_DATA_FILE","Integration.dat");
  std::string processdat=p_dataread->GetValue<string>
    ("PROCESSFILE",string("Processes.dat"));
  std::string selectordat=p_dataread->
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
  if (cf.RereadInFile()) processdat=fname+"|(processes){|}(processes)";
  cf.ClearFileBegin(); cf.ClearFileEnd();
  cf.SetInputFile(fname+"|(selector){|}(selector)");
  if (cf.RereadInFile()) selectordat=fname+"|(selector){|}(selector)";
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

  rpa.gen.SetVariable("MODEL_DATA_FILE",m_modeldat);
  rpa.gen.SetVariable("ME_DATA_FILE",m_medat);
  rpa.gen.SetVariable("SHOWER_DATA_FILE",m_showerdat);
  rpa.gen.SetVariable("INTEGRATION_DATA_FILE",integrationdat);
  rpa.gen.SetVariable("PROCESSFILE",processdat);
  rpa.gen.SetVariable("SELECTORFILE",selectordat);
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
  if (p_analysis)      { delete p_analysis;      p_analysis      = NULL; }
  if (p_dataread)      { delete p_dataread;      p_dataread      = NULL; }
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
}

void Initialization_Handler::LoadLibraries() const
{
  Data_Reader read(" ",";","!","=");
  std::vector<std::string> ldadd;
  if (!read.VectorFromFile(ldadd,"SHERPA_LDADD")) return;
  for (size_t i(0);i<ldadd.size();++i) 
    if (!s_loader->LoadLibrary(ldadd[i])) 
      THROW(fatal_error,"Cannot load extra library.");
}

void Initialization_Handler::ShowParameterSyntax() const
{
  Data_Reader read(" ",";","!","=");
  int helpi(0);
  if (!read.ReadFromFile(helpi,"SHOW_MODEL_SYNTAX")) helpi=0;
  if (helpi>0) {
    msg->SetLevel(2);
    MODEL::Model_Base::ShowSyntax(helpi);
    THROW(normal_exit,"Syntax shown.");
  }
  if (!read.ReadFromFile(helpi,"SHOW_ANALYSIS_SYNTAX")) helpi=0;
  if (helpi>0) {
    msg->SetLevel(2);
    ATOOLS::Variable_Base<double>::ShowVariables(helpi);
    ATOOLS::Particle_Qualifier_Base::ShowQualifiers(helpi);
    ANALYSIS::Analysis_Handler::ShowSyntax(helpi);
    THROW(normal_exit,"Syntax shown.");
  }
  if (!read.ReadFromFile(helpi,"SHOW_QUALIFIER_SYNTAX")) helpi=0;
  if (helpi>0) {
    msg->SetLevel(2);
    ATOOLS::Particle_Qualifier_Base::ShowQualifiers(helpi);
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
  CopyFile(m_path+StripSectionTags(m_decaydat),path+StripSectionTags(m_decaydat));
  CopyFile(m_path+StripSectionTags(m_showerdat),path+StripSectionTags(m_showerdat));
  CopyFile(m_path+StripSectionTags(m_beamremnantdat),path+StripSectionTags(m_beamremnantdat));
  CopyFile(m_path+StripSectionTags(m_fragmentationdat),path+StripSectionTags(m_fragmentationdat));
  CopyFile(m_path+StripSectionTags(m_hadrondecaysdat),path+StripSectionTags(m_hadrondecaysdat));
  CopyFile(m_path+StripSectionTags(m_analysisdat),path+StripSectionTags(m_analysisdat));
  CopyFile(m_path+StripSectionTags(rpa.gen.Variable("SELECTORFILE")),
	   path+StripSectionTags(rpa.gen.Variable("SELECTORFILE")));
  CopyFile(m_path+StripSectionTags(rpa.gen.Variable("PROCESSFILE")),
	   path+StripSectionTags(rpa.gen.Variable("PROCESSFILE")));
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
  okay = okay && InitializeTheMatrixElements();
  //  only if events:
  if (rpa.gen.NumberOfEvents()>0) {
    okay = okay && InitializeTheShowers();
    okay = okay && InitializeTheFragmentation();
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
  p_model->FillDecayTables();
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
  for (size_t i=0;i<2;++i) {
    isr::id id=(isr::id)(i+1);
    if (m_isrhandlers.find(id)!=m_isrhandlers.end()) 
      delete m_isrhandlers[id]; 
    Data_Reader dataread(" ",";","!","=");
    dataread.AddWordSeparator("\t");
    dataread.SetInputPath(m_path);
    dataread.SetInputFile(m_isrdat[i]);
    PDF_Handler pdfhandler;
    PDF_Base * pdfbase;
    ISR_Base ** isrbases = new ISR_Base*[2];
    double m_bunch_splimits[2];
    for (int j=0;j<2;++j) {
      pdfbase = pdfhandler.GetPDFLib(&dataread,m_bunch_particles[j],j);
      if (m_bunch_particles[j].IsHadron() && pdfbase==NULL)
	THROW(fatal_error,"ISR must be enabled in 'ISR.dat' for "
	      +ToString(m_bunch_particles[j])+" bunch.");
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
    msg_Info()<<"Initialized the ISR["<<id<<"] : "<<m_isrhandlers[id]->Type()<<endl;
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
  p_harddecays = new Hard_Decay_Handler(m_path,m_decaydat,m_medat,p_model);
  if (p_harddecays->GetMEHandler()!=NULL) {
    msg_Info()<<"Initialized the Hard_Decay_Handler. Its ME_Handler is : "
		      <<p_harddecays->GetMEHandler()->Name()<<"/"
		      <<p_harddecays->GetMEHandler()<<std::endl;
    m_mehandlers.insert(std::make_pair(std::string("HardDecays"),p_harddecays->GetMEHandler()));
  }
  return 1;
}

bool Initialization_Handler::InitializeTheMatrixElements()
{
  Matrix_Element_Handler * me = NULL;
  if (p_harddecays) {
    me = new Matrix_Element_Handler(m_path,m_medat,p_model,p_beamspectra,
				    m_isrhandlers[isr::hard_process],
				    p_harddecays->GetMEHandler());
  }
  else {
    me = new Matrix_Element_Handler(m_path,m_medat,p_model,p_beamspectra,
				    m_isrhandlers[isr::hard_process],NULL);
  }
  me->SetSpinCorrelations(m_spincorrelations);
  MEHandlersMap::iterator it=m_mehandlers.find("SignalMEs");
  if (it!=m_mehandlers.end()) delete it->second;
  m_mehandlers["SignalMEs"]=me; 
  msg_Info()<<"Initialized the Matrix_Element_Handler for the hard processes : "<<me->Name()<<endl;
  if (p_analysis) {
    int weighted=1-me->EventGenerationMode();
    msg_Info()<<"Initialization_Handler::InitializeTheMatrixElements(): "
			  <<"Setting analysis mode "<<(weighted?"weighted":"unweighted")<<endl;
    p_analysis->SetWeighted(weighted);
  }
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
  Matrix_Element_Handler *mehandler = p_mihandler->HardMEHandler();
  if (mehandler!=NULL) {
    m_mehandlers.insert(std::make_pair(std::string("MIMEs"),mehandler)); 
    msg_Info()<<"Added the Matrix_Element_Handler for the u.e. :"<<mehandler->Name()<<endl;
  }
  else {
    ISR_Handler_Map::iterator iit=m_isrhandlers.find(isr::hard_subprocess);
    delete iit->second;
    m_isrhandlers.erase(iit);
  }
  return true;
}

bool Initialization_Handler::InitializeTheShowers()
{
  if (p_showerhandler) delete p_showerhandler;
  int maxjets     = GetMatrixElementHandler(std::string("SignalMEs"))->MaxJets();
  p_showerhandler = new Shower_Handler(m_path,m_showerdat,p_model,
				       m_isrhandlers[isr::hard_process],maxjets);
  APACIC::Apacic *apacic=p_showerhandler->GetApacic();
  if (apacic!=NULL && apacic->IniShower()!=NULL) 
    p_beamremnants->SetScale(-apacic->IniShower()->CutOff());
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
  int test = p_dataread->GetValue<int>("TEST_DETECTOR",0);

  int helpi=p_dataread->GetValue<int>("ANALYSIS",0);
  if (!helpi&&test==0) return true;
  std::string outpath=p_dataread->GetValue<std::string>("ANALYSIS_OUTPUT","Analysis/");
  p_analysis = new ANALYSIS::Analysis_Handler();
  p_analysis->SetInputPath(m_path);
  p_analysis->SetInputFile(m_analysisdat);
  p_analysis->SetOutputPath(outpath);
  if (test>0) {
    p_analysis->Test(test);
    THROW(normal_exit,"Tested detector.");
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
  int scalechoice = 0;
  if (p_showerhandler) {
    if (p_showerhandler->ISROn()) scalechoice += 1;
    if (p_showerhandler->FSROn()) scalechoice += 2;
  }
  Matrix_Element_Handler * me = GetMatrixElementHandler(std::string("SignalMEs"));
  msg_Events()<<"=========================================================================="<<std::endl
              <<"Start calculating the hard cross sections. This may take some time.       "<<std::endl;
  int ok = me->CalculateTotalXSecs(scalechoice);
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
  if (rpa.gen.NumberOfEvents()==0 || 
      rpa.gen.Variable("SUDAKOV_WEIGHT","1")!="1") return;
  Data_Reader reader(" ",";","!","=");
  reader.AddWordSeparator("\t");
  reader.SetInputPath(rpa.gen.Variable("SHERPA_DAT_PATH")+m_path+"/");
  reader.SetInputFile(m_showerdat);
  bool changed(false);
  double fac(1.0);
  if (!reader.ReadFromFile(fac,"IS_CPL_SCALE_FACTOR")) fac=1.0;
  else changed=true;
  rpa.gen.SetVariable("IS_CPL_SCALE_FACTOR",ToString(fac));
  if (!reader.ReadFromFile(fac,"FS_CPL_SCALE_FACTOR")) fac=1.0;
  else changed=true;
  rpa.gen.SetVariable("FS_CPL_SCALE_FACTOR",ToString(fac));
  int scheme(0);
  if (!reader.ReadFromFile(scheme,"S_KFACTOR_SCHEME")) scheme=0;
  else changed=true;
  rpa.gen.SetVariable("S_KFACTOR_SCHEME",ToString(scheme));
  if (changed)
    msg_Error()<<om::bold<<METHOD<<"(): WARNING {\n"<<om::reset<<om::red
	       <<"  Scale- and K-factors for Matrix Element weighting\n"
	       <<"  are set to account for Parton Shower settings.\n"
	       <<"  If re-using integration results, please make sure\n"
	       <<"  that the integration was performed this way.\n"
	       <<om::reset<<om::bold<<"}"<<om::reset<<std::endl;
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
    else if (par=="--version" || par=="-v"){
      msg_Out()<<" Sherpa Version "<<SHERPA_VERSION<<"."<<SHERPA_SUBVERSION<<endl;
      msg_Out()<<"   employing: "<<endl;
      msg_Out()<<"    * AMEGIC++ Version 2."<<SHERPA_SUBVERSION<<endl;
      msg_Out()<<"    * APACIC++ Version 2."<<SHERPA_SUBVERSION<<endl;
      msg_Out()<<"    * Pythia Version 6.214"<<endl;
      exit(0);
    }
    else {
      msg_Out()<<" Usage: "<<endl;
      msg_Out()<<" Sherpa [<variable>=<value>] "<<endl;
      msg_Out()<<endl;
      msg_Out()<<" Possible options: "<<endl;
      msg_Out()<<"  -v,--version   prints the Version number"<<endl;
      msg_Out()<<"  -?,--help      prints this help message"<<endl;
      exit(0);
    }
  }
  if (m_file.find("|")==std::string::npos) {
    Read_Write_Base cf(1,0," ",";","!","=");
    cf.SetInputPath(m_path);
    cf.SetInputFile(m_file+"|(run){|}(run)");
    if (cf.OpenInFile()) m_file+="|(run){|}(run)";
  }

  // Add parameters from Run.dat to command line
  // (this makes it possible to overwrite particle properties in Run.dat)
  Data_Reader dr(" ",";","!");
  dr.AddWordSeparator("\t");
  dr.AddComment("#");
  dr.SetInputPath(m_path);
  dr.SetInputFile(m_file);
  std::vector<std::vector<std::string> > helpsvv;
  dr.MatrixFromFile(helpsvv,"");
  std::vector<std::string> helpsv2(helpsvv.size());
  for (size_t i(0);i<helpsvv.size();++i) {
    for (size_t j(0);j<helpsvv[i].size();++j) helpsv2[i]+=helpsvv[i][j];
  }
  helpsv2.insert(helpsv2.end(),helpsv.begin(),helpsv.end());
  for (size_t i(0);i<helpsv2.size();++i) {
    string par = helpsv2[i];
    string key,value;
    size_t equal=par.find("=");
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


void Initialization_Handler::CheckFlagConsistency()
{
  Data_Reader dr(" ",";","!","=");
  dr.AddWordSeparator("\t");
  dr.SetInputPath(m_path);
  dr.SetInputFile(m_medat);
  int  sudweight = dr.GetValue<int>("SUDAKOV_WEIGHT",1);
  rpa.gen.SetVariable("SUDAKOV_WEIGHT",ToString(sudweight));

  // if SUDAKOV_WEIGHT=On
  if (sudweight>0) {
    //  Run.dat
    long nevt = p_dataread->GetValue<long>("EVENTS",0);
    if (nevt<=0) {
      Read_Write_Base::AddCommandLine("EVENTS = 1; ");
    }

    //  ME.dat 
    if (dr.GetValue<std::string>("SCALE_SCHEME","CKKW")!="CKKW" ||
        dr.GetValue<std::string>("KFACTOR_SCHEME","1")!="1" ||
        dr.GetValue<std::string>("COUPLING_SCHEME","Running_alpha_S")!="Running_alpha_S") {
      msg_Error()<<om::bold<<METHOD<<"(): WARNING {\n"<<om::reset<<om::red
                 <<"  CKKW is switched on by 'SUDAKOV_WEIGHT = 1'.\n"
                 <<"  This causes the following settings:\n"
                 <<"    SCALE_SCHEME = CKKW\n"
                 <<"    KFACTOR_SCHEME = 1\n"
                 <<"    COUPLING_SCHEME = Running_alpha_S\n"<<om::reset
                 <<om::bold<<"}"<<om::reset<<std::endl;
    }
    Read_Write_Base::AddCommandLine("SCALE_SCHEME = CKKW; ");
    Read_Write_Base::AddCommandLine("KFACTOR_SCHEME = 1; ");
    Read_Write_Base::AddCommandLine("COUPLING_SCHEME = Running_alpha_S; ");

    //  Shower.dat
    Read_Write_Base::AddCommandLine("FSR_SHOWER = 1; ");
  }
  else {
    Read_Write_Base::AddCommandLine("JET_VETO_SCHEME = 0; ");
    Read_Write_Base::AddCommandLine("LOSE_JET_SCHEME = 0; ");
  }


  // check if all MI.dat / Decays.dat
  

}
