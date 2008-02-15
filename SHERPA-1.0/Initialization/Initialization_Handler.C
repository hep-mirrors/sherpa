#include <time.h>
#include "Initialization_Handler.H"

#include "Input_Output_Handler.H"
#include "Model_Handler.H"
#include "Structure_Function.H"
#include "Intact.H"
#include "PDF_Handler.H"
#include "PDF_Base.H"
#include "Initial_State_Shower.H"
#include "MI_Base.H"
#include "Data_Read.H"
#include "Data_Reader.H"
#include "Message.H"
#include "Scaling.H"
#include "Shell_Tools.H"
#include "Particle_Qualifier.H"
#include "Data_Collector.H"
#include "Variable.H"
#include "Lund_Interface.H"
#include "Data_Writer.H"

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

Initialization_Handler::Initialization_Handler(string _path,string _file) : 
  m_path(_path), m_file(_file), m_mode(0), m_savestatus(false),
  p_model(NULL), p_beamspectra(NULL), p_harddecays(NULL), 
  p_showerhandler(NULL), p_beamremnants(NULL), p_fragmentation(NULL), 
  p_mihandler(NULL), p_iohandler(NULL), p_pythia(NULL), 
  p_evtreader(NULL),
  p_analysis(NULL)
{
  std::vector<std::string> names(4);
  names[0]="Decaydata";
  names[1]="Particle.dat";
  names[2]="Hadron.dat";
  names[3]="Run.dat";
  My_In_File::SetNoComplains(names);

  p_dataread         = new Data_Read(m_path+m_file);
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
  m_analysisdat      = p_dataread->GetValue<string>("ANALYSIS_DATA_FILE",string("Analysis.dat"));
  rpa.gen.SetVariable("ME_DATA_FILE",m_medat);
  rpa.gen.SetVariable("SHOWER_DATA_FILE",m_showerdat);
  m_spincorrelations = bool(p_dataread->GetValue<int>("SPIN_CORRELATIONS",0));
  rpa.gen.SetSpinCorrelation(m_spincorrelations);
  exh->AddTerminatorObject(this);
}

Initialization_Handler::Initialization_Handler(int argc,char * argv[]) : 
  m_mode(0), m_savestatus(false), p_model(NULL), p_beamspectra(NULL), 
  p_harddecays(NULL), p_showerhandler(NULL), p_beamremnants(NULL), 
  p_fragmentation(NULL), p_mihandler(NULL),
  p_iohandler(NULL), p_pythia(NULL), p_evtreader(NULL), 
  p_analysis(NULL)
{
  m_path=std::string("./");
  m_file=std::string("Run.dat");

  std::vector<std::string> names(4);
  names[0]="Decaydata";
  names[1]="Particle.dat";
  names[2]="Hadron.dat";
  names[3]="Run.dat";
  My_In_File::SetNoComplains(names);

  ExtractCommandLineParameters(argc, argv);

  if (m_mode==9999) {
    p_evtreader   = new Event_Reader(m_path,m_evtfile);
    p_dataread    = new Data_Read(m_path+m_file);
    m_analysisdat = p_dataread->GetValue<string>("ANALYSIS_DATA_FILE",string("Analysis.dat"));
    ShowParameterSyntax();
    rpa.Init(m_path,m_file,argc,argv);
    return;
  }  
  ShowParameterSyntax();
  rpa.Init(m_path,m_file,argc,argv);

  p_dataread         = new Data_Read(m_path+m_file);
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
  m_analysisdat      = p_dataread->GetValue<string>("ANALYSIS_DATA_FILE",string("Analysis.dat"));
  rpa.gen.SetVariable("ME_DATA_FILE",m_medat);
  rpa.gen.SetVariable("SHOWER_DATA_FILE",m_showerdat);

  CheckFlagConsistency();
  
  m_spincorrelations = bool(p_dataread->GetValue<int>("SPIN_CORRELATIONS",0));
  rpa.gen.SetSpinCorrelation(m_spincorrelations);
  exh->AddTerminatorObject(this);
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

void Initialization_Handler::ShowParameterSyntax() const
{
  Data_Reader read(" ",";","!","=");
  int helpi(0);
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

void Initialization_Handler::PrepareTerminate()
{
  std::string path(rpa.gen.Variable("SHERPA_STATUS_PATH")+"/");
  if (path=="/") return;
  CopyFile(m_path+m_file,path+m_file);
  CopyFile(m_path+m_modeldat,path+m_modeldat);
  CopyFile(m_path+m_beamdat,path+m_beamdat);
  CopyFile(m_path+m_isrdat[0],path+m_isrdat[0]);
  CopyFile(m_path+m_isrdat[1],path+m_isrdat[1]);
  CopyFile(m_path+m_medat,path+m_medat);
  CopyFile(m_path+m_midat,path+m_midat);
  CopyFile(m_path+m_decaydat,path+m_decaydat);
  CopyFile(m_path+m_showerdat,path+m_showerdat);
  CopyFile(m_path+m_beamremnantdat,path+m_beamremnantdat);
  CopyFile(m_path+m_fragmentationdat,path+m_fragmentationdat);
  CopyFile(m_path+m_hadrondecaysdat,path+m_hadrondecaysdat);
  CopyFile(m_path+m_analysisdat,path+m_analysisdat);
  CopyFile(m_path+rpa.gen.Variable("SELECTORFILE"),
	   path+rpa.gen.Variable("SELECTORFILE"));
  CopyFile(m_path+rpa.gen.Variable("PROCESSFILE"),
	   path+rpa.gen.Variable("PROCESSFILE"));
  CopyFile(m_path+"Integration.dat",path+"Integration.dat");
  CopyFile(m_path+"Particle.dat",path+"Particle.dat");
  CopyFile(m_path+"Hadron.dat",path+"Hadron.dat");
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
  if (nr<=0) {
    ATOOLS::ParticleInit(m_path); 
  }
  bool okay = InitializeTheIO();
  if (m_mode==9999) {
    msg_Events()<<"SHERPA will read in the events."<<std::endl
	     <<"   The full framework is not needed."<<std::endl;
    InitializeTheAnalyses();
    return true;
  }
  if (rpa.gen.NumberOfEvents()>0) SetScaleFactors();
  okay = okay && InitializeTheBeams();
  okay = okay && InitializeThePDFs();
  okay = okay && InitializeTheModel();  
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
  std::vector<std::string> infiles, outfiles;
  infiles.push_back(p_dataread->GetValue<string>("SHERPA_INPUT",string("")));
  infiles.push_back(p_dataread->GetValue<string>("HEPMC_INPUT",string("")));
  infiles.push_back(p_dataread->GetValue<string>("HEPEVT_INPUT",string("")));
  infiles.push_back(p_dataread->GetValue<string>("D0_HEPEVT_INPUT",string("")));
  outfiles.push_back(p_dataread->GetValue<string>("SHERPA_OUTPUT",string("")));
  outfiles.push_back(p_dataread->GetValue<string>("HEPMC_OUTPUT",string("")));
  outfiles.push_back(p_dataread->GetValue<string>("OLD_HEPMC_OUTPUT",string("")));
  outfiles.push_back(p_dataread->GetValue<string>("HEPEVT_OUTPUT",string("")));
  outfiles.push_back(p_dataread->GetValue<string>("D0_HEPEVT_OUTPUT",string("")));
  outfiles.push_back(p_dataread->GetValue<string>("HEPMC2_OUTPUT",string("")));
  std::string evtpath = p_dataread->GetValue<string>("EVT_FILE_PATH",m_path);
  int filesize        = p_dataread->GetValue<int>("FILE_SIZE",1000);
  int precision       = p_dataread->GetValue<int>("OUTPUT_PRECISION",6);
  std::string outmode = p_dataread->GetValue<string>("EVENT_MODE",string("Sherpa"));

  p_iohandler = new Input_Output_Handler(outmode,outfiles,infiles,evtpath,filesize,precision);

  return true;
}

bool Initialization_Handler::InitializeTheExternalMC()
{
  std::string file;
  switch (m_mode) {
  case 9000: 
    p_pythia  = new Lund_Interface(m_path,m_evtfile,false);
    return true;
  default: 
    m_mode = 9999;
    msg_Info()<<"Initialization_Handler::InitializeTheExternalMC :"<<std::endl
	      <<"   SHERPA will read in the events, the full framework is not needed."<<std::endl;
  }
  return false;
}

bool Initialization_Handler::InitializeTheModel()
{
  if (p_model) { delete p_model; p_model = NULL; }
  Data_Read     * dataread     = new Data_Read(m_path+m_modeldat);
  // - add commandline parameter - !!
  Model_Handler * modelhandler = new Model_Handler();
  p_model                      = modelhandler->GetModel(dataread,m_path,m_modeldat);

  delete modelhandler;
  delete dataread;  
  return 1;
}


bool Initialization_Handler::InitializeTheBeams()
{
  if (p_beamspectra) { delete p_beamspectra; p_beamspectra = NULL; }
  Data_Read * dataread = new Data_Read(m_path+m_beamdat);
  p_beamspectra        = new Beam_Spectra_Handler(dataread);
  msg_Info()<<"Initialized the beams "<<p_beamspectra->Type()<<endl;
  delete dataread;  
  return 1;
}


bool Initialization_Handler::InitializeThePDFs()
{
  for (size_t i=0;i<2;++i) {
    isr::id id=(isr::id)(i+1);
    if (m_isrhandlers.find(id)!=m_isrhandlers.end()) 
      delete m_isrhandlers[id]; 
    Data_Read dataread(m_path+m_isrdat[i]);
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
    m_bunch_splimits[0] = dataread.GetValue<double>("ISR_SMIN",0.);
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
  Lund_Interface::SetISRHandler(m_isrhandlers[isr::hard_process]);
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
  Data_Read dr(m_path+m_hadrondecaysdat);
  std::string frag=dr.GetValue<string>("FRAGMENTATION",string("Off"));
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
    Fl_Iter fli;
    for (Flavour flav=fli.first();flav!=Flavour(kf::none);flav = fli.next()) {
      if (flav.IsOn() && flav.IsHadron() && !flav.IsStable() &&
          0.197e-12>max_propertime*flav.Width() && flav.Kfcode()!=kf::K)
      {
        flav.SetStable(true);
      }
    }
  }
  
  bool needextra = true; set<kf::code>* hadrons_cans=NULL;
  Hadron_Decay_Handler * hdhandler = NULL;
  string decmodel = dr.GetValue<string>("DECAYMODEL",string("Lund"));
  msg_Tracking()<<"Decaymodel = "<<decmodel<<std::endl;
#ifdef USING__Hadrons
  if (decmodel==std::string("Hadrons")) {
    string decaypath       = dr.GetValue<string>("DECAYPATH",string("Decaydata/"));
    struct stat fst;
    if (stat(decaypath.c_str(),&fst)==-1 || (fst.st_mode&S_IFMT)!=S_IFDIR)
      decaypath=rpa.gen.Variable("SHERPA_SHARE_PATH")+"/Decaydata/";
    string decayfile       = dr.GetValue<string>("DECAYFILE",string("HadronDecays.dat"));
    string decayconstfile  = dr.GetValue<string>("DECAYCONSTFILE",string("HadronConstants.dat"));
    HADRONS::Hadrons* hadrons = new HADRONS::Hadrons(decaypath,decayfile,decayconstfile);
    hadrons->SetSpinCorrelations(m_spincorrelations);
    hdhandler              = new Hadron_Decay_Handler(hadrons);
    hadrons_cans = hdhandler->GetCans();
    hdhandler->SetMassSmearing(dr.GetValue<int>("MASS_SMEARING",1));
    m_hdhandlers["Hadrons"] = hdhandler;
  }
#endif
  if ((decmodel==string("Lund") || needextra) ) {
    Lund_Interface * lund(NULL);
    if (p_fragmentation->GetLundInterface()==NULL) {
      string lfile = dr.GetValue<std::string>("LUND_FILE",std::string("Lund.dat"));
      lund         = new Lund_Interface(m_path,lfile,true);
    }
    else lund      = p_fragmentation->GetLundInterface();
    if(hadrons_cans) {
      for(set<kf::code>::iterator cankf=hadrons_cans->begin();cankf!=hadrons_cans->end();cankf++) {
        lund->SwitchOffDecays((*cankf));
      }
    }
    hdhandler      = new Hadron_Decay_Handler(lund);
    hdhandler->SetMassSmearing(dr.GetValue<int>("MASS_SMEARING",1));
    m_hdhandlers["Lund"]   = hdhandler;
  }
  if (decmodel!=std::string("Hadrons") && decmodel!=string("Lund")) {
    THROW(critical_error,"Fragmentation model not implemented.");
    abort();
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


bool Initialization_Handler::InitializeTheAnalyses()
{
  int test = p_dataread->GetValue<int>("TEST_DETECTOR",0);

  int helpi=p_dataread->GetValue<int>("ANALYSIS",0);
  if (!helpi&&test==0) return true;
  std::string outpath=p_dataread->GetValue<std::string>("ANALYSIS_OUTPUT","./");
  if (outpath==NotDefined<std::string>()) outpath="";
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
  if (ok) 
    msg_Events()<<"Calculating the hard cross sections has been successful.                  "<<std::endl
	     <<"=========================================================================="<<std::endl;
  else
    msg_Events()<<"Calculating the hard cross sections failed. Check this carefully.         "<<std::endl
	     <<"=========================================================================="<<std::endl;
  return ok;
}

void Initialization_Handler::SetScaleFactors() 
{
  if (rpa.gen.Variable("SUDAKOV_WEIGHT","0")!="1") return;
  Data_Reader reader(" ",";","!","=");
  reader.AddWordSeparator("\t");
  reader.SetInputPath(m_path+"/");
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
  map<std::string,int> special_options;
  special_options["-V"]=12;
  special_options["--version"]=12;
  special_options["-?"]=13;
  special_options["-h"]=13;
  special_options["--help"]=13;

  special_options["PATH"]=101;
  special_options["RUNDATA"]=102;
  special_options["STATUS_PATH"]=110;
  special_options["SAVE_STATUS"]=111;
  special_options["PYTHIA"]=9000;
  special_options["EVTDATA"]=9999;
  
  std::string datpath;
  

  for (int i=1; i<argc;++i) {
    int mode   = 0;
    string par = string(argv[i]);
    string key,value;
    int equal  = par.find("=");
    if (equal!=-1) {
      mode=1;
      value = par.substr(equal+1);
      key   = par = par.substr(0,equal);
    }


    if (special_options.find(par)!=special_options.end()) 
      mode = special_options[par];

    // variables in dat files
    if (equal!=-1 && mode==1) {
      if (key[key.length()-1]==':') {
	key.erase(key.length()-1,1);
	Data_Read::SetGlobalTag(key,value);
	Read_Write_Base::AddGlobalTag(key,value);
      }
      else {
	Data_Read::SetCommandLine(key,value);
	Read_Write_Base::AddCommandLine(key+" = "+value+"; ");
      }
    }
    
    // special variables
    if (mode>=100) {
      MyStrStream s;
      switch (mode) {
      case 101:
	if (value[value.length()-1]=='/') value.erase(value.length()-1,1);
	m_path=value;
	break;
      case 102:
	m_file=value;
	break;
      case 110:
	if (value[value.length()-1]!='/') value+=std::string("/");
	datpath=value;
	break;
      case 111:
	if (value[value.length()-1]!='/') value+=std::string("/");
	rpa.gen.SetVariable
	  ("SHERPA_STATUS_PATH",rpa.gen.Variable("SHERPA_RUN_PATH")+"/"+value);
	m_savestatus=true;
	break;
      case 9000:
	m_mode       = 9000;
	m_evtfile    = value;
	msg_Out()<<" Sherpa will produce Pythia events according to "<<value<<endl;
	break;
      case 9001:
	m_mode       = 9001;
	m_evtfile    = value;
	msg_Out()<<" Sherpa will produce Herwig events according to "<<value<<endl;
	break;
      case 9002:
	m_mode       = 9002;
	m_evtfile    = value;
	msg_Out()<<" Sherpa will produce MCatNLO events according to "<<value<<endl;
	break;
      case 9999:
	m_mode       = 9999;
	m_evtfile    = value;
	msg_Out()<<" Sherpa will read in events from : "<<value<<endl;
	break;
      case 100:
	m_options[key] = value;
      }
    }
    else {
      // other option
      switch (mode) {
      case 12:
	{
	  // should call a version roution
	  msg_Out()<<" Sherpa Version "<<SHERPA_VERSION<<"."<<SHERPA_SUBVERSION<<endl;
	  msg_Out()<<"   employing: "<<endl;
	  msg_Out()<<"    * AMEGIC++ Version 2."<<SHERPA_SUBVERSION<<endl;
	  msg_Out()<<"    * APACIC++ Version 2."<<SHERPA_SUBVERSION<<endl;

	  string pyver("6.214");
	  msg_Out()<<"    * Pythia Version "<<pyver<<endl;
	}
	exit(0);
      case 13:
	msg_Out()<<" Help: "<<endl;
	msg_Out()<<" Sherpa [options] [<variable>=<value>] "<<endl;
	msg_Out()<<endl;
	msg_Out()<<" Possible options: "<<endl;
	msg_Out()<<"  -V,--version   prints the Version number"<<endl;
	msg_Out()<<"  -?,--help      prints this help message"<<endl;
	exit(0);
      case 15:
	// xsout
	break;
      case 16:
	// eventout
	break;
      }
    }
  }
  if (datpath!="") m_path=datpath;
  std::vector<std::string> searchpaths;
  searchpaths.push_back(rpa.gen.Variable("SHERPA_RUN_PATH")+"/"+m_path);
  My_Out_File::SetSearchPaths(searchpaths);
  searchpaths.push_back(rpa.gen.Variable("SHERPA_DAT_PATH")+"/"+m_path);
  searchpaths.push_back(rpa.gen.Variable("SHERPA_DAT_PATH"));
  searchpaths.push_back(SHERPA_SHARE_PATH+std::string("/")+m_path);
  searchpaths.push_back(SHERPA_SHARE_PATH);
  My_In_File::SetSearchPaths(searchpaths);
  rpa.gen.SetVariable("PATH_PIECE",m_path);
  m_path="";
  return m_mode;
}


void Initialization_Handler::CheckFlagConsistency()
{
  Data_Read dr(m_path+m_medat);
  int  sudweight = dr.GetValue<int>("SUDAKOV_WEIGHT",0);
  rpa.gen.SetVariable("SUDAKOV_WEIGHT",ToString(sudweight));

  // if SUDAKOV_WEIGHT=On
  if (sudweight>0) {
    //  Run.dat
    long nevt = p_dataread->GetValue<long>("EVENTS",0);
    if (nevt<=0) {
      Data_Read::SetCommandLine("EVENTS","1");
      Read_Write_Base::AddCommandLine("EVENTS = 1; ");
    }

    //  ME.dat 
    Data_Read::SetCommandLine("SCALE_SCHEME","CKKW");
    Data_Read::SetCommandLine("KFACTOR_SCHEME","1");
    Data_Read::SetCommandLine("COUPLING_SCHEME","Running_alpha_S");
    Read_Write_Base::AddCommandLine("SCALE_SCHEME = CKKW; ");
    Read_Write_Base::AddCommandLine("KFACTOR_SCHEME = 1; ");
    Read_Write_Base::AddCommandLine("COUPLING_SCHEME = Running_alpha_S; ");

    //  Shower.dat
    Data_Read::SetCommandLine("FSR_SHOWER","1");
    Read_Write_Base::AddCommandLine("FSR_SHOWER = 1; ");
  }
  else {
    Data_Read::SetCommandLine("JET_VETO_SCHEME","0");
    Data_Read::SetCommandLine("LOSE_JET_SCHEME","0");
    Read_Write_Base::AddCommandLine("JET_VETO_SCHEME = 0; ");
    Read_Write_Base::AddCommandLine("LOSE_JET_SCHEME = 0; ");
  }


  // check if all MI.dat / Decays.dat
  

}
