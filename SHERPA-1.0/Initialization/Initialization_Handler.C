#include "Initialization_Handler.H"

#include "IO_Handler.H"
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
#ifdef USING__MCatNLO
#include "MCatNLO_Wrapper.H"
#endif
#include "Data_Collector.H"
#include "Variable.H"

using namespace SHERPA;
using namespace MODEL;
using namespace BEAM;
using namespace PDF;
using namespace ATOOLS;
using namespace std;

struct Char_Array40 {
 char  s[40];
} ;

extern "C" {
  Char_Array40 visaje_();
}

Initialization_Handler::Initialization_Handler(string _path,string _file) : 
  m_path(_path), m_file(_file), m_mode(0),
  p_model(NULL), p_beamspectra(NULL), p_harddecays(NULL), 
  p_showerhandler(NULL), p_beamremnants(NULL), p_fragmentation(NULL), 
  p_hadrondecays(NULL), p_mihandler(NULL), p_iohandler(NULL), p_pythia(NULL), 
  p_evtreader(NULL),
#ifdef USING__MCatNLO
  p_herwig(NULL), p_mcatnlo(NULL),
#endif
  p_analysis(NULL) 
{
  m_scan_istep=-1;  

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
}


Initialization_Handler::Initialization_Handler(int argc,char * argv[]) : 
  m_mode(0), p_model(NULL), p_beamspectra(NULL), 
  p_harddecays(NULL), p_showerhandler(NULL), p_beamremnants(NULL), 
  p_fragmentation(NULL), p_hadrondecays(NULL), p_mihandler(NULL),
  p_iohandler(NULL), p_pythia(NULL), p_evtreader(NULL), 
#ifdef USING__MCatNLO
  p_herwig(NULL), p_mcatnlo(NULL),
#endif
  p_analysis(NULL)
{
  m_path=std::string("./");
  m_file=std::string("Run.dat");

  m_scan_istep=-1;

  ExtractCommandLineParameters(argc, argv);

  if (m_mode==9999) {
    p_evtreader   = new Event_Reader(m_path,m_evtfile);
    p_dataread    = new Data_Read(m_path+m_file);
    m_analysisdat = p_dataread->GetValue<string>("ANALYSIS_DATA_FILE",string("Analysis.dat"));
    rpa.Init(m_path,m_file,argc,argv);
    return;
  }
  
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

  CheckFlagConsitency();

  rpa.Init(m_path,m_file,argc,argv);

}


Initialization_Handler::~Initialization_Handler()
{
  if (p_evtreader)     { delete p_evtreader;     p_evtreader     = NULL; }
  if (p_iohandler)     { delete p_iohandler;     p_iohandler     = NULL; }
  if (p_hadrondecays)  { delete p_hadrondecays;  p_hadrondecays  = NULL; }
  if (p_fragmentation) { delete p_fragmentation; p_fragmentation = NULL; }
  if (p_beamremnants)  { delete p_beamremnants;  p_beamremnants  = NULL; }
  if (p_showerhandler) { delete p_showerhandler; p_showerhandler = NULL; }
  if (p_harddecays)    { delete p_harddecays;    p_harddecays    = NULL; }
  if (p_mihandler)     { delete p_mihandler;     p_mihandler     = NULL; }
  if (p_beamspectra)   { delete p_beamspectra;   p_beamspectra   = NULL; }
  if (p_model)         { delete p_model;         p_model         = NULL; }
  if (p_pythia)        { delete p_pythia;        p_pythia        = NULL; }
#ifdef USING__MCatNLO
  if (p_herwig)        { delete p_herwig;        p_herwig        = NULL; }
  if (p_mcatnlo)       { delete p_mcatnlo;       p_mcatnlo       = NULL; }
#endif
  if (p_analysis)      { delete p_analysis;      p_analysis      = NULL; }
  if (p_dataread)      { delete p_dataread;      p_dataread      = NULL; }
  std::set<Matrix_Element_Handler*> deleted;
  while (m_mehandlers.size()>0) {
    if (deleted.find(m_mehandlers.begin()->second)==deleted.end()) {
      deleted.insert(m_mehandlers.begin()->second);
      delete m_mehandlers.begin()->second;
    }
    m_mehandlers.erase(m_mehandlers.begin());
  }
  while (m_isrhandlers.size()>0) {
    delete m_isrhandlers.begin()->second;
    m_isrhandlers.erase(m_isrhandlers.begin());
  }
  PHASIC::Phase_Space_Handler::DeleteInfo();
}

bool Initialization_Handler::InitializeTheFramework(int nr)
{
  if (nr<=0) {
    ATOOLS::ParticleInit(m_path); 
  }

  bool okay = InitializeTheIO();
 
  if (m_mode==9999) {
    msg.Out()<<"SHERPA will read in the events."<<std::endl
	     <<"   The full framework is not needed."<<std::endl;
    InitializeTheAnalyses();
    return true;
  }
     
  okay = okay && InitializeTheModel();  

  //  set masses and widths from command line
  SetParameter(nr);
  UpdateParameters();
    
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
  }
  return okay;
}

bool Initialization_Handler::CheckBeamISRConsistency()
{
  if (p_model->Name()==std::string("ADD")) {
    double ms = p_model->ScalarConstant("M_s");
    if (ms<rpa.gen.Ecms()) {
      msg.Out()<<"WARNING in Initialization_Handler::CheckBeamISRConsistency :"<<std::endl
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
    msg.Error()<<"Error in Initialization of the Sherpa framework : "<<endl
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
  outfiles.push_back(p_dataread->GetValue<string>("SHERPA_OUTPUT",string("")));
  outfiles.push_back(p_dataread->GetValue<string>("HEPMC_OUTPUT",string("")));
  outfiles.push_back(p_dataread->GetValue<string>("HEPEVT_OUTPUT",string("")));
  std::string evtpath = p_dataread->GetValue<string>("EVT_FILE_PATH",m_path);
  int filesize        = p_dataread->GetValue<int>("FILE_SIZE",1000);
  for (int i=0;i<3;i++) {
    if (infiles[i]!=string("") || outfiles[i]!=string("")) {
      p_iohandler     = new IO_Handler(outfiles,infiles,evtpath,filesize);
      return true;
    }
  }
  std::string outmode = p_dataread->GetValue<string>("EVENT_MODE",string("Sherpa"));
  p_iohandler         = new IO_Handler(outmode);
  return true;
}

bool Initialization_Handler::InitializeTheExternalMC()
{
  std::string file;
  switch (m_mode) {
  case 9000: 
    p_pythia  = new Lund_Interface(m_path,m_evtfile,false);
    return true;
#ifdef USING__MCatNLO
  case 9001: 
    p_herwig  = new Herwig_Interface(m_path,m_evtfile,false);
    return true;
  case 9002: 
    p_mcatnlo = new MCatNLO_Interface(m_path,m_evtfile,false);
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
  if (p_model) { delete p_model; p_model = NULL; }
  Data_Read     * dataread     = new Data_Read(m_path+m_modeldat);
  // - add commandline parameter - !!
  Model_Handler * modelhandler = new Model_Handler();
  p_model                      = modelhandler->GetModel(dataread,m_path,m_modeldat);

  if (!p_model->RunSpectrumGenerator()) {
    msg.Error()<<"Error in Model_Initialization::Model_Initialization."<<endl
	       <<"    RunSpectrumGenerator() delivered false. Abort()."<<endl;
    abort();
  }

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
    if (m_isrhandlers.find(id)!=m_isrhandlers.end()) delete m_isrhandlers[id]; 
    Data_Read * dataread = new Data_Read(m_path+m_isrdat[i]);
    PDF_Handler * pdfhandler = new PDF_Handler();
    PDF_Base * pdfbase;
    ISR_Base ** isrbases = new ISR_Base*[2];
    double m_bunch_splimits[2];
    for (int j=0;j<2;++j) {
      pdfbase = pdfhandler->GetPDFLib(dataread,m_bunch_particles[j],j);
      if (pdfbase==NULL) isrbases[j] = new Intact(m_bunch_particles[j]);     
      else isrbases[j] = new Structure_Function(pdfbase,m_bunch_particles[j]);
      ATOOLS::rpa.gen.SetBunch(m_bunch_particles[j],j);
    }
    m_bunch_splimits[0] = dataread->GetValue<double>("ISR_SMIN",0.);
    m_bunch_splimits[1] = dataread->GetValue<double>("ISR_SMAX",1.);
    double kplimits[2];
    kplimits[0] = dataread->GetValue<double>("ISR_KPMIN",m_bunch_splimits[0]);
    kplimits[1] = dataread->GetValue<double>("ISR_KPMAX",m_bunch_splimits[1]);
    m_isrhandlers[id] = new ISR_Handler(isrbases,m_bunch_splimits,kplimits);
    m_isrhandlers[id]->SetBeamEnergy(p_beamspectra->GetBeam(0)->Energy(),0);
    m_isrhandlers[id]->SetBeamEnergy(p_beamspectra->GetBeam(1)->Energy(),1);
    msg_Info()<<"Initialized the ISR["<<id<<"] : "<<m_isrhandlers[id]->Type()<<endl;
    delete pdfhandler;
    delete dataread;
    if (!(p_beamspectra->CheckConsistency(m_bunch_particles))) {
      msg.Error()<<"Error in Environment::InitializeThePDFs()"<<endl
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
  MEHandlersMap::iterator it=m_mehandlers.find("SignalMEs");
  if (it!=m_mehandlers.end()) delete it->second;
  m_mehandlers["SignalMEs"]=me; 
  msg_Info()<<"Initialized the Matrix_Element_Handler for the hard processes :"<<me->Name()<<endl;
  if (p_analysis) {
    int weighted=1-me->EventGenerationMode();
    msg_Info()<<"Initialization_Handler::InitializeTheMatrixElements(): "
	      <<"Setting analysis mode "<<(weighted?"weighted":"unweighted")<<endl;
    p_analysis->SetWeighted(weighted);
  }
  return 1;
}

Matrix_Element_Handler * Initialization_Handler::GetMatrixElementHandler(std::string _key) { 
  MEHandlerIter pos = m_mehandlers.find(_key);
  if (pos!=m_mehandlers.end()) return pos->second;
  msg.Error()<<"Error in Initialization_Handler::GetMatrixElementHandler("<<_key<<") :"
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
  if (p_hadrondecays)  { delete p_hadrondecays;  p_hadrondecays  = NULL; }
  p_hadrondecays  = new Hadron_Decay_Handler(m_path,m_hadrondecaysdat,
  					     p_fragmentation->GetLundInterface());
  return 1;
}

bool Initialization_Handler::InitializeTheAnalyses()
{
  int helpi=p_dataread->GetValue<int>("SHOW_ANALYSIS_SYNTAX",0);
  if (helpi>0) {
    ATOOLS::msg.SetLevel(2);
    ATOOLS::Particle_Qualifier_Base::ShowQualifiers(helpi);
    ANALYSIS::Analysis_Handler::ShowSyntax(helpi);
    THROW(normal_exit,"Syntax shown.");
  }
  helpi=p_dataread->GetValue<int>("SHOW_QUALIFIER_SYNTAX",0);
  if (helpi>0) {
    ATOOLS::msg.SetLevel(2);
    ATOOLS::Particle_Qualifier_Base::ShowQualifiers(helpi);
    THROW(normal_exit,"Syntax shown.");
  }
  helpi=p_dataread->GetValue<int>("SHOW_VARIABLE_SYNTAX",0);
  if (helpi>0) {
    ATOOLS::msg.SetLevel(2);
    ATOOLS::Variable_Base<double>::ShowVariables(helpi);
    THROW(normal_exit,"Syntax shown.");
  }
  helpi=p_dataread->GetValue<int>("ANALYSIS",0);
  if (!helpi) return true;
  std::string outpath=p_dataread->GetValue<std::string>("ANALYSIS_OUTPUT","./");
  if (outpath==NotDefined<std::string>()) outpath="";
  p_analysis = new ANALYSIS::Analysis_Handler();
  p_analysis->SetInputPath(m_path);
  p_analysis->SetInputFile(m_analysisdat);
  p_analysis->SetOutputPath(outpath);
  return true;
}

bool Initialization_Handler::CalculateTheHardProcesses()
{
  if (m_mode>8999) {
    switch (m_mode) {
    case 9000:
      msg.Out()<<"SHERPA will generate the events through Pythia."<<std::endl
	       <<"   No cross sections for hard processes to be calculated."<<std::endl;
      return true;
#ifdef USING__MCatNLO
    case 9001:
      msg.Out()<<"SHERPA will generate the events through Herwig."<<std::endl
             <<"   No cross sections for hard processes to be calculated."<<std::endl;
      p_herwig->Initialize();
      return true;
    case 9002:
      msg.Out()<<"SHERPA will generate the events through MCatNLO."<<std::endl;
      p_mcatnlo->Initialize();
      return true;
#endif
    case 9999:
      msg.Out()<<"SHERPA will read in the events."<<std::endl
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
  msg.Out()<<"=========================================================================="<<std::endl
	   <<"Start calculating the hard cross sections. This may take some time.       "<<std::endl;
  int ok = me->CalculateTotalXSecs(scalechoice);
  if (ok && m_scan_istep!=-1) {
    AMEGIC::Process_Base * procs= me->GetAmegic()->Processes();
    msg.Out()<<ParameterValue()<<" ";
    for (size_t i=0; i<procs->Size();++i) {
      double xstot = (*procs)[i]->TotalXS()*rpa.Picobarn();
      msg.Out()<<xstot<<" ";
    }
    for (size_t i=0; i<procs->Size();++i) {
      msg.Out()<<"###"<<(*procs)[i]->Name();
    }
    msg.Out()<<endl;
  }
  if (ok) 
    msg.Out()<<"Calculating the hard cross sections has been successful.                  "<<std::endl
	     <<"=========================================================================="<<std::endl;
  else
    msg.Out()<<"Calculating the hard cross sections failed. Check this carefully.         "<<std::endl
	     <<"=========================================================================="<<std::endl;
  return ok;
}

void Initialization_Handler::SetParameter(int nr) {
  if (nr<0) return;

  if (nr!=m_scan_istep) 
    msg.Out()<<"WARNING: internal and external scan counter do not coincide "<<nr<<" vs. "<<m_scan_istep<<endl;

  bool logmode=false;
  if (m_scan_variable==string("YCUT")) {
    logmode=true;
  }

  if (logmode) {
    m_scan_value=m_scan_begin*exp(log(m_scan_end/m_scan_begin)*double(m_scan_istep)/double(m_scan_nsteps));
  }
  else {
    m_scan_value=m_scan_begin+(m_scan_end-m_scan_begin)*double(m_scan_istep)/double(m_scan_nsteps);
  }

  double value=m_scan_value;

  MyStrStream s;
  string sval;
  if (m_scan_variable==string("ECMS")) {
    s<<value/2.;
    s>>sval;
    msg.Out()<<" Setting Ecms/2 to : "<<sval<<endl;
    Data_Read::SetCommandLine("BEAM_ENERGY_1",sval);
    Data_Read::SetCommandLine("BEAM_ENERGY_2",sval);
    Read_Write_Base::AddCommandLine("BEAM_ENERGY_1 = "+sval+"; ");
    Read_Write_Base::AddCommandLine("BEAM_ENERGY_2 = "+sval+"; ");
  }
  else if (m_scan_variable.find("MASS(")!=string::npos || m_scan_variable.find("WIDTH(")!=string::npos ) {
    s<<value;
    s>>sval;
    m_options[m_scan_variable]=sval;
    // make sure UpdateParameters() is called
  }   
  else if (m_scan_variable==string("YCUT")) {
    rpa.gen.SetYcut(value);
  }
  else  {
    msg.Out()<<" Unknown Variable "<< m_scan_variable<<" in scan modus "<<endl;
    msg.Out()<<"  setting "<<m_scan_variable<<" = "<<value<<endl;
    s<<value; 
    s>>sval;
    m_options[m_scan_variable]=sval;    
    //    exit(1);
  }

  s.clear();
  double vmax=m_scan_end, vmin=m_scan_begin, vstep=(m_scan_end-m_scan_begin)/double(m_scan_nsteps);
  const  double ln10=log(10.);
  if (logmode) {
    value=log(value)/ln10;
    vmax=log(vmax)/ln10;
    vmin=log(vmin)/ln10;
    vstep=(vmax-vmin)/double(m_scan_nsteps);
  }
  double fac=exp(-ln10*int(log(vstep)/ln10-.999));
  value=dabs(value*fac);vmax=dabs(vmax*fac);vmin=dabs(vmin*fac);
  int nprec=int(log(Max(vmax,vmin))/ln10+.9999);
  
  s.width(nprec);
  s.fill('0');
  s<<int(value);
  s>>sval;

  std::string resdir="none";
  Parameter_Iterator it=m_options.find("#RESULT_DIRECTORY");
  if (it==m_options.end()) {
    Data_Read dataread(m_path+m_medat);
    resdir=dataread.GetValue<string>("RESULT_DIRECTORY",string("./Results"));
    m_options["#RESULT_DIRECTORY"]=resdir;
  }
  else {
    resdir=it->second;
  }
  size_t pos=resdir.find("$");
  if (pos!=std::string::npos) {
    if (pos!=resdir.size()) resdir=resdir.substr(0,pos)+sval+resdir.substr(pos+1);
    else resdir=resdir.substr(0,pos)+sval;
    std::cout<<"new RESULT_DIRECTORY = #"<<resdir<<"#\n";
    m_options["RESULT_DIRECTORY"]=resdir;
    ATOOLS::MakeDir(resdir,0755);
  }

  ++m_scan_istep;
}

int Initialization_Handler::ExtractCommandLineParameters(int argc,char * argv[])
{
  map<std::string,int> special_options;
  special_options["-V"]=12;
  special_options["--version"]=12;
  special_options["-?"]=13;
  special_options["-h"]=13;
  special_options["--help"]=13;
  special_options["-scan"]=14;
  special_options["-xsout"]=15;
  special_options["-eventout"]=16;

  special_options["PATH"]=101;
  special_options["RUNDATA"]=102;
  special_options["ECMS"]=103;
  special_options["PYTHIA"]=9000;
  special_options["EVTDATA"]=9999;
#ifdef USING__MCatNLO
  special_options["HERWIG"]=9001;
  special_options["MCatNLO"]=9002;
#endif
  

  

  for (int i=1; i<argc;++i) {
    int mode   = 0;
    string par = string(argv[i]);
    string key,value;
    int equal  = par.find("=");
    if (equal!=-1) {
      mode=1;
      value = par.substr(equal+1);
      key   = par = par.substr(0,equal);
      if (key.find("MASS")!=string::npos || key.find("WIDTH")!=string::npos) mode=100;
    }


    if (special_options.find(par)!=special_options.end()) 
      mode = special_options[par];

    // variables in dat files
    if (equal!=-1 && mode==1) {
      // perhaps check varible name first

      Data_Read::SetCommandLine(key,value);
      Read_Write_Base::AddCommandLine(key+" = "+value+"; ");
    }
    
    // special variables
    if (mode>=100) {
      MyStrStream s;
      switch (mode) {
       case 101:
	if (value[value.length()-1]!='/') value+=std::string("/");
	m_path=value;
	break;
      case 102:
	m_file=value;
	break;
      case 103:
	s<<value;
	double ecms;
	s>>ecms;
	s.clear();
	s<<ecms/2.;
	s>>value;
	msg.Out()<<" Setting ECMS/2 to : "<<value<<endl;
	Data_Read::SetCommandLine("BEAM_ENERGY_1",value);
	Data_Read::SetCommandLine("BEAM_ENERGY_2",value);
	Read_Write_Base::AddCommandLine("BEAM_ENERGY_1 = "+value+"; ");
	Read_Write_Base::AddCommandLine("BEAM_ENERGY_2 = "+value+"; ");
	break;
      case 9000:
	m_mode       = 9000;
	m_evtfile    = value;
	msg.Out()<<" Sherpa will produce Pythia events according to "<<value<<endl;
	break;
      case 9001:
	m_mode       = 9001;
	m_evtfile    = value;
	msg.Out()<<" Sherpa will produce Herwig events according to "<<value<<endl;
	break;
      case 9002:
	m_mode       = 9002;
	m_evtfile    = value;
	msg.Out()<<" Sherpa will produce MCatNLO events according to "<<value<<endl;
	break;
      case 9999:
	m_mode       = 9999;
	m_evtfile    = value;
	msg.Out()<<" Sherpa will read in events from : "<<value<<endl;
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
	  msg.Out()<<" Sherpa Version 1.0.5"<<endl;
	  msg.Out()<<"   employing: "<<endl;
	  msg.Out()<<"    * AMEGIC++ Version 2.0.5 "<<endl;
	  msg.Out()<<"    * APACIC++ Version 2.0.5 "<<endl;

	  string pyver("6.214");
	  msg.Out()<<"    * Pythia Version "<<pyver<<endl;
	  Char_Array40 cid=visaje_();
	  msg.Out()<<"    * IsaJet Version "<<cid.s<<endl;
	}
	exit(0);
      case 13:
	msg.Out()<<" Help: "<<endl;
	msg.Out()<<" Sherpa [options] [<variable>=<value>] "<<endl;
	msg.Out()<<endl;
	msg.Out()<<" Possible options: "<<endl;
	msg.Out()<<"  -V,--version   prints the Version number"<<endl;
	msg.Out()<<"  -?,--help      prints this help message"<<endl;
// 	msg.Out()<<"  -xsout <filename> "<<endl;
// 	msg.Out()<<"                 sets a file where calculated cross sections should be printed to"<<endl;
// 	msg.Out()<<"  -eventout <filename> "<<endl;
// 	msg.Out()<<"                 sets a file where events should be printed to"<<endl;	
	msg.Out()<<"  -scan <variable> <startvalue> <stopvalue> <number of steps>"<<endl;
	msg.Out()<<"                 performs a parameter scan"<<endl<<endl;
	msg.Out()<<"        <variable> ... in addition to all parameters in configuration files"<<endl
		 <<"                       PATH, RUNDATA, ECMS, MASS(<kfcode>) can be used"<<endl;
	exit(0);
      case 14: 
	// scan
	m_mode=14;
	Data_Read::SetCommandLine("EVENTS","0");
	Read_Write_Base::AddCommandLine("EVENTS = 0; ");
	if (i+4<argc) {
	  m_scan_variable=argv[++i];
	  MyStrStream s;
	  s<<argv[++i];
	  s>>m_scan_begin;
	  s.clear();
	  s<<argv[++i];
	  s>>m_scan_end;
	  s.clear();
	  s<<argv[++i];
	  s>>m_scan_nsteps;
	  m_scan_istep=0;
	  msg.Out()<<" scanning "<<m_scan_variable
	      <<" from "<<m_scan_begin<<" to "<<m_scan_end
	      <<" in "<<m_scan_nsteps<<" steps"<<endl;
	}
	else {
	  msg.Error()<<"ERROR:  missing scan parameter -scan"<<endl;
	  msg.Error()<<"       try Sherpa -? for more information "<<endl;
	  exit(1);
	}
	break;
      case 15:
	// xsout
	break;
      case 16:
	// eventout
	break;
      }
    }
  }
  return m_mode;
}


void Initialization_Handler::CheckFlagConsitency()
{
  Data_Read dr(m_path+m_medat);
  int  sudweight = dr.GetValue<int>("SUDAKOV_WEIGHT",0);

  // if SUDAKOV_WEIGHT=On
  if (sudweight>0) {
    //  Run.dat
    long nevt = p_dataread->GetValue<long>("EVENTS",0);
    if (nevt<=0) {
      Data_Read::SetCommandLine("EVENTS","1");
      Read_Write_Base::AddCommandLine("EVENTS = 1; ");
    }

    //  ME.dat 
    Data_Read::SetCommandLine("SCALE_SCHEME","65");
    Data_Read::SetCommandLine("KFACTOR_SCHEME","65");
    Data_Read::SetCommandLine("COUPLING_SCHEME","Running_alpha_S");
    Read_Write_Base::AddCommandLine("SCALE_SCHEME = 65; ");
    Read_Write_Base::AddCommandLine("KFACTOR_SCHEME = 65; ");
    Read_Write_Base::AddCommandLine("COUPLING_SCHEME = Running_alpha_S; ");

    //  Shower.dat
    Data_Read::SetCommandLine("FSR_SHOWER","1");
    Data_Read::SetCommandLine("FS_LOSEJETVETO","1");
    Read_Write_Base::AddCommandLine("FSR_SHOWER = 1; ");
    Read_Write_Base::AddCommandLine("FS_LOSEJETVETO = 1; ");
  }


  // check if all MI.dat / Decays.dat
  

}

int Initialization_Handler::UpdateParameters() 
{
  for (Parameter_Iterator it = m_options.begin(); it!=m_options.end() ; ++it) {
    MyStrStream s;
    string key=it->first;
    string value=it->second;
    if (key[0]=='#') continue;
    if (key.find("MASS(")!=string::npos || key.find("WIDTH(")!=string::npos ) {
      int a=key.find("(")+1;
      int b=key.find(")")-a;
      msg_Tracking()<<"Flavour "<<key.substr(a,b);
      s<<key.substr(a,b);
      int kfc;
      s>>kfc;
      s.clear();
      Flavour fl((kf::code)kfc);
      msg_Tracking()<<" : "<<fl<<endl;
      if (key.find("MASS(")!=string::npos) {
	double mass=fl.Mass();
	msg_Tracking()<<" old mass = "<<mass<<endl;
	s<<value;
	s>>mass;
	s.clear();
	msg_Tracking()<<" new mass = "<<mass<<endl;
	fl.SetMass(mass);
      }
      if (key.find(string("WIDTH("))!=string::npos) {
	msg_Tracking()<<"key:"<<key<<endl;
	double width=fl.Width();
	msg_Tracking()<<" old width = "<<width<<endl;
	s<<value;
	s>>width;
	s.clear();
	msg_Tracking()<<" new width = "<<width<<endl;
	fl.SetWidth(width);
      }
    }
    else {
      msg_Info()<<" update "<<key<<" = "<<value<<endl;
      Data_Read::SetCommandLine(key,value);
      Read_Write_Base::AddCommandLine(key+" = "+value+"; ");
    }
  }
  return 1;
}
