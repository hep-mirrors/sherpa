#include "Initialization_Handler.H"

#include "Model_Handler.H"
#include "Structure_Function.H"
#include "Intact.H"
#include "PDF_Handler.H"
#include "PDF_Base.H"

#include "Initial_State_Shower.H"
#include "LL_Branching.H"

#include "Data_Read.H"
#include "Message.H"
#include "Scaling.H"

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
  p_model(NULL), p_beamspectra(NULL), 
  p_harddecays(NULL), p_showerhandler(NULL), p_beamremnants(NULL), 
  p_fragmentation(NULL), p_hadrondecays(NULL), p_mihandler(NULL),
  p_pythia(NULL)
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
  p_pythia(NULL)
{
  m_path=std::string("./");
  m_file=std::string("Run.dat");

  m_scan_istep=-1;

  ExtractCommandLineParameters(argc, argv);

  if (m_mode>8999) {
    p_dataread         = new Data_Read(m_path+m_file);
    m_analysisdat      = p_dataread->GetValue<string>("ANALYSIS_DATA_FILE",string("Analysis.dat"));
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
}


Initialization_Handler::~Initialization_Handler()
{
  while (m_analyses.size()>0) {
    delete m_analyses.begin()->second;
    m_analyses.erase(m_analyses.begin());
  }
  if (p_outputhandler) { delete p_outputhandler; p_outputhandler = NULL; }
  if (p_hadrondecays)  { delete p_hadrondecays;  p_hadrondecays  = NULL; }
  if (p_fragmentation) { delete p_fragmentation; p_fragmentation = NULL; }
  if (p_beamremnants)  { delete p_beamremnants;  p_beamremnants  = NULL; }
  if (p_showerhandler) { delete p_showerhandler; p_showerhandler = NULL; }
  if (p_harddecays)    { delete p_harddecays;    p_harddecays    = NULL; }
  if (p_mihandler)     { delete p_mihandler;     p_mihandler     = NULL; }
  if (p_beamspectra)   { delete p_beamspectra;   p_beamspectra   = NULL; }
  if (p_model)         { delete p_model;         p_model         = NULL; }
  if (p_pythia)        { delete p_pythia;        p_pythia        = NULL; }
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
  PDF::LL_Branching::DeleteSplittings();
}


bool Initialization_Handler::InitializeTheFramework(int nr)
{
  if (nr<=0) {
    rpa.Init(m_path,m_file);
    ATOOLS::ParticleInit(m_path); 
  }

  bool okay = InitializeTheIO();
  std::cout<<"m_mode = "<<m_mode<<std::endl;
  if (m_mode>8999) {
    okay &= InitializeTheExternalMC();
    InitializeTheAnalyses();
    return true;
  }

  okay      = okay && InitializeTheModel();  

  //  set masses and widths from command line
  SetParameter(nr);
  UpdateParameters();
    
  okay = okay && InitializeTheBeams();
  okay = okay && InitializeThePDFs();

  ATOOLS::Integration_Info *info=PHASIC::Phase_Space_Handler::GetInfo();
  m_isrhandlers[isr::hard_process]->AssignKeys(info);
  if (m_isrhandlers.find(isr::hard_subprocess)!=m_isrhandlers.end()) {
    m_isrhandlers[isr::hard_subprocess]->AssignKeys(info);
  }

  if (!CheckBeamISRConsistency()) return 0.;
  
  okay = okay && InitializeTheHardDecays();
  okay = okay && InitializeTheMatrixElements();
  //  only if events:
  if (rpa.gen.NumberOfEvents()>0) {
    okay = okay && InitializeTheShowers();
    okay = okay && InitializeTheFragmentation();
    okay = okay && InitializeTheUnderlyingEvents();
    okay = okay && InitializeTheHadronDecays();
    okay = okay && InitializeTheBeamRemnants();
    okay = okay && InitializeTheAnalyses();
  }
  return okay;
}

bool Initialization_Handler::CheckBeamISRConsistency()
{
  if (p_model->Name()==std::string("ADD")) {
    double ms = p_model->ScalarConstant("M_s");
    if (ms<rpa.gen.Ecms()) {
      msg.Out()<<" WARNING : You are using the ADD model beyond its valid range ! "<<endl;
      msg.Out()<<" WARNING : You are using the ADD model beyond its valid range ! "<<endl;
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
  int filesize;
  std::vector<std::string> infiles, outfiles;
  infiles.push_back(p_dataread->GetValue<string>("SHERPA_INPUT",string("")));
  infiles.push_back(p_dataread->GetValue<string>("HEPMC_INPUT",string("")));
  infiles.push_back(p_dataread->GetValue<string>("HEPEVT_INPUT",string("")));
  outfiles.push_back(p_dataread->GetValue<string>("SHERPA_OUTPUT",string("")));
  outfiles.push_back(p_dataread->GetValue<string>("HEPMC_OUTPUT",string("")));
  outfiles.push_back(p_dataread->GetValue<string>("HEPEVT_OUTPUT",string("")));
  filesize  = p_dataread->GetValue<int>("FILE_SIZE",1000);
  p_outputhandler = new Output_Handler(outfiles,infiles);
  return true;
}

bool Initialization_Handler::InitializeTheExternalMC()
{
  std::string file;
  switch (m_mode) {
  case 9000: 
    p_pythia = new Pythia_Interface(m_path,m_evtfile);
    return true;
  default: 
    m_mode = 9999;
    msg.Out()<<"SHERPA will read in the events."<<std::endl
	     <<"   The full framework is not needed."<<std::endl;
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
  msg.Events()<<"Initialized the model : "<<p_model->Name()<<endl;

  delete modelhandler;
  delete dataread;  
  return 1;
}


bool Initialization_Handler::InitializeTheBeams()
{
  if (p_beamspectra) { delete p_beamspectra; p_beamspectra = NULL; }
  Data_Read * dataread = new Data_Read(m_path+m_beamdat);
  p_beamspectra        = new Beam_Spectra_Handler(dataread);
  msg.Events()<<"Initialized the beams "<<p_beamspectra->Type()<<endl;
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
    delete pdfhandler;
    delete dataread;
    if (!(p_beamspectra->CheckConsistency(m_bunch_particles))) {
      msg.Error()<<"Error in Environment::InitializeThePDFs()"<<endl
		 <<"   Inconsistent ISR & Beam:"<<endl
		 <<"   Abort program."<<endl;
      abort();
    }
  }
  return 1;
}

bool Initialization_Handler::InitializeTheHardDecays()
{
  if (p_harddecays)    { delete p_harddecays;    p_harddecays    = NULL; }
  p_harddecays = new Hard_Decay_Handler(m_path,m_decaydat,m_medat,p_model);
  if (p_harddecays->GetMEHandler()!=NULL) {
    ATOOLS::msg.Tracking()<<"Initialized the Hard_Decay_Handler. Its ME_Handler is : "
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
  m_mehandlers.insert(std::make_pair(std::string("SignalMEs"),me)); 
  ATOOLS::msg.Tracking()<<"Initialized the Hard_Decay_Handler. Its ME_Handler is : "
			<<me->Name()<<"/"<<me<<endl;
  return 1;
}

Matrix_Element_Handler * Initialization_Handler::GetMatrixElementHandler(std::string _key) { 
  MEHandlerIter pos = m_mehandlers.find(_key);
  if (pos!=m_mehandlers.end()) return pos->second;
  ATOOLS::msg.Error()<<"Error in Initialization_Handler::GetMatrixElementHandler("<<_key<<") :"
		     <<"   Key not found. Return Null pointer."<<endl;
  return NULL;
}


bool Initialization_Handler::InitializeTheUnderlyingEvents()
{
  p_mihandler = new MI_Handler(m_path,m_midat,p_model,p_beamspectra,
			       m_isrhandlers[isr::hard_subprocess]);
  Matrix_Element_Handler *mehandler;
  if ((mehandler=p_mihandler->HardMEHandler())!=NULL) {
    m_mehandlers.insert(std::make_pair(std::string("MIMEs"),mehandler)); 
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
  if (p_showerhandler) { delete p_showerhandler; p_showerhandler = NULL; }
  int maxjets     = GetMatrixElementHandler(std::string("SignalMEs"))->MaxJets();
  p_showerhandler = new Shower_Handler(m_path,m_showerdat,p_model,
				       m_isrhandlers[isr::hard_process],maxjets);
  return 1;
}


bool Initialization_Handler::InitializeTheBeamRemnants() 
{
  if (p_beamremnants)  { delete p_beamremnants;  p_beamremnants  = NULL; }
  double scale=-4.0;
  if (p_showerhandler!=NULL && p_showerhandler->GetApacic()) {
    if (p_showerhandler->GetApacic()->IniShower()) {
      scale=p_showerhandler->GetApacic()->IniShower()->CutOff();
    }
  }
  p_beamremnants = new Beam_Remnant_Handler(m_path,m_beamremnantdat,
					    m_isrhandlers[isr::hard_process],
					    p_beamspectra,scale);
  return 1;
}

bool Initialization_Handler::InitializeTheFragmentation() 
{
  if (p_fragmentation) { delete p_fragmentation; p_fragmentation = NULL; }
  p_fragmentation = new Fragmentation_Handler(m_path,m_fragmentationdat);
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
  if (rpa.gen.Analysis()<=0) {
    msg.Error()<<"Warning in Initialization_Handler::InitializeTheAnalyses()."<<std::endl
	       <<"   Analysis is switched off - continue run."<<std::endl;
    return 1;
  } 
  Data_Read dr(m_path+"/Run.dat");
  std::string prefix=dr.GetValue<std::string>("ANALYSIS_OUTPUT");
  if (prefix==NotDefined<std::string>()) prefix="";


  ifstream * from = new ifstream((m_path+m_analysisdat).c_str());
  if (!from->good()) {
    msg.Error()<<"Error in Initialization_Handler::InitializeTheAnalyses()."<<std::endl
	       <<"   File : "<<(m_path+m_analysisdat).c_str()<<" not found ! Abort program execution."<<endl;
    abort();
  }
  std::string buffer, phase;
  Sample_Analysis * sa = NULL;
  unsigned int pos;
  bool add;
  for (;;) {
    if (from->eof()) break;
    getline(*from,buffer);
    buffer += std::string(" ");
    if (buffer[0] != '%' && buffer.length()>0) {
      if (buffer.find("ANALYSIS_PHASE =")!=std::string::npos) {
	pos    = buffer.find("ANALYSIS_PHASE =");
	buffer = buffer.substr(pos+16);
	while(buffer.length()>0) {
	  if (buffer[0]==' ') buffer = buffer.substr(1);
	  else {
	    pos = buffer.find(string(" "));
	    if (pos>0) phase = buffer.substr(0,pos);
	    break;
	  }
	}
	add = false;
	if (phase!=std::string("")) {
	  sa  = new Sample_Analysis(from,phase,prefix);
	  if (sa->On()) add = true;
	  else delete sa;
	}
	if (add) { 
	  int i=0;
	  std::string name;
	  do {
	    name=phase+std::string("_")+ATOOLS::ToString(++i);
	  } while (m_analyses.find(name)!=m_analyses.end());
	  m_analyses.insert(std::make_pair(name,sa)); 
	  ATOOLS::msg.Tracking()<<"Initialized Sample_Analysis "<<name<<std::endl;
	}
      }
    }
  }
  return 1;
}

bool Initialization_Handler::CalculateTheHardProcesses()
{
  if (m_mode>8999) {
    switch (m_mode) {
    case 9000:
      msg.Out()<<"SHERPA will generate the events through Pyrthia."<<std::endl
	       <<"   No cross sections for hard processes to be calculated."<<std::endl;
      return true;
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
  int ok = me->CalculateTotalXSecs(scalechoice);
  if (ok && m_scan_istep!=-1) {
    AMEGIC::Process_Base * procs= me->GetAmegic()->Processes();
    ATOOLS::msg.Out()<<ParameterValue()<<" ";
    for (size_t i=0; i<procs->Size();++i) {
      double xstot = (*procs)[i]->TotalXS()*rpa.Picobarn();
      ATOOLS::msg.Out()<<xstot<<" ";
    }
    for (size_t i=0; i<procs->Size();++i) {
      ATOOLS::msg.Out()<<"###"<<(*procs)[i]->Name();
    }
    ATOOLS::msg.Out()<<endl;
  }
  return ok;
}

void Initialization_Handler::SetParameter(int nr) {
  if (nr<0) return;

  if (nr!=m_scan_istep) 
    ATOOLS::msg.Out()<<"WARNING: internal and external scan counter do not coincide "<<nr<<" vs. "<<m_scan_istep<<endl;

  double value=m_scan_value=m_scan_begin+(m_scan_end-m_scan_begin)*double(m_scan_istep)/double(m_scan_nsteps);

  MyStrStream s;
  string sval;
  if (m_scan_variable==string("ECMS")) {
    s<<value/2.;
    s>>sval;
    ATOOLS::msg.Out()<<" Setting Ecms/2 to : "<<sval<<endl;
    Data_Read::SetCommandLine("BEAM_ENERGY_1",sval);
    Data_Read::SetCommandLine("BEAM_ENERGY_2",sval);
  }
  else if (m_scan_variable.find("MASS(")!=string::npos || m_scan_variable.find("WIDTH(")!=string::npos ) {
    s<<value;
    s>>sval;
    m_options[m_scan_variable]=sval;
    // make sure UpdateParameters() is called
  }   
  else {
    ATOOLS::msg.Out()<<" Unknown Variable "<< m_scan_variable<<" in scan modus "<<endl;
    ATOOLS::msg.Out()<<"  setting "<<m_scan_variable<<" = "<<value<<endl;
    s<<value; 
    s>>sval;
    m_options[m_scan_variable]=sval;    
    //    exit(1);
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
    ATOOLS::msg.Out()<<i<<" : "<<argv[i]<<" ->"<<mode<<endl;

    // variables in dat files
    if (equal!=-1 && mode==1) {
      ATOOLS::msg.Out()<<equal<<":"<<key<<" = "<<value<<" ("<<par<<")"<<endl;
      Data_Read::SetCommandLine(key,value);
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
	ATOOLS::msg.Out()<<" Setting Ecms to : "<<ecms<<endl;
	s<<ecms/2.;
	s>>value;
	ATOOLS::msg.Out()<<" Setting Ecms/2 to : "<<value<<endl;
	Data_Read::SetCommandLine("BEAM_ENERGY_1",value);
	Data_Read::SetCommandLine("BEAM_ENERGY_2",value);
	break;
      case 9000:
	m_mode       = 9000;
	m_evtfile    = value;
	ATOOLS::msg.Out()<<" Sherpa will produce Pythia events according to "<<value<<endl;
	break;
      case 9999:
	m_mode       = 9999;
	m_evtfile    = value;
	ATOOLS::msg.Out()<<" Sherpa will read in events from : "<<value<<endl;
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
	  ATOOLS::msg.Out()<<" Sherpa Version 1.0.3"<<endl;
	  ATOOLS::msg.Out()<<"   employing: "<<endl;
	  ATOOLS::msg.Out()<<"    * AMEGIC++ Version 2.0.3 "<<endl;
	  ATOOLS::msg.Out()<<"    * APACIC++ Version 2.0.3 "<<endl;

	  string pyver("6.214");
	  ATOOLS::msg.Out()<<"    * Pythia Version "<<pyver<<endl;
	  Char_Array40 cid=visaje_();
	  ATOOLS::msg.Out()<<"    * IsaJet Version "<<cid.s<<endl;
	}
	exit(0);
      case 13:
	ATOOLS::msg.Out()<<" Help: "<<endl;
	ATOOLS::msg.Out()<<" Sherpa [options] [<variable>=<value>] "<<endl;
	ATOOLS::msg.Out()<<endl;
	ATOOLS::msg.Out()<<" Possible options: "<<endl;
	ATOOLS::msg.Out()<<"  -V,--version   prints the Version number"<<endl;
	ATOOLS::msg.Out()<<"  -?,--help      prints this help message"<<endl;
	ATOOLS::msg.Out()<<"  -xsout <filename> "<<endl;
	ATOOLS::msg.Out()<<"                 sets a file where calculated cross sections should be printed to"<<endl;
	ATOOLS::msg.Out()<<"  -eventout <filename> "<<endl;
	ATOOLS::msg.Out()<<"                 sets a file where events should be printed to"<<endl;	
	ATOOLS::msg.Out()<<"  -scan <variable> <startvalue> <stopvalue> <number of steps>"<<endl;
	ATOOLS::msg.Out()<<"                 performs a parameter scan"<<endl;
	exit(0);
      case 14: 
	// scan
	m_mode=14;
	Data_Read::SetCommandLine("EVENTS","0");
	if (i+4<argc) {
	  m_scan_variable=argv[++i];
	  MyStrStream s;
	  s<<argv[++i];
	  s>>m_scan_begin;
	  s<<argv[++i];
	  s>>m_scan_end;
	  s<<argv[++i];
	  s>>m_scan_nsteps;
	  m_scan_istep=0;
	  ATOOLS::msg.Out()<<" scanning "<<m_scan_variable
	      <<" from "<<m_scan_begin<<" to "<<m_scan_end
	      <<" in "<<m_scan_nsteps<<" steps"<<endl;
	}
	else {
	  ATOOLS::msg.Out()<<"ERROR: missing scan parameter -scan"<<endl;
	  ATOOLS::msg.Out()<<"       try Sherpa -? for more information "<<endl;
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

int Initialization_Handler::UpdateParameters() 
{
  for (Parameter_Iterator it = m_options.begin(); it!=m_options.end() ; ++it) {
    MyStrStream s;
    string key=it->first;
    string value=it->second;
    ATOOLS::msg.Out()<<" "<<key<<" = "<<value<<endl;
    int a=key.find("(")+1;
    int b=key.find(")")-a;
    ATOOLS::msg.Out()<<"Flavour "<<key.substr(a,b);
    s<<key.substr(a,b);
    int kfc;
    s>>kfc;
    Flavour fl((kf::code)kfc);
    ATOOLS::msg.Out()<<" : "<<fl<<endl;
    if (key.find("MASS")!=string::npos) {
      double mass=fl.Mass();
      ATOOLS::msg.Out()<<" old mass = "<<mass<<endl;
      s<<value;
      s>>mass;
      ATOOLS::msg.Out()<<" new mass = "<<mass<<endl;
      fl.SetMass(mass);
    }
    if (key.find(string("WIDTH"))!=string::npos) {
      ATOOLS::msg.Out()<<"key:"<<key<<endl;
      double width=fl.Width();
      ATOOLS::msg.Out()<<" old width = "<<width<<endl;
      s<<value;
      s>>width;
      ATOOLS::msg.Out()<<" new width = "<<width<<endl;
      fl.SetWidth(width);
    }
  }
  return 1;
}
