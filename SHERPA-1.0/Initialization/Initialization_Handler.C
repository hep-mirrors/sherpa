#include "Initialization_Handler.H"

#include "Model_Handler.H"
#include "Structure_Function.H"
#include "Intact.H"
#include "PDF_Handler.H"
#include "PDF_Base.H"

#include "Initial_State_Shower.H"

#include "Data_Read.H"
#include "Message.H"

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
  m_path(_path), m_file(_file),
  p_model(NULL), p_beamspectra(NULL), p_isrhandler(NULL),
  p_harddecays(NULL), p_showerhandler(NULL), p_beamremnants(NULL), 
  p_fragmentation(NULL), p_hadrondecays(NULL),  p_mihandler(NULL)
{
  m_scan_istep=-1;  

  p_dataread         = new Data_Read(m_path+m_file);
  m_modeldat         = p_dataread->GetValue<string>("MODEL_DATA_FILE",string("Model.dat"));
  m_beamdat          = p_dataread->GetValue<string>("BEAM_DATA_FILE",string("Beam.dat"));
  m_isrdat           = p_dataread->GetValue<string>("ISR_DATA_FILE",string("ISR.dat"));
  m_medat            = p_dataread->GetValue<string>("ME_DATA_FILE",string("ME.dat"));
  m_midat            = p_dataread->GetValue<string>("MI_DATA_FILE",string("MI.dat"));
  m_decaydat         = p_dataread->GetValue<string>("DECAY_DATA_FILE",string("Decays.dat"));
  m_showerdat        = p_dataread->GetValue<string>("SHOWER_DATA_FILE",string("Shower.dat"));
  m_beamremnantdat   = p_dataread->GetValue<string>("BEAMREMNANT_DATA_FILE",string("Beam.dat"));
  m_fragmentationdat = p_dataread->GetValue<string>("FRAGMENTATION_DATA_FILE",string("Fragmentation.dat"));
  m_hadrondecaysdat  = p_dataread->GetValue<string>("FRAGMENTATION_DATA_FILE",string("Fragmentation.dat"));
}


Initialization_Handler::Initialization_Handler(int argc,char * argv[]) : 
  p_model(NULL), p_beamspectra(NULL), p_isrhandler(NULL),
  p_harddecays(NULL), p_showerhandler(NULL), p_beamremnants(NULL), 
  p_fragmentation(NULL), p_hadrondecays(NULL), p_mihandler(NULL)
{
  m_path=std::string("./");
  m_file=std::string("Run.dat");

  m_scan_istep=-1;

  ExtractCommandLineParameters(argc, argv);

  p_dataread         = new Data_Read(m_path+m_file);
  m_modeldat         = p_dataread->GetValue<string>("MODEL_DATA_FILE",string("Model.dat"));
  m_beamdat          = p_dataread->GetValue<string>("BEAM_DATA_FILE",string("Beam.dat"));
  m_isrdat           = p_dataread->GetValue<string>("ISR_DATA_FILE",string("ISR.dat"));
  m_medat            = p_dataread->GetValue<string>("ME_DATA_FILE",string("ME.dat"));
  m_midat            = p_dataread->GetValue<string>("MI_DATA_FILE",string("MI.dat"));
  m_decaydat         = p_dataread->GetValue<string>("DECAY_DATA_FILE",string("Decays.dat"));
  m_showerdat        = p_dataread->GetValue<string>("SHOWER_DATA_FILE",string("Shower.dat"));
  m_beamremnantdat   = p_dataread->GetValue<string>("BEAMREMNANT_DATA_FILE",string("Beam.dat"));
  m_fragmentationdat = p_dataread->GetValue<string>("FRAGMENTATION_DATA_FILE",string("Fragmentation.dat"));
  m_hadrondecaysdat  = p_dataread->GetValue<string>("FRAGMENTATION_DATA_FILE",string("Fragmentation.dat"));
}


Initialization_Handler::~Initialization_Handler()
{
  if (p_hadrondecays)  { delete p_hadrondecays;  p_hadrondecays  = NULL; }
  if (p_fragmentation) { delete p_fragmentation; p_fragmentation = NULL; }
  if (p_beamremnants)  { delete p_beamremnants;  p_beamremnants  = NULL; }
  if (p_showerhandler) { delete p_showerhandler; p_showerhandler = NULL; }
  if (p_harddecays)    { delete p_harddecays;    p_harddecays    = NULL; }
  if (p_mihandler)     { delete p_mihandler;     p_mihandler     = NULL; }
  if (p_isrhandler)    { delete p_isrhandler;    p_isrhandler    = NULL; }
  if (p_beamspectra)   { delete p_beamspectra;   p_beamspectra   = NULL; }
  if (p_model)         { delete p_model;         p_model         = NULL; }
}


bool Initialization_Handler::InitializeTheFramework(int nr)
{
  if (nr<=0) {
    ATOOLS::ParticleInit(m_path); 
    rpa.Init(m_path,m_file);
  }
  bool okay = InitializeTheModel();  

  //  set masses and widths from command line
  SetParameter(nr);
  UpdateParameters();
    
  okay      = okay && InitializeTheBeams();
  okay      = okay && InitializeThePDFs();

  if (!CheckBeamISRConsistency()) return 0.;
  
  //  okay = okay && InitializeTheHardDecays();
  okay = okay && InitializeTheMatrixElements();
  //  only if events:
  if (rpa.gen.NumberOfEvents()>0) {
    okay = okay && InitializeTheShowers();
    okay = okay && InitializeTheFragmentation();
    okay = okay && InitializeTheUnderlyingEvents();
    okay = okay && InitializeTheHadronDecays();
    okay = okay && InitializeTheBeamRemnants();
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
  if (p_isrhandler->On()) {
    smin = Max(smin,p_isrhandler->SprimeMin());
    smax = Min(smax,p_isrhandler->SprimeMax());
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
    if (p_isrhandler->On()) {
      p_isrhandler->SetFixedSprimeMax(smax);
    } 
    else if (p_beamspectra->On()) {
      p_beamspectra->SetSprimeMax(smax);
    }
  }

  if (!(p_beamspectra->CheckConsistency(m_bunch_particles))) {
    msg.Error()<<"Error in Initialization of the Sherpa framework : "<<endl
	       <<"    Detected a mismatch of flavours from beams to bunches : "<<endl
	       <<"    "<<p_beamspectra->GetBeam(0)<<" -> "<<p_isrhandler->Flav(0)<<" and "
	       <<p_beamspectra->GetBeam(1)<<" -> "<<p_isrhandler->Flav(1)<<endl;
    return 0;
  }

  return 1;
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
  if (p_beamspectra) { 
    delete p_beamspectra; p_beamspectra = NULL; 

    // look for memory leaks
    /*
    for (int i=0; i<10;++i) {
      Data_Read * dataread = new Data_Read(m_path+m_beamdat);
      // - add commandline parameter - !!
      p_beamspectra        = new Beam_Spectra_Handler(dataread);
      delete p_beamspectra; p_beamspectra = NULL; 
      cout<<" Press Enter "<<endl;
      char key = cin.get();
    }
    */
  }
  Data_Read * dataread = new Data_Read(m_path+m_beamdat);
  // - add commandline parameter - !!
  p_beamspectra        = new Beam_Spectra_Handler(dataread);
  msg.Events()<<"Initialized the beams "<<p_beamspectra->Type()<<endl;

  delete dataread;  
  return 1;
}


bool Initialization_Handler::InitializeThePDFs()
{
  if (p_isrhandler)  { delete p_isrhandler; p_isrhandler = NULL; }
  Data_Read * dataread     = new Data_Read(m_path+m_isrdat);
  // - add commandline parameter - !!
  PDF_Handler * pdfhandler = new PDF_Handler();
  PDF_Base *  pdfbase;
  ISR_Base ** isrbases     = new ISR_Base*[2];
  double  m_bunch_splimits[2];

  for (int i=0;i<2;++i) {
    pdfbase = pdfhandler->GetPDFLib(dataread,m_bunch_particles[i],i);
    if (pdfbase==NULL) isrbases[i] = new Intact(m_bunch_particles[i]);     
                  else isrbases[i] = new Structure_Function(pdfbase,m_bunch_particles[i]);
  }
  m_bunch_splimits[0]      = dataread->GetValue<double>("ISR_SMIN",0.);
  m_bunch_splimits[1]      = dataread->GetValue<double>("ISR_SMAX",1.);
  p_isrhandler             = new ISR_Handler(isrbases,m_bunch_splimits);

  delete pdfhandler;
  delete dataread;

  if (!(p_beamspectra->CheckConsistency(m_bunch_particles))) {
    msg.Error()<<"Error in Environment::InitializeThePDFs()"<<endl
	       <<"   Inconsistent ISR & Beam:"<<endl
	       <<"   Abort program."<<endl;
    abort();
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
    me = new Matrix_Element_Handler(m_path,m_medat,p_model,p_beamspectra,p_isrhandler,
				    p_harddecays->GetMEHandler());
  }
  else {
    me = new Matrix_Element_Handler(m_path,m_medat,p_model,p_beamspectra,p_isrhandler,NULL);
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
  p_mihandler = new MI_Handler(m_path,m_midat,p_model,p_beamspectra,p_isrhandler);
  Matrix_Element_Handler *mehandler;
  if ((mehandler=p_mihandler->HardMEHandler())!=NULL) {
    m_mehandlers.insert(std::make_pair(std::string("MIMEs"),mehandler)); 
  }
  return true;
}

bool Initialization_Handler::InitializeTheShowers()
{
  if (p_showerhandler) { delete p_showerhandler; p_showerhandler = NULL; }
  int maxjets     = GetMatrixElementHandler(std::string("SignalMEs"))->MaxJets();
  p_showerhandler = new Shower_Handler(m_path,m_showerdat,p_model,p_isrhandler,maxjets);
  return 1;
}


bool Initialization_Handler::InitializeTheBeamRemnants() 
{
  if (p_beamremnants)  { delete p_beamremnants;  p_beamremnants  = NULL; }
  double scale=-4.0;
  if (p_showerhandler!=NULL) {
    if (p_showerhandler->GetApacic()->IniShower()) {
      scale=p_showerhandler->GetApacic()->IniShower()->CutOff();
    }
  }
  p_beamremnants = new Beam_Remnant_Handler(m_path,m_beamremnantdat,
					    p_isrhandler,p_beamspectra,scale);
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
					     p_fragmentation->GetLundFortranInterface());
  return 1;
}

bool Initialization_Handler::CalculateTheHardProcesses()
{
  int scalechoice = 0;
  if (p_showerhandler) {
    if (p_showerhandler->ISROn()) scalechoice += 1;
    if (p_showerhandler->FSROn()) scalechoice += 2;
  }
  Matrix_Element_Handler * me = GetMatrixElementHandler(std::string("SignalMEs"));
  int ok = me->CalculateTotalXSecs(scalechoice);
  if (ok && m_scan_istep!=-1) {
    AMEGIC::Process_Base * procs= me->GetAmegic()->Processes();
    cout<<ParameterValue()<<" ";
    for (int i=0; i<procs->Size();++i) {
      double xstot = (*procs)[i]->Total()*rpa.Picobarn();
      cout<<xstot<<" ";
    }
    for (int i=0; i<procs->Size();++i) {
      cout<<"###"<<(*procs)[i]->Name();
    }
    cout<<endl;
  }
  return ok;
}

bool Initialization_Handler::InitializeAllHardDecays()
{
  return 1;
  return p_harddecays->InitializeAllHardDecays(m_medat,p_model);
}

void Initialization_Handler::SetParameter(int nr) {
  if (nr<0) return;

  if (nr!=m_scan_istep) 
    cout<<"WARNING: internal and external scan counter do not coincide "<<nr<<" vs. "<<m_scan_istep<<endl;

  double value=m_scan_value=m_scan_begin+(m_scan_end-m_scan_begin)*double(m_scan_istep)/double(m_scan_nsteps);

  MyStrStream s;
  string sval;
  if (m_scan_variable==string("ECMS")) {
    s<<value/2.;
    s>>sval;
    cout<<" Setting Ecms/2 to : "<<sval<<endl;
    Data_Read::SetCommandLine("BEAM_ENERGY_1",sval);
    Data_Read::SetCommandLine("BEAM_ENERGY_2",sval);
  }
  else if (m_scan_variable.find("MASS")!=string::npos || m_scan_variable.find("WIDTH")!=string::npos ) {
    s<<value;
    s>>sval;
    m_options[m_scan_variable]=sval;
    // make sure UpdateParameters() is called
  }   
  else {
    cout<<" Unknown Variable "<< m_scan_variable<<" in scan modus "<<endl;
    cout<<"  setting "<<m_scan_variable<<" = "<<value<<endl;
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
    cout<<i<<" : "<<argv[i]<<" ->"<<mode<<endl;

    // variables in dat files
    if (equal!=-1 && mode==0) {
      cout<<equal<<":"<<key<<" = "<<value<<" ("<<par<<")"<<endl;
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
	cout<<" Setting Ecms to : "<<ecms<<endl;
	s<<ecms/2.;
	s>>value;
	cout<<" Setting Ecms/2 to : "<<value<<endl;
	Data_Read::SetCommandLine("BEAM_ENERGY_1",value);
	Data_Read::SetCommandLine("BEAM_ENERGY_2",value);
	break;
      case 100:
	m_options[key]=value;
      }
    }
    else {
      // other option
      switch (mode) {
      case 12:
	{
	  // should call a version roution
	  cout<<" Sherpa Version 1.0.2"<<endl;
	  cout<<"   employing: "<<endl;
	  cout<<"    * AMEGIC++ Version 2.0.2 "<<endl;
	  cout<<"    * APACIC++ Version 2.0.2 "<<endl;

	  string pyver("6.214");
	  cout<<"    * Pythia Version "<<pyver<<endl;
	  Char_Array40 cid=visaje_();
	  cout<<"    * IsaJet Version "<<cid.s<<endl;
	}
	exit(0);
      case 13:
	cout<<" Help: "<<endl;
	cout<<" Sherpa [options] [<variable>=<value>] "<<endl;
	cout<<endl;
	cout<<" Possible options: "<<endl;
	cout<<"  -V,--version   prints the Version number"<<endl;
	cout<<"  -?,--help      prints this help message"<<endl;
	cout<<"  -xsout <filename> "<<endl;
	cout<<"                 sets a file where calculated cross sections should be printed to"<<endl;
	cout<<"  -eventout <filename> "<<endl;
	cout<<"                 sets a file where events should be printed to"<<endl;	
	cout<<"  -scan <variable> <startvalue> <stopvalue> <number of steps>"<<endl;
	cout<<"                 performs a parameter scan"<<endl;
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
	  cout<<" scanning "<<m_scan_variable
	      <<" from "<<m_scan_begin<<" to "<<m_scan_end
	      <<" in "<<m_scan_nsteps<<" steps"<<endl;
	}
	else {
	  cout<<"ERROR: missing scan parameter -scan"<<endl;
	  cout<<"       try Sherpa -? for more information "<<endl;
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
    cout<<" "<<key<<" = "<<value<<endl;
    int a=key.find("(")+1;
    int b=key.find(")")-a;
    cout<<"Flavour "<<key.substr(a,b);
    s<<key.substr(a,b);
    int kfc;
    s>>kfc;
    Flavour fl((kf::code)kfc);
    cout<<" : "<<fl<<endl;
    if (key.find("MASS")!=string::npos) {
      double mass=fl.Mass();
      cout<<" old mass = "<<mass<<endl;
      s<<value;
      s>>mass;
      cout<<" new mass = "<<mass<<endl;
      fl.SetMass(mass);
    }
    if (key.find(string("WIDTH"))!=string::npos) {
      cout<<"key:"<<key<<endl;
      double width=fl.Width();
      cout<<" old width = "<<width<<endl;
      s<<value;
      s>>width;
      cout<<" new width = "<<width<<endl;
      fl.SetWidth(width);
    }
  }
  return 1;
}
