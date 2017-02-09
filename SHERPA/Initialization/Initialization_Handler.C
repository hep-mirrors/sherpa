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
#include "PHASIC++/Scales/Core_Scale_Setter.H"
#include "MODEL/Main/Model_Base.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "METOOLS/Currents/C_Spinor.H"
#include "PDF/Main/Structure_Function.H"
#include "PDF/Main/Intact.H"
#include "PDF/Main/PDF_Base.H"
#include "AMISIC++/Main/MI_Base.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Default_Reader.H"
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
#include "ATOOLS/Org/Command_Line_Interface.H"

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
  m_mode(eventtype::StandardPerturbative), 
  m_savestatus(false), p_model(NULL), p_beamspectra(NULL), 
  p_mehandler(NULL), p_harddecays(NULL), p_beamremnants(NULL),
  p_fragmentation(NULL), p_softcollisions(NULL), p_hdhandler(NULL), 
  p_mihandler(NULL), p_softphotons(NULL), p_evtreader(NULL),
  p_variations(NULL)
{
  m_path=std::string("");
  m_file=std::string("Run.dat");

  ExtractCommandLineParameters(argc, argv);

  SetFileNames();

  if (p_dataread->Read<std::string>(m_evtform,"EVENT_INPUT", "")) {
    m_mode=eventtype::EventReader;
    msg_Out()<<" Sherpa will read in events as "<<m_evtform<<endl;
  }

  ATOOLS::s_loader->SetCheck(p_dataread->Get<int>("CHECK_LIBLOCK",0));

  rpa->Init(m_path,m_file,argc,argv);
  CheckVersion();
  LoadLibraries();
  ShowParameterSyntax();
  ran->InitExternal(m_path,m_file);

  rpa->gen.SetSoftSC(p_dataread->Get<int>("SOFT_SPIN_CORRELATIONS",0));
  std::string hdstr(p_dataread->GetStringNormalisingNoneLikeValues("HARD_DECAYS","None"));
  int defhsc = !(hdstr=="None");
  rpa->gen.SetHardSC(p_dataread->Get<int>("HARD_SPIN_CORRELATIONS",defhsc));
  exh->AddTerminatorObject(this);
}

void Initialization_Handler::SetFileNames()
{
  p_dataread    = new Default_Reader;
  p_dataread->SetInputPath(m_path);
  p_dataread->SetInputFile(m_file);
  std::string fname(m_file);
  if (fname.find("|")!=std::string::npos) 
    fname=fname.substr(0,fname.find("|"));
  m_modeldat         = p_dataread->Get<string>("MODEL_DATA_FILE",fname+"|(model){|}(model)");
  m_beamdat          = p_dataread->Get<string>("BEAM_DATA_FILE",fname+"|(beam){|}(beam)");
  m_isrdat[0]        = p_dataread->Get<string>("ISR_DATA_FILE",fname+"|(isr){|}(isr)");
  m_isrdat[1]        = p_dataread->Get<string>("MI_ISR_DATA_FILE",m_isrdat[0]);
  m_medat            = p_dataread->Get<string>("ME_DATA_FILE",fname+"|(me){|}(me)");
  m_midat            = p_dataread->Get<string>("MI_DATA_FILE",fname+"|(mi){|}(mi)");
  m_showerdat        = p_dataread->Get<string>("SHOWER_DATA_FILE",fname+"|(shower){|}(shower)");
  m_beamremnantdat   = p_dataread->Get<string>("BEAMREMNANT_DATA_FILE",fname+"|(beam){|}(beam)");
  m_fragmentationdat = p_dataread->Get<string>("FRAGMENTATION_DATA_FILE",fname+"|(fragmentation){|}(fragmentation)");
  m_softcollisiondat = p_dataread->Get<string>("SOFTCOLLISIONS_DATA_FILE",string("SoftCollisions.dat"));
  m_hadrondecaysdat  = p_dataread->Get<string>("FRAGMENTATION_DATA_FILE",fname+"|(fragmentation){|}(fragmentation)");
  m_softphotonsdat   = p_dataread->Get<string>("SOFT_PHOTON_DATA_FILE",fname+"|(fragmentation){|}(fragmentation)");
  m_processesdat     = p_dataread->Get<string>("PROCESSFILE",fname+"|(processes){|}(processes)");
  m_selectordat      = p_dataread->Get<string>("SELECTORFILE",fname+"|(selector){|}(selector)");
  m_analysisdat      = p_dataread->Get<string>("ANALYSIS_DATA_FILE",FileExists("Analysis.dat")?
						    "Analysis.dat":fname+"|(analysis){|}(analysis)");
  std::string integrationdat = p_dataread->Get<string>("INTEGRATION_DATA_FILE",fname+"|(integration){|}(integration)");
  std::string momentadat     = p_dataread->Get<string>("MOMENTA_DATA_FILE",fname+"|(momenta){|}(momenta)");
  if (FileExists("Momenta.dat")) momentadat="Momenta.dat";

  rpa->gen.SetVariable("MODEL_DATA_FILE",m_modeldat);
  rpa->gen.SetVariable("ME_DATA_FILE",m_medat);
  rpa->gen.SetVariable("MODEL_DATA_FILE",m_modeldat);
  rpa->gen.SetVariable("SHOWER_DATA_FILE",m_showerdat);
  rpa->gen.SetVariable("INTEGRATION_DATA_FILE",integrationdat);
  rpa->gen.SetVariable("FRAGMENTATION_DATA_FILE",m_fragmentationdat);
  rpa->gen.SetVariable("MOMENTA_DATA_FILE",momentadat);
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
  if (p_beamspectra)   { delete p_beamspectra;   p_beamspectra   = NULL; }
  if (p_model)         { delete p_model;         p_model         = NULL; }
  if (p_dataread)      { delete p_dataread;      p_dataread      = NULL; }
  if (p_variations)    { delete p_variations;    p_variations    = NULL; }
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
    if (exit==NULL) THROW(fatal_error,"Cannot unload PDF library "+*pdflib);
    ((PDF_Exit_Function)exit)();
  }
}

void Initialization_Handler::CheckVersion()
{
  std::vector<std::string> versioninfo;
  p_dataread->ReadVector(versioninfo,"SHERPA_VERSION");
  if (!versioninfo.size()) return;
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

void Initialization_Handler::LoadLibraries() const
{
  Default_Reader reader;
  reader.SetInputPath(m_path);
  reader.SetInputFile(m_file);
  std::vector<std::string> ldadd;
  if (!reader.ReadVector(ldadd,"SHERPA_LDADD")) return;
  for (size_t i(0);i<ldadd.size();++i) {
    if (!s_loader->LoadLibrary(ldadd[i])) {
      THROW(fatal_error,"Cannot load extra library.");
    }
    else msg_Info()<<METHOD<<"(): Library lib"<<ldadd[i]<<".so loaded.\n";
  }
}

void Initialization_Handler::ShowParameterSyntax()
{
  Default_Reader reader;
  int helpi(reader.Get("SHOW_ME_GENERATORS", 0));
  if (helpi>0) {
    msg->SetLevel(2);
    PHASIC::ME_Generator_Base::ShowSyntax(helpi);
    THROW(normal_exit,"Syntax shown.");
  }
  helpi = reader.Get("SHOW_PS_GENERATORS", 0);
  if (helpi>0) {
    msg->SetLevel(2);
    PHASIC::Channel_Generator::ShowSyntax(helpi);
    THROW(normal_exit,"Syntax shown.");
  }
  helpi = reader.Get("SHOW_NLOMC_GENERATORS", 0);
  if (helpi>0) {
    msg->SetLevel(2);
    PDF::NLOMC_Base::ShowSyntax(helpi);
    THROW(normal_exit,"Syntax shown.");
  }
  helpi = reader.Get("SHOW_SHOWER_GENERATORS", 0);
  if (helpi>0) {
    msg->SetLevel(2);
    PDF::Shower_Base::ShowSyntax(helpi);
    THROW(normal_exit,"Syntax shown.");
  }
  helpi = reader.Get("SHOW_SCALE_SYNTAX", 0);
  if (helpi>0) {
    msg->SetLevel(2);
    if (helpi&1) PHASIC::Scale_Setter_Base::ShowSyntax(helpi);
    if (helpi&2) PHASIC::Core_Scale_Setter::ShowSyntax(helpi);
    THROW(normal_exit,"Syntax shown.");
  }
  helpi = reader.Get("SHOW_SELECTOR_SYNTAX", 0);
  if (helpi>0) {
    msg->SetLevel(2);
    PHASIC::Selector_Base::ShowSyntax(helpi);
    THROW(normal_exit,"Syntax shown.");
  }
  helpi = reader.Get("SHOW_MODEL_SYNTAX", 0);
  if (helpi>0) {
    msg->SetLevel(2);
    MODEL::Model_Base::ShowSyntax(helpi);
    THROW(normal_exit,"Syntax shown.");
  }
  helpi = reader.Get("SHOW_ANALYSIS_SYNTAX", 0);
  if (helpi>0) {
    msg->SetLevel(2);
    InitializeTheAnalyses();
    for (Analysis_Vector::iterator it=m_analyses.begin(); it!=m_analyses.end(); ++it)
      (*it)->ShowSyntax(helpi);
    THROW(normal_exit,"Syntax shown.");
  }
  helpi = reader.Get("SHOW_VARIABLE_SYNTAX", 0);
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
  
  if (m_mode==eventtype::StandardPerturbative) {
  std::string eventtype;
  eventtype = p_dataread->Get<std::string>("EVENT_TYPE", "StandardPerturbative");
  if (eventtype=="StandardPerturbative") 
    m_mode=eventtype::StandardPerturbative;
  else if (eventtype=="MinimumBias") {
    m_mode=eventtype::MinimumBias;
    Read_Write_Base::AddCommandLine("MI_HANDLER None;");
  }
  else if (eventtype=="HadronDecay") {
    m_mode=eventtype::HadronDecay;
    Read_Write_Base::AddCommandLine("MI_HANDLER None;");
  }
  else {
    THROW(not_implemented,"Unknown event type '"+eventtype+"'");
  }
  }
  okay = okay && InitializeTheBeams();
  okay = okay && InitializeThePDFs();
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
      (m_evtform,Input_Arguments(m_path,infile,p_dataread,
				 p_model,m_isrhandlers[isr::hard_process]));
    if (p_evtreader==NULL) THROW(fatal_error,"Event reader not found");
    msg_Events()<<"SHERPA will read in the events."<<std::endl
  		<<"   The full framework is not needed."<<std::endl;
    InitializeTheBeamRemnants();
    InitializeTheIO();
    return true;
  }
  PHASIC::Phase_Space_Handler::GetInfo();
  if (rpa->gen.NumberOfEvents()>0) {
  okay = okay && InitializeTheFragmentation();
  okay = okay && InitializeTheSoftCollisions();
  }
  okay = okay && InitializeTheShowers();
  okay = okay && InitializeTheMatrixElements();
  okay = okay && InitializeTheBeamRemnants();
  okay = okay && InitializeTheHardDecays();
  //  only if events:
  if (rpa->gen.NumberOfEvents()>0) {
    okay = okay && InitializeTheHadronDecays();
    okay = okay && InitializeTheUnderlyingEvents();
    okay = okay && InitializeTheSoftPhotons();
    okay = okay && InitializeTheIO();
    okay = okay && InitializeTheReweighting();
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
  std::string outpath=p_dataread->Get<std::string>("EVT_FILE_PATH",".");
  std::string format=p_dataread->GetStringNormalisingNoneLikeValues("EVENT_OUTPUT","None");
  std::vector<std::string> outputs;
  Data_Reader readline(",",";","#","");
  std::string stag(rpa->gen.Variable("RNG_SEED"));
  while (stag.find(' ')!=std::string::npos) stag.replace(stag.find(' '),1,"-");
  readline.AddTag("RNG_SEED",stag);
  readline.SetString(format);
  readline.StringVectorFromStringNormalisingNoneLikeValues(outputs);
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
      (outputs[i],Output_Arguments(outpath,outfile,p_dataread));
    if (out==NULL) {
      if (!s_loader->LoadLibrary("Sherpa"+libname+"Output")) 
	THROW(missing_module,"Cannot load output library Sherpa"+libname+"Output.");
      out=Output_Base::Getter_Function::GetObject
	(outputs[i],Output_Arguments(outpath,outfile,p_dataread));
    }
    if (out==NULL) THROW(fatal_error,"Cannot initialize "+outputs[i]+" output");
    m_outputs.push_back(out);
  }
  return true;
}

bool Initialization_Handler::InitializeTheModel()
{
  if (p_model) delete p_model;
  Default_Reader reader;
  reader.SetInputPath(m_path);
  reader.SetInputFile(m_modeldat);
  std::string name(reader.Get<std::string>("MODEL", "SM"));
  p_model=Model_Base::Model_Getter_Function::
    GetObject(name,Model_Arguments(m_path,m_modeldat,true));
  if (p_model==NULL) {
    if (!s_loader->LoadLibrary("Sherpa"+name))
      THROW(missing_module,"Cannot load model library Sherpa"+name+".");
    p_model=Model_Base::Model_Getter_Function::
      GetObject(name,Model_Arguments(m_path,m_modeldat,true));
  }
  if (p_model==NULL) THROW(not_implemented,"Model not implemented");
  MODEL::s_model=p_model;
  return 1;
}


bool Initialization_Handler::InitializeTheBeams()
{
  if (p_beamspectra) { delete p_beamspectra; p_beamspectra = NULL; }
  Default_Reader reader;
  reader.SetInputPath(m_path);
  reader.SetInputFile(m_beamdat);
  p_beamspectra        = new Beam_Spectra_Handler(&reader);
  p_beamspectra->Output();
  
  return 1;
}

bool Initialization_Handler::InitializeThePDFs()
{
  Default_Reader reader;
  reader.SetInputPath(m_path);
  reader.SetInputFile(m_isrdat[0]);

  // load PDF libraries
  std::string defset[2];
  for (int beam(0);beam<=1;++beam) {
    std::string deflib("None");
    if (p_beamspectra->GetBeam(beam)->Bunch().Kfcode()==kf_p_plus) {
      deflib="NNPDFSherpa";
      defset[beam]="NNPDF30NNLO";
    }
    else if (p_beamspectra->GetBeam(beam)->Bunch().Kfcode()==kf_e) {
      deflib="PDFESherpa";
      defset[beam]="PDFe";
    }
    else if (p_beamspectra->GetBeam(beam)->Bunch().IsPhoton()) {
      deflib="GRVSherpa";
      defset[beam]="GRV";
    }
    std::vector<std::string> pdflibs;
    std::string mpilib, beamlib;
    reader.ReadStringVectorNormalisingNoneLikeValues(pdflibs,"PDF_LIBRARY");
    if (pdflibs.size()==0) m_pdflibs.insert(deflib);
    for (size_t i(0);i<pdflibs.size();++i) m_pdflibs.insert(pdflibs[i]);
    if (reader.Read<std::string>(mpilib,"MI_PDF_LIBRARY","")) {
      m_pdflibs.insert(mpilib);
    }
    if (reader.Read<std::string>(beamlib,"PDF_LIBRARY_"+ToString(beam+1),"")) {
      m_pdflibs.insert(beamlib);
    }
  }
  if (Variations::NeedsLHAPDF6Interface(m_path)) {
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
  int helpi(reader.Get("SHOW_PDF_SETS", 0));
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
    reader.SetInputFile(m_isrdat[i]);
    PDF_Base * pdfbase;
    ISR_Base ** isrbases = new ISR_Base*[2];
    double m_bunch_splimits[2];
    for (int j=0;j<2;++j) {
      std::string indextag("_" + ToString(j + 1));
      int defaultflav(p_beamspectra->GetBeam(j)->Bunch());
      int flav = reader.Get<int>("BUNCH_"+ToString(j+1),defaultflav);
      m_bunch_particles[j] = Flavour((kf_code)abs(flav));
      if (flav<0) m_bunch_particles[j] = m_bunch_particles[j].Bar();

      std::string set;
      int version(0);
      bool specializedformi(false);

      // Read PDF set and version
      set = reader.Get<std::string>("PDF_SET", defset[j]);
      set = reader.Get<std::string>("PDF_SET" + indextag, set);
      version = reader.Get<int>("PDF_SET_MEMBER", 0);
      version = reader.Get<int>("PDF_SET_VERSION", version);
      version = reader.Get<int>("PDF_SET_MEMBER" + indextag, version);
      version = reader.Get<int>("PDF_SET_VERSION" + indextag, version);

      // Read PDF set and version for hard subprocesses (multiple interactions)
      if (id==isr::hard_subprocess) {
        std::string mpiset;
        if (reader.Read<std::string>(mpiset, "MI_PDF_SET" + indextag, set)
            || reader.Read<std::string>(mpiset, "MI_PDF_SET", set)) {
          set = mpiset;
          // If using a special MI PDF set, then do not use normal PDF set
          // version
          version = 0;
          specializedformi = true;
        }
        int mpiversion;
        if (reader.Read<int>(mpiversion, "MI_PDF_SET_VERSION" + indextag, version)
            || reader.Read<int>(mpiversion, "MI_PDF_SET_MEMBER" + indextag, version)
            || reader.Read<int>(mpiversion, "MI_PDF_SET_VERSION", version)
            || reader.Read<int>(mpiversion, "MI_PDF_SET_MEMBER", version)) {
          version = mpiversion;
          specializedformi = true;
        }
      }

      // Load PDF
      pdfbase = PDF_Base::PDF_Getter_Function::GetObject
        (set,PDF_Arguments(m_bunch_particles[j],&reader, j, set, version));
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
    m_bunch_splimits[0] = reader.Get<double>("ISR_SMIN",1e-10);
    m_bunch_splimits[1] = reader.Get<double>("ISR_SMAX",1.);
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

bool Initialization_Handler::InitializeTheHardDecays()
{
  Default_Reader reader;
  reader.SetInputPath(m_path);
  reader.SetInputFile(m_medat);
  std::string decays=reader.GetStringNormalisingNoneLikeValues("HARD_DECAYS","None");
  if (decays=="None") return true;

  if (p_harddecays)    { delete p_harddecays;    p_harddecays    = NULL; }
  p_harddecays = new Hard_Decay_Handler(m_path,m_medat);
  return 1;
}

bool Initialization_Handler::InitializeTheMatrixElements()
{
  if (p_mehandler) delete p_mehandler;
  p_mehandler = new Matrix_Element_Handler(m_path,m_medat,m_processesdat,m_selectordat);
  p_mehandler->SetShowerHandler(m_showerhandlers[isr::hard_process]);
  int ret(p_mehandler->InitializeProcesses(p_model,p_beamspectra,m_isrhandlers[isr::hard_process]));
  msg_Info()<<"Initialized the Matrix_Element_Handler for the hard processes."
            <<endl;
  return ret==1;
}

bool Initialization_Handler::InitializeTheUnderlyingEvents()
{
  as->SetActiveAs(isr::hard_subprocess);
  p_mihandler = new MI_Handler(m_path,m_midat,p_model,p_beamspectra,
			       m_isrhandlers[isr::hard_subprocess]);
  p_mihandler->SetShowerHandler(m_showerhandlers[isr::hard_process]);
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
      (m_path, m_showerdat, p_model, m_isrhandlers[isrtypes[i]],i);
  }
  as->SetActiveAs(isr::hard_process);
  msg_Info()<<"Initialized the Shower_Handler."<<endl;
  return 1;
}


bool Initialization_Handler::InitializeTheSoftCollisions() 
{
  if (p_softcollisions) { delete p_softcollisions; p_softcollisions = NULL; }
  p_softcollisions = new Soft_Collision_Handler(m_path,m_softcollisiondat,
						p_beamspectra,
						m_isrhandlers[isr::hard_process]);
  msg_Info()<<"Initialized the Soft_Collision_Handler."<<endl;
  return 1;
}

bool Initialization_Handler::InitializeTheBeamRemnants() 
{
  if (p_beamremnants)  delete p_beamremnants;
  p_beamremnants = 
    new Beam_Remnant_Handler(m_path,m_beamremnantdat,
			     p_beamspectra,
			     m_isrhandlers[isr::hard_process],
			     p_softcollisions);
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
  Default_Reader reader;
  reader.SetInputPath(m_path);
  reader.SetInputFile(m_hadrondecaysdat);
  std::string frag=reader.GetStringNormalisingNoneLikeValues("FRAGMENTATION", string("Ahadic"));
  if (frag=="None") return true;

  std::string defdecmodel("Hadrons");
  if (frag=="Lund") defdecmodel="Lund";
  string decmodel = reader.GetStringNormalisingNoneLikeValues("DECAYMODEL",defdecmodel);
  msg_Tracking()<<"Decaymodel = "<<decmodel<<std::endl;
  if (decmodel=="None") return true;
  else if (decmodel==std::string("Hadrons")) {
    as->SetActiveAs(isr::hard_subprocess);
    Hadron_Decay_Handler* hd=new Hadron_Decay_Handler(m_path,m_hadrondecaysdat);
    as->SetActiveAs(isr::hard_process);
    p_hdhandler=hd;
  }
  else if ((decmodel==string("Lund")) ) {
#ifdef USING__PYTHIA
    as->SetActiveAs(isr::hard_subprocess);
    Lund_Interface * lund(NULL);
    if (p_fragmentation->GetLundInterface()==NULL) {
      string lfile = reader.Get<string>("LUND_FILE","Lund.dat");
      lund = new Lund_Interface(m_path,lfile);
    }
    else lund = p_fragmentation->GetLundInterface();
    Lund_Decay_Handler* hd=new Lund_Decay_Handler(lund,m_path,m_hadrondecaysdat);
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
  p_softphotons = new Soft_Photon_Handler(m_path,m_softphotonsdat,p_mehandler);
  if (p_harddecays) p_harddecays->SetSoftPhotonHandler(p_softphotons);
  if (p_hdhandler)  p_hdhandler->SetSoftPhotonHandler(p_softphotons);
  msg_Info()<<"Initialized the Soft_Photon_Handler."<<endl;
  return true;
}

bool Initialization_Handler::InitializeTheAnalyses()
{
  std::string outpath=p_dataread->Get<std::string>("ANALYSIS_OUTPUT","Analysis/");
  std::vector<std::string> analyses;
  p_dataread->ReadStringVectorNormalisingNoneLikeValues(analyses, "ANALYSIS");
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
    Analysis_Interface* ana=Analysis_Interface::Analysis_Getter_Function::GetObject
                            (analyses[i],Analysis_Arguments(m_path,m_analysisdat,outpath));
    if (ana==NULL) {
      if (!s_loader->LoadLibrary("Sherpa"+analyses[i]+"Analysis")) 
	THROW(missing_module,"Cannot load Analysis library '"+analyses[i]+"'.");
      ana=Analysis_Interface::Analysis_Getter_Function::GetObject
	(analyses[i],Analysis_Arguments(m_path,m_analysisdat,outpath));
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
  Default_Reader dataread;
  dataread.SetInputPath(m_path);
  Variations::CheckConsistencyWithBeamSpectra(p_beamspectra);
  p_variations = new Variations(&dataread);
  msg_Info()<<"Initialized the Reweighting."<<endl;
  return true;
}

bool Initialization_Handler::CalculateTheHardProcesses()
{
  if (m_mode!=eventtype::StandardPerturbative) return true;
  
  msg_Events()<<"===================================================================\n"
              <<"Start calculating the hard cross sections. This may take some time.\n";
  Default_Reader reader;
  ATOOLS::msg->SetLevel(reader.Get<int>("INT_OUTPUT",ATOOLS::msg->Level()));
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
  Default_Reader reader;
  reader.SetInputPath(m_path);
  reader.SetInputFile(m_medat);
  double sf(reader.Get<double>("SCALE_FACTOR",1.));
  rpa->gen.SetVariable("FACTORIZATION_SCALE_FACTOR",
		      ToString(sf*reader.Get<double>("FACTORIZATION_SCALE_FACTOR",1.0)));
  rpa->gen.SetVariable("RENORMALIZATION_SCALE_FACTOR",
		      ToString(sf*reader.Get<double>("RENORMALIZATION_SCALE_FACTOR",1.0)));
  msg_Debugging()<<METHOD<<"(): Set scale factors {\n"
		 <<"  fac scale: "<<rpa->gen.Variable("FACTORIZATION_SCALE_FACTOR")<<"\n"
		 <<"  ren scale: "<<rpa->gen.Variable("RENORMALIZATION_SCALE_FACTOR")<<"\n}\n";
  int cmode=reader.Get<int>("METS_CLUSTER_MODE",0);
  rpa->gen.SetVariable("METS_CLUSTER_MODE",ToString(cmode));
  if (cmode!=0) msg_Info()<<METHOD<<"(): Set cluster mode "<<cmode<<".\n";
  Default_Reader css_reader;
  css_reader.SetInputPath(m_path);
  css_reader.SetInputFile(m_showerdat);
  int evol          = css_reader.Get<int>("CSS_EVOLUTION_SCHEME",1);
  int kfmode        = css_reader.Get<int>("CSS_KFACTOR_SCHEME",1);
  int scs           = css_reader.Get<int>("CSS_SCALE_SCHEME",0);
  int svmode        = css_reader.Get<double>("CSS_SCALE_VARIATION_SCHEME",1);
  double k0sqf      = css_reader.Get<double>("CSS_FS_PT2MIN",1.0);
  double k0sqi      = css_reader.Get<double>("CSS_IS_PT2MIN",2.00);
  double fs_as_fac  = css_reader.Get<double>("CSS_FS_AS_FAC",1.0);
  double is_as_fac  = css_reader.Get<double>("CSS_IS_AS_FAC",0.5);
  double as_var_fac = css_reader.Get<double>("CSS_SCALE_FACTOR",1.);
  double mth        = css_reader.Get<double>("CSS_MASS_THRESHOLD",0.0);
  rpa->gen.SetVariable("CSS_EVOLUTION_SCHEME",ToString(evol));
  rpa->gen.SetVariable("CSS_KFACTOR_SCHEME",ToString(kfmode));
  rpa->gen.SetVariable("CSS_SCALE_SCHEME",ToString(scs));
  rpa->gen.SetVariable("CSS_SCALE_VARIATION_SCHEME",ToString(svmode));
  rpa->gen.SetVariable("CSS_FS_PT2MIN",ToString(k0sqf));
  rpa->gen.SetVariable("CSS_IS_PT2MIN",ToString(k0sqi));
  rpa->gen.SetVariable("CSS_FS_AS_FAC",ToString(fs_as_fac));
  rpa->gen.SetVariable("CSS_IS_AS_FAC",ToString(is_as_fac));
  rpa->gen.SetVariable("CSS_SCALE_FACTOR",ToString(as_var_fac));
  rpa->gen.SetVariable("CSS_MASS_THRESHOLD",ToString(mth));
}

void Initialization_Handler::ExtractCommandLineParameters(int argc,char * argv[])
{
  // NOTE: The implementation here is essentially glue code between the old
  // input handling and the new Command_Line_Interface class. The old input
  // handling will essentially be superseded and this implementation will
  // simplify greatly.
  Command_Line_Interface cli(argc, argv);


  // handle parameters with an immediate effect on the Initialization_Handler
  // and delete them after use so they are not re-processed below
  std::string datpath;
  if (cli.GetParameterValue("STATUS_PATH") != "") {
    datpath = cli.GetParameterValue("STATUS_PATH");
    if (datpath[datpath.size() - 1] != '/') datpath += "/";
    cli.SetParameterValue("STATUS_PATH", "");
  }

  if (cli.GetParameterValue("PATH") != "") {
    m_path = cli.GetParameterValue("PATH");
    if (m_path[m_path.size() - 1] != '/') m_path += "/";
    cli.SetParameterValue("PATH", "");
  }

  if (cli.GetParameterValue("RUNDATA") != "") {
    m_file = cli.GetParameterValue("RUNDATA");
    cli.SetParameterValue("RUNDATA", "");
  }

  if (cli.GetParameterValue("SAVE_STATUS") != "") {
    std::string savestatus(cli.GetParameterValue("SAVE_STATUS"));
    if (savestatus[savestatus.size() - 1] != '/') savestatus += "/";
    rpa->gen.SetVariable("SHERPA_STATUS_PATH",
                         rpa->gen.Variable("SHERPA_RUN_PATH") + "/" + savestatus);
    m_savestatus=true;
    cli.SetParameterValue("SAVE_STATUS", "");
  }

  // set default section, i.e. run section, of Run card
  if (m_file.find("|")==std::string::npos) {
    Read_Write_Base cf(1,0," ",";","!","=");
    cf.SetInputPath(m_path);
    cf.SetInputFile(m_file+"|(run){|}(run)");
    if (cf.OpenInFile()) m_file+="|(run){|}(run)";
  }

  // store status path
  if (datpath!="") m_path = datpath;
  rpa->gen.SetVariable("PATH_PIECE", m_path);
  m_path = "";

  // contract unparsed parameters and tags into a vector of strings, such that
  // they can be fed below into Read_Write_Base
  std::vector<std::string> helpsv;
  typedef std::map<std::string, std::string>::const_iterator It;
  for (It it(cli.GetParameters().begin());
       it != cli.GetParameters().end();
       ++it) {
    helpsv.push_back(it->first + "=" + it->second);
  }
  for (It it(cli.GetTags().begin());
       it != cli.GetTags().end();
       ++it) {
    helpsv.push_back(it->first + ":=" + it->second);
  }

  // add parameters from the default section in Run.dat
  // (this makes it possible to overwrite particle properties in Run.dat)
  // helpsv2 will then holds all Run.dat/(run) settings and the unparsed
  // command line arguments, in that order
  std::vector<std::string> helpsv2;
  Data_Reader reader;
  reader.SetInputPath(m_path);
  reader.SetInputFile(m_file);
  std::vector<std::vector<std::string> > helpsvv;
  if (reader.MatrixFromFile(helpsvv,"")) {
    size_t oldsize(helpsv2.size());
    helpsv2.resize(oldsize+helpsvv.size());
    for (size_t i(0);i<helpsvv.size();++i) {
      helpsv2[oldsize+i]=helpsvv[i][0];
      for (size_t j(1);j<helpsvv[i].size();++j)
        helpsv2[oldsize+i]+=" "+helpsvv[i][j];
    }
  }
  helpsv2.insert(helpsv2.end(),helpsv.begin(),helpsv.end());

  // add all elements in helpsv2 either to the Read_Write_Base command line
  // or to its global tags (if the ':=' delimiter is found)
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
      if (key=="TUNE") {
        THROW(not_implemented,"Currently TUNE is not supported.");
        //SetTuneParameters(value);
      }
    }
    else {
      Read_Write_Base::AddCommandLine(par+";");
    }
  }
  rpa->gen.SetVariable("RUN_DATA_FILE",m_file);
}

/// Disabled for release 2.2.0
void Initialization_Handler::SetTuneParameters(const std::string tune)
{
  std::vector<std::string> tuneparams;
  if (tune == "NNPDF23" ||
      tune == "NNPDF23_UEup" || tune == "NNPDF23_UEdown") {
    THROW(fatal_error,"Currently there is no such tune.");
    tuneparams.push_back("PDF_LIBRARY                  = LHAPDFSherpa");
    tuneparams.push_back("PDF_SET                      = NNPDF23_nlo_as_0119.LHgrid");
    tuneparams.push_back("K_PERP_MEAN_1                = 1.08");
    tuneparams.push_back("K_PERP_MEAN_2                = 1.08");
    tuneparams.push_back("K_PERP_SIGMA_1               = 1.10");
    tuneparams.push_back("K_PERP_SIGMA_2               = 1.10");
    tuneparams.push_back("PROFILE_PARAMETERS           = 0.44 0.93");
    tuneparams.push_back("RESCALE_EXPONENT             = 0.208");
    tuneparams.push_back("SCALE_MIN                    = 2.63");
    if (tune == "NNPDF23_UEup") {
      tuneparams.push_back("SIGMA_ND_FACTOR              = 0.358");
      Read_Write_Base::AddCommandLine("MI_RESULT_DIRECTORY_SUFFIX _up;");
    }
    else if (tune == "NNPDF23_UEdown") {
      tuneparams.push_back("SIGMA_ND_FACTOR              = 0.418");
      Read_Write_Base::AddCommandLine("MI_RESULT_DIRECTORY_SUFFIX _down;");
    }
    else {
      tuneparams.push_back("SIGMA_ND_FACTOR              = 0.388");
    }
    tuneparams.push_back("CSS_IS_AS_FAC                = 0.872");
    tuneparams.push_back("CSS_IS_PT2MIN                = 2.21");
    tuneparams.push_back("COLOUR_RECONNECTION_STRENGTH = 0.25");
  }
  else if (tune == "CT10" ||
           tune == "CT10_UEup" || tune == "CT10_UEdown") {
    tuneparams.push_back("PDF_LIBRARY                  = CT10Sherpa");
    tuneparams.push_back("PDF_SET                      = ct10");
    tuneparams.push_back("K_PERP_MEAN_1                = 1.10");
    tuneparams.push_back("K_PERP_MEAN_2                = 1.10");
    tuneparams.push_back("K_PERP_SIGMA_1               = 0.85");
    tuneparams.push_back("K_PERP_SIGMA_2               = 0.85");
    tuneparams.push_back("PROFILE_PARAMETERS           = 0.76 0.58");
    tuneparams.push_back("RESCALE_EXPONENT             = 0.244");
    tuneparams.push_back("SCALE_MIN                    = 2.44");
    if (tune == "CT10_UEup") {
      tuneparams.push_back("SIGMA_ND_FACTOR              = 0.4104909");
      Read_Write_Base::AddCommandLine("MI_RESULT_DIRECTORY_SUFFIX _up;");
    }
    else if (tune == "CT10_UEdown") {
      tuneparams.push_back("SIGMA_ND_FACTOR              = 0.4882269");
      Read_Write_Base::AddCommandLine("MI_RESULT_DIRECTORY_SUFFIX _down;");
    }
    else {
      tuneparams.push_back("SIGMA_ND_FACTOR              = 0.4452459");
    }
    tuneparams.push_back("CSS_IS_AS_FAC                = 0.50");
    tuneparams.push_back("CSS_IS_PT2MIN                = 4.78");
    tuneparams.push_back("COLOUR_RECONNECTION_STRENGTH = 0.23");
  }
  else {
    msg_Error()<<"Ignoring unknown tune name \"" << tune << "\"" << std::endl;
    return;
  }
  msg_Out()<<"******************************************************" << std::endl;
  msg_Out()<<"****" << std::endl;
  msg_Out()<<"**** Setting tune parameters for " << tune << std::endl;
  msg_Out()<<"****" << std::endl;
  for (size_t i=0; i<tuneparams.size(); i++) {
    msg_Out()<<"**** " << tuneparams[i] << std::endl;
    Read_Write_Base::AddCommandLine(tuneparams[i]);
  }
  msg_Out()<<"****" << std::endl;
  msg_Out()<<"**** Note that these parameters might get overwritten on the command line" << std::endl;
  msg_Out()<<"**** or by parameters set appearing later in the run card." << std::endl;
  msg_Out()<<"****" << std::endl;
  msg_Out()<<"******************************************************" << std::endl;
}
