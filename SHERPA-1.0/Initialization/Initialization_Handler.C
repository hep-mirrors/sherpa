#include "Initialization_Handler.H"

#include "Model_Handler.H"
#include "Structure_Function.H"
#include "Intact.H"
#include "PDF_Handler.H"
#include "PDF_Base.H"

#include "Data_Read.H"
#include "Message.H"

using namespace SHERPA;
using namespace MODEL;
using namespace BEAM;
using namespace PDF;
using namespace AORGTOOLS;
using namespace APHYTOOLS;
using namespace AMATOOLS;
using namespace std;



Initialization_Handler::Initialization_Handler(string _path,string _file) : 
  m_path(_path), m_file(_file),
  p_model(NULL), p_beamspectra(NULL), p_isrhandler(NULL), p_mehandler(NULL),
  p_showerhandler(NULL), p_beamremnants(NULL), p_fragmentation(NULL),
  p_hadrondecays(NULL)
{
  p_dataread         = new Data_Read(m_path+m_file);
  m_modeldat         = p_dataread->GetValue<string>("MODEL_DATA_FILE",string("Model.dat"));
  m_beamdat          = p_dataread->GetValue<string>("BEAM_DATA_FILE",string("Beam.dat"));
  m_isrdat           = p_dataread->GetValue<string>("ISR_DATA_FILE",string("ISR.dat"));
  m_medat            = p_dataread->GetValue<string>("ME_DATA_FILE",string("ME.dat"));
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
  if (p_mehandler)     { delete p_mehandler;     p_mehandler     = NULL; }
  if (p_isrhandler)    { delete p_isrhandler;    p_isrhandler    = NULL; }
  if (p_beamspectra)   { delete p_beamspectra;   p_beamspectra   = NULL; }
  if (p_model)         { delete p_model;         p_model         = NULL; }
}


bool Initialization_Handler::InitializeTheFramework()
{
  APHYTOOLS::ParticleInit(m_path); 
  rpa.Init(m_path);

  bool okay = InitializeTheModel();  
  okay      = okay && InitializeTheBeams();
  okay      = okay && InitializeThePDFs();

  if (!(p_beamspectra->CheckConsistency(m_bunch_particles))) {
    msg.Error()<<"Error in Initialization of the Sherpa framework : "<<endl
	       <<"    Detected a mismatch of flavours from beams to bunches : "<<endl
	       <<"    "<<p_beamspectra->GetBeam(0)<<" -> "<<p_isrhandler->Flav(0)<<" and "
	       <<p_beamspectra->GetBeam(1)<<" -> "<<p_isrhandler->Flav(1)<<endl;
    return 0;
  }
  
  okay = okay && InitializeTheMatrixElements();
  okay = okay && InitializeTheShowers();
  okay = okay && InitializeTheBeamRemnants();
  okay = okay && InitializeTheFragmentation();
  okay = okay && InitializeTheHadronDecays();

  return okay;
}



bool Initialization_Handler::InitializeTheModel()
{
  Data_Read     * dataread     = new Data_Read(m_path+m_modeldat);
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
  Data_Read * dataread = new Data_Read(m_path+m_beamdat);
  p_beamspectra        = new Beam_Spectra_Handler(dataread);
  msg.Events()<<"Initialized the beams "<<p_beamspectra->Type()<<endl;

  delete dataread;  
  return 1;
}


bool Initialization_Handler::InitializeThePDFs()
{
  Data_Read * dataread     = new Data_Read(m_path+m_isrdat);
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


bool Initialization_Handler::InitializeTheMatrixElements()
{
  p_mehandler = new Matrix_Element_Handler(m_path,m_medat,
					   p_model,p_beamspectra,p_isrhandler);
  return 1;
}

bool Initialization_Handler::InitializeTheShowers()
{
  p_showerhandler = new Shower_Handler(m_path,m_showerdat,p_model,
				       p_isrhandler,p_mehandler->MaxJets());
  return 1;
}


bool Initialization_Handler::InitializeTheBeamRemnants() 
{
  p_beamremnants = new Beam_Remnant_Handler(m_path,m_beamremnantdat,
					    p_isrhandler,p_beamspectra);
  return 1;
}

bool Initialization_Handler::InitializeTheFragmentation() 
{
  p_fragmentation = new Fragmentation_Handler(m_path,m_fragmentationdat);
  return 1;
}

bool Initialization_Handler::InitializeTheHadronDecays() 
{
  p_hadrondecays  = new Hadron_Decay_Handler(m_path,m_hadrondecaysdat,
					     p_fragmentation->GetLundFortranInterface());
  return 1;
}

bool Initialization_Handler::CalculateTheHardProcesses()
{
  return p_mehandler->CalculateTotalXSecs();
}
