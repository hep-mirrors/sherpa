#include "Environment.H"
#include "Message.H"
#include "Run_Parameter.H"

#include "Model_Handler.H"
#include "PDF_Handler.H"
#include "PDF_Base.H"
#include "Structure_Function.H"
#include "Intact.H"

using namespace AMEGIC;
using namespace ATOOLS;
using namespace MODEL;
using namespace BEAM;
using namespace PDF;


Environment::Environment(std::string _path,std::string _file) : 
  m_path(_path), m_file(_file)
{
  msg.Tracking()<<"Initialize Initialization_Handler for "<<m_path+m_file<<std::endl;
  p_dataread = new Data_Read(m_path+m_file);
  m_modeldat = p_dataread->GetValue("MODEL_DATA_FILE",std::string("Model.dat"));
  m_beamdat  = p_dataread->GetValue("BEAM_DATA_FILE",std::string("Beam.dat"));
  m_isrdat   = p_dataread->GetValue("ISR_DATA_FILE",std::string("ISR.dat"));
  m_medat    = p_dataread->GetValue("ME_DATA_FILE",std::string("ME.dat"));

  p_model = NULL;
}


Environment::~Environment() {
  if (p_model)       { delete p_model;       p_model       = NULL; }
  if (p_beamspectra) { delete p_beamspectra; p_beamspectra = NULL; }
  if (p_isrhandler)  { delete p_isrhandler;  p_isrhandler  = NULL; }
}

void Environment::InitializeTheEnvironment() {
  ATOOLS::ParticleInit(m_path); 
  rpa.Init(m_path,m_file);

  bool okay =         InitializeTheModel();  
  okay      = okay && InitializeTheBeams();  
  okay      = okay && InitializeThePDFs();  
}



bool Environment::InitializeTheModel()
{
  msg.Debugging()<<"Initialized Model_Initialization for "<<m_path<<m_modeldat<<std::endl;
  Data_Read     * dataread     = new Data_Read(m_path+m_modeldat);
  Model_Handler * modelhandler = new MODEL::Model_Handler();
  p_model                      = modelhandler->GetModel(dataread,m_path,m_modeldat);

  if (!p_model->RunSpectrumGenerator()) {
    msg.Error()<<"Error in Model_Initialization::Model_Initialization."<<std::endl
	       <<"    RunSpectrumGenerator() delivered false. Abort()."<<std::endl;
    abort();
  }
  delete modelhandler;
  delete dataread;  

  return 1;

}

bool Environment::InitializeTheBeams() 
{
  msg.Debugging()<<"Initialized Beam_Initialization for "<<m_path<<m_beamdat<<std::endl;
  Data_Read * dataread = new Data_Read(m_path+m_beamdat);
  p_beamspectra        = new Beam_Spectra_Handler(dataread);
  for (short int i=0;i<2;i++) m_beam_particles[i] = p_beamspectra->GetBeam(i)->Beam();
  delete dataread;  
  
  return 1;

}

bool Environment::InitializeThePDFs() 
{
  msg.Debugging()<<"Initialize ISR_Initialization for "<<m_path<<m_isrdat<<std::endl;
  Data_Read * dataread     = new Data_Read(m_path+m_isrdat);
  p_isrhandler             = NULL;
  PDF_Handler * pdfhandler = new PDF_Handler();
  PDF_Base *  pdfbase;
  ISR_Base ** isrbases     = new ISR_Base*[2];
  Flavour bunch_particles[2];
  double  bunch_splimits[2];

  for (int i=0;i<2;++i) {
    pdfbase = pdfhandler->GetPDFLib(dataread,bunch_particles[i],i);
    if (pdfbase==NULL) {
      msg.Debugging()<<"No ISR for beam "<<i+1<<" : Initialize Intact for "<<bunch_particles[i]<<std::endl;
      isrbases[i]          = new Intact(bunch_particles[i]);     
    }
    else {
      msg.Debugging()<<"ISR for beam "<<i+1<<" : Initialize SF for "<<bunch_particles[i]<<std::endl;
      isrbases[i]          = new Structure_Function(pdfbase,bunch_particles[i]);
    }
  }
  m_bunch_splimits[0] = dataread->GetValue<double>("ISR_SMIN",0.);
  m_bunch_splimits[1] = dataread->GetValue<double>("ISR_SMAX",1.);
  double kplimits[2];
  kplimits[0] = dataread->GetValue<double>("ISR_KPMIN",m_bunch_splimits[0]);
  kplimits[1] = dataread->GetValue<double>("ISR_KPMAX",m_bunch_splimits[1]);
  p_isrhandler = new ISR_Handler(isrbases,m_bunch_splimits,kplimits);
  delete dataread;

  if (!(p_beamspectra->CheckConsistency(bunch_particles))) {
    msg.Error()<<"Error in Environment::InitializeThePDFs()"<<std::endl
	       <<"   Inconsistent ISR & Beam:"<<std::endl
	       <<"   Abort program."<<std::endl;
    abort();
  }
  return 1;  
}
