#include "Environment.H"

#include "Model_Handler.H"
#include "PDF_Handler.H"
#include "Structure_Function.H"
#include "Intact.H"
#include "Run_Parameter.H"
#include "Phase_Space_Handler.H"

using namespace AMISIC;

Environment::Environment(const std::string path,const std::string file): 
  m_path(path), m_file(file),
  p_beamhandler(NULL), p_isrhandler(NULL), p_model(NULL)
{
  ATOOLS::msg->Init(2,"");
  ATOOLS::Data_Read dataread(m_path+m_file);
  m_modeldat=dataread.GetValue("MODEL_DATA_FILE",std::string("Model.dat"));
  m_beamdat=dataread.GetValue("BEAM_DATA_FILE",std::string("Beam.dat"));
  m_isrdat=dataread.GetValue("ISR_DATA_FILE",std::string("ISR.dat"));
}

Environment::~Environment() 
{
  if (p_isrhandler) delete p_isrhandler;
  if (p_beamhandler) delete p_beamhandler;
  if (p_model) delete p_model;
}

bool Environment::InitializeTheEnvironment() 
{
  char *dummy="Sherpa";
  ATOOLS::rpa.Init(m_path,m_file,1,&dummy);
  ATOOLS::ParticleInit(m_path); 
  if (!InitializeTheModel()) return false;  
  if (!InitializeTheBeams()) return false;  
  if (!InitializeThePDFs()) return false;
  ATOOLS::Integration_Info *info=PHASIC::Phase_Space_Handler::GetInfo();
  p_isrhandler->AssignKeys(info);
  return true;
}

bool Environment::InitializeTheModel()
{
  ATOOLS::Data_Read dataread(m_path+m_modeldat);
  MODEL::Model_Handler *modelhandler = new MODEL::Model_Handler();
  p_model = modelhandler->GetModel(&dataread,m_path,m_modeldat);
  if (!p_model->RunSpectrumGenerator()) {
    msg_Error()<<"Environment::InitializeTheModel(): "
		       <<"RunSpectrumGenerator() failed. Abort."<<std::endl;
    abort();
  }
  delete modelhandler;
  return true;
}

bool Environment::InitializeTheBeams() 
{
  ATOOLS::Data_Read dataread(m_path+m_beamdat);
  p_beamhandler = new BEAM::Beam_Spectra_Handler(&dataread);
  for (short int i=0;i<2;i++) m_bunch[i]=p_beamhandler->GetBeam(i)->Beam();
  ATOOLS::rpa.gen.SetBeam1(m_beam[0]);
  ATOOLS::rpa.gen.SetBeam2(m_beam[1]);
  return true;
}

bool Environment::InitializeThePDFs() 
{
  ATOOLS::Data_Read dataread(m_path+m_isrdat);
  PDF::PDF_Base *pdfbase;
  PDF::ISR_Base **isrbases = new PDF::ISR_Base*[2];
  PDF::PDF_Handler pdfhandler;
  for (int i=0;i<2;++i) {
    pdfbase = pdfhandler.GetPDFLib(&dataread,m_beam[i],i);
    if (pdfbase==NULL) {
      msg_Info()<<"No ISR for beam "<<i+1
		<<" : Initialize Intact for "<<m_beam[i]<<std::endl;
      isrbases[i] = new PDF::Intact(m_beam[i]);     
    }
    else {
      msg_Info()<<"ISR for beam "<<i+1
		<<" : Initialize SF for "<<m_beam[i]<<std::endl;
      isrbases[i] = new PDF::Structure_Function(pdfbase,m_beam[i]);
    }
    ATOOLS::rpa.gen.SetBunch(m_beam[i],i);
  }
  double splimits[2];
  splimits[0]=dataread.GetValue<double>("ISR_SMIN",0.);
  splimits[1]=dataread.GetValue<double>("ISR_SMAX",1.);
  double kplimits[2];
  kplimits[0]=dataread.GetValue<double>("ISR_KPMIN",ATOOLS::Accu());
  kplimits[1]=dataread.GetValue<double>("ISR_KPMAX",splimits[1]);
  p_isrhandler = new PDF::ISR_Handler(isrbases);
  p_isrhandler->SetBeam(p_beamhandler->GetBeam(0),0);
  p_isrhandler->SetBeam(p_beamhandler->GetBeam(1),1);
  p_isrhandler->Init(splimits,kplimits);
  if (!(p_beamhandler->CheckConsistency(m_beam))) {
    msg_Error()<<"Error in Environment::InitializeThePDFs() \n"
		       <<"   Inconsistent ISR & Beam:"<<std::endl
		       <<"   Abort program."<<std::endl;
    abort();
  }
  return true;  
}
