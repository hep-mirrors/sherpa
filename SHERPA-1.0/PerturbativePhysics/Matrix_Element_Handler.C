#include"Matrix_Element_Handler.H"

#include "Data_Read.H"
#include "Message.H"

using namespace SHERPA;
using namespace MODEL;
using namespace BEAM;
using namespace ISR;
using namespace PDF;
using namespace AORGTOOLS;
using namespace APHYTOOLS;
using namespace AMATOOLS;
using namespace std;

Matrix_Element_Handler::Matrix_Element_Handler(std::string _dir,std::string _file,
					       MODEL::Model_Base * _model,
					       BEAM::Beam_Spectra_Handler * _beam,
					       ISR::ISR_Handler * _isr) :
  p_amegic(NULL), p_simplexs(NULL), m_dir(_dir), m_file(_file), m_mode(0)
{
  p_dataread = new Data_Read(m_dir+m_file);
  m_signalgenerator  = p_dataread->GetValue<string>("ME_SIGNAL_GENERATOR",std::string("Amegic"));
  if (m_signalgenerator==string("Amegic"))   m_mode = InitializeAmegic(_model,_beam,_isr);
  if (m_signalgenerator==string("Internal")) m_mode = InitializeSimpleXS(_model,_beam,_isr);

  if (p_dataread->GetValue<string>("EVENT_GENERATION_MODE",std::string("Unweighted"))==string("Unweighted"))
    m_eventmode=1;

  msg.Debugging()<<"Run Matrix_Element_Handler in mode :"<<m_mode
		 <<" and event generation mode : "<<m_eventmode<<endl;
  if (m_mode>0) return;

  msg.Error()<<"Error in Matrix_Element_Handler::Matrix_Element_Handler."<<endl
	     <<"   Failed to initialize "<<m_signalgenerator<<" for hard interactions."<<endl
	     <<"   will abort the run."<<endl;
  abort();
}

Matrix_Element_Handler::~Matrix_Element_Handler()
{
  if (p_dataread) { delete p_dataread; p_dataread = NULL; }
  if (p_amegic)   { delete p_amegic;   p_amegic   = NULL; }
  if (p_simplexs) { delete p_simplexs; p_simplexs = NULL; }
}


int Matrix_Element_Handler::InitializeAmegic(MODEL::Model_Base * _model,
					      BEAM::Beam_Spectra_Handler * _beam,
					      ISR::ISR_Handler * _isr) 
{
  m_name    = string("Amegic");
  p_amegic  = new AMEGIC::Amegic(m_dir,m_file,_model);
  if (p_amegic->InitializeProcesses(_beam,_isr)) return 1;
  return 0;
}

int Matrix_Element_Handler::InitializeSimpleXS(MODEL::Model_Base * _model,
						BEAM::Beam_Spectra_Handler * _beam,
						ISR::ISR_Handler * _isr) 
{
  m_name    = string("Simple X-section");
  p_simplexs = new EXTRAXS::SimpleXSecs(m_dir,m_file,_model);
  if (p_simplexs->InitializeProcesses(_beam,_isr)) return 2;
  return 0;
}

bool Matrix_Element_Handler::CalculateTotalXSecs() 
{
  switch (m_mode) {
  case 1: 
    m_readin = p_dataread->GetValue<string>("RESULT DIRECTORY",string(""));
    cout<<" fac "<<rpa.test.FactorYcut()<<endl;
    cout<<" ycut "<<(rpa.integ.Ycut()*sqr(rpa.gen.Ecms()))<<endl;
    p_amegic->Processes()->SetScale(rpa.test.FactorYcut()*rpa.integ.Ycut()*sqr(rpa.gen.Ecms()));
    if (p_amegic->CalculateTotalXSec(m_readin)) {
      RescaleJetrates();
      return 1;
    }
    msg.Error()<<"Error in Matrix_Element_Handler::CalculateTotalXSecs()."<<endl
	       <<"   Failed to Calculate total XSec through Amegic. Abort."<<endl;
    abort();
  case 2:
    if (p_simplexs->CalculateTotalXSec()) return 1;
    msg.Error()<<"Error in Matrix_Element_Handler::CalculateTotalXSecs()."<<endl
	       <<"   Failed to Calculate total XSec through SimpleXS. Abort."<<endl;
    abort();
  }
  msg.Error()<<"Error in Matrix_Element_Handler::CalculateTotalXSecs()."<<endl
	     <<"   Failed to Calculate total XSec. m_mode = "<<m_mode<<" Abort."<<endl;
  abort();
}

bool Matrix_Element_Handler::RescaleJetrates() 
{
  // processes not rescaled in the moment only status printed
  AMEGIC::Process_Base * procs = p_amegic->Processes();
  for (int i=0; i<procs->Size();++i) {
    double xstot = (*procs)[i]->Total()*rpa.Picobarn();
    double njet  = (*procs)[i]->Nout();
    cout<<" "<<njet<<" : "<<xstot<<endl;
  }
}

bool Matrix_Element_Handler::PrepareXSecTables() {};
bool Matrix_Element_Handler::LookUpXSec(double,bool,std::string) {};


bool Matrix_Element_Handler::GenerateOneEvent() 
{
  if (m_eventmode) return UnweightedEvent();
  return WeightedEvent();
}

bool Matrix_Element_Handler::GenerateSameEvent() {}

bool Matrix_Element_Handler::UnweightedEvent() 
{
  switch (m_mode) {
  case 1: return p_amegic->UnweightedEvent();
  case 2: return p_simplexs->UnweightedEvent();
  }
  return 0;
}



bool Matrix_Element_Handler::WeightedEvent() {};

int Matrix_Element_Handler::MaxJets() {
  switch (m_mode) {
  case 1: return p_amegic->MaxJets();
  case 2: return p_simplexs->MaxJets();
  }
  return 0;
}


std::string Matrix_Element_Handler::Name() { return m_name; }

int Matrix_Element_Handler::Nin() {
  switch (m_mode) {
  case 1: return p_amegic->Nin();
  case 2: return p_simplexs->Nin();
  }
  return 0;
}

int Matrix_Element_Handler::Nout() {
  switch (m_mode) {
  case 1: return p_amegic->Nout();
  case 2: return p_simplexs->Nout();
  }
  return 0;
}

std::string Matrix_Element_Handler::SignalGenerator() { return m_signalgenerator; }

std::string Matrix_Element_Handler::ProcessName() 
{
  switch (m_mode) {
  case 1: return p_amegic->ProcessName();
  case 2: return p_simplexs->ProcessName();
  }
  return string("Error - no process name!");
}

AMATOOLS::Vec4D * Matrix_Element_Handler::Momenta() {
  switch (m_mode) {
  case 1: return p_amegic->Momenta();
  case 2: return p_simplexs->Momenta();
  }
  return NULL;
}

APHYTOOLS::Flavour * Matrix_Element_Handler::Flavs() {
  switch (m_mode) {
  case 1: return p_amegic->Flavs();
  case 2: return p_simplexs->Flavs();
  }
  return NULL;
}


int Matrix_Element_Handler::NumberOfDiagrams() 
{
  if (m_mode==1) return p_amegic->NumberOfDiagrams();
  msg.Error()<<"Error in Matrix_Element_Handler::NumberOfDiagrams()."<<endl
	     <<"   Run in mode for "<<m_signalgenerator<<", abort."<<endl;
  abort();
}

AMEGIC::Point * Matrix_Element_Handler::GetDiagram(int _diag) 
{
  if (m_mode==1) return p_amegic->Diagram(_diag);
  msg.Error()<<"Error in Matrix_Element_Handler::GetDiagram("<<_diag<<")."<<endl
	     <<"   Run in mode for "<<m_signalgenerator<<", abort."<<endl;
  abort();
}


EXTRAXS::XS_Base * Matrix_Element_Handler::GetXS() 
{
  if (m_mode==2) return p_simplexs->Selected();
  msg.Error()<<"Error in Matrix_Element_Handler::GetXS()."<<endl
	     <<"   Run in mode for "<<m_signalgenerator<<", abort."<<endl;
  abort();
}

