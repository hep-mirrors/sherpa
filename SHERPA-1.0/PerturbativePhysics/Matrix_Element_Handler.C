#include"Matrix_Element_Handler.H"
#include "Data_Read.H"
#include "Message.H"
#include "Amegic.H"
#include "SimpleXSecs.H"
#include <iomanip>

using namespace SHERPA;
using namespace MODEL;
using namespace BEAM;
using namespace PDF;
using namespace ATOOLS;
using namespace std;

Matrix_Element_Handler::Matrix_Element_Handler() :
  m_dir("./"), m_file(""), p_amegic(NULL), p_simplexs(NULL),
  p_isr(NULL), m_mode(0), m_weight(1.), m_name(""), m_eventmode(1),
  p_dataread(NULL) {}

Matrix_Element_Handler::Matrix_Element_Handler(std::string _dir,std::string _file,
					       MODEL::Model_Base * _model,
					       Matrix_Element_Handler * _me) :
  m_dir(_dir), m_file(_file), p_amegic(NULL), p_simplexs(NULL),
  p_isr(NULL), m_mode(0), m_weight(1.), m_name(""), m_eventmode(1),
  p_dataread(NULL) 
{
  if (_me) p_amegic = _me->GetAmegic(); 
  m_mode      = InitializeAmegic(_model,NULL,NULL);
  msg.Debugging()<<"Run Matrix_Element_Handler in mode :"<<m_mode
		 <<" and event generation mode : "<<m_eventmode<<endl;
  if (m_mode>0) return;

  msg.Error()<<"Error in Matrix_Element_Handler::Matrix_Element_Handler."<<endl
	     <<"   Failed to initialize "<<m_signalgenerator<<" for hard interactions."<<endl
	     <<"   will abort the run."<<endl;
  abort();
}

Matrix_Element_Handler::Matrix_Element_Handler(std::string _dir,std::string _file,
					       MODEL::Model_Base * _model,
					       BEAM::Beam_Spectra_Handler * _beam,
					       PDF::ISR_Handler * _isr,
					       Matrix_Element_Handler * _me) :
  m_dir(_dir), m_file(_file), p_amegic(NULL), p_simplexs(NULL),
  p_isr(_isr), m_mode(0), m_weight(1.)
{
  p_dataread        = new Data_Read(m_dir+m_file);
  m_signalgenerator = p_dataread->GetValue<string>("ME_SIGNAL_GENERATOR",std::string("Amegic"));
  m_sudakovon       = p_dataread->GetValue<int>("SUDAKOV WEIGHT",0);
  if (m_signalgenerator==string("Amegic")) {
    if (_me) p_amegic = _me->GetAmegic(); 
    m_mode = InitializeAmegic(_model,_beam,_isr);
  }
  if (m_signalgenerator==string("Internal")) m_mode = InitializeSimpleXS(_model,_beam,_isr);

  if (p_dataread->GetValue<string>("EVENT_GENERATION_MODE",std::string("Unweighted"))==string("Unweighted"))
    m_eventmode=1;
  else
    m_eventmode=0;

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
					     PDF::ISR_Handler * _isr) 
{
  m_name    = string("Amegic");
  if (!p_amegic) p_amegic = new AMEGIC::Amegic(m_dir,m_file,_model);
  if (_beam!=NULL || _isr!=NULL) {
    if (p_amegic->InitializeProcesses(_beam,_isr)) return 1;
  }
  else {
    if (p_amegic->InitializeDecays(0)) return 1;
  }
  return 0;
}

int Matrix_Element_Handler::InitializeSimpleXS(MODEL::Model_Base * _model,
					       BEAM::Beam_Spectra_Handler * _beam,
					       PDF::ISR_Handler * _isr) 
{
  m_name     = string("SimpleXS");
  p_simplexs = new EXTRAXS::SimpleXSecs(m_dir,m_file,_model);
  if (p_simplexs->InitializeProcesses(_beam,_isr)) return 2;
  return 0;
}

bool Matrix_Element_Handler::AddToDecays(ATOOLS::Flavour & _flav) 
{
  switch (m_mode) {
  case 1 : return p_amegic->GetAllDecays()->AddToDecays(_flav);
  }
  msg.Error()<<"Error in Matrix_Element_Handler::AddToDecays("<<_flav<<") : "<<endl
	     <<"   m_mode = "<<m_mode<<" Abort."<<endl;
  abort();

}

bool Matrix_Element_Handler::InitializeDecayTables()
{
  switch (m_mode) {
  case 1 : return p_amegic->GetAllDecays()->InitializeDecayTables();
  }
  msg.Error()<<"Error in Matrix_Element_Handler::InitializeDecayTables() : "<<endl
	     <<"   m_mode = "<<m_mode<<" Abort."<<endl;
  abort();
}

bool Matrix_Element_Handler::CalculateWidths() 
{
  switch (m_mode) {
  case 1: 
    return p_amegic->GetAllDecays()->CalculateWidths();
  }
  msg.Error()<<"Error in Matrix_Element_Handler::CalculateWidths() : "<<endl
	     <<"   m_mode = "<<m_mode<<" Abort."<<endl;
  abort();
}

bool Matrix_Element_Handler::FillDecayTable(ATOOLS::Decay_Table * _dt,bool _ow) 
{
  switch (m_mode) {
  case 1: 
    AMEGIC::Full_Decay_Table * fdt;
    fdt = p_amegic->GetAllDecays()->GetFullDecayTable(_dt->Flav());
    for (int i=0;i<fdt->NumberOfChannels();i++) _dt->AddDecayChannel(fdt->GetChannel(i));
    if (_ow) _dt->Flav().SetWidth(fdt->Width());
    return 1;
  }
  msg.Error()<<"Error in Matrix_Element_Handler::FillDecayTable() : "<<endl
	     <<"   m_mode = "<<m_mode<<" Abort."<<endl;
  abort();
}


bool Matrix_Element_Handler::CalculateTotalXSecs(int scalechoice) 
{
  cout<<" scalechoice "<<scalechoice<<endl;
  switch (m_mode) { 
  case 1: 
    m_readin = p_dataread->GetValue<string>("RESULT DIRECTORY",string(""));
    if (scalechoice>0) p_amegic->Processes()->SetScale(rpa.gen.Ycut()*sqr(rpa.gen.Ecms()));
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

  double errsum=0;
  for (int i=0; i<procs->Size();++i) {
    errsum+= (*procs)[i]->TotalError();
  }

  if (errsum!=0.) {
    MyStrStream sstr;
    int ecms = int(rpa.gen.Ecms()*10.);
    double ycut=log(rpa.gen.Ycut())/log(10.);
    sstr<<"xsections_"<<ecms<<".dat"<<endl;
    std::string filename;
    sstr>>filename;
    cout<<" looking for "<<filename<<endl;
    std::ofstream  rfile(filename.c_str(),std::ios::app);
    rfile<<"# ";
    for (int i=0; i<procs->Size();++i) {
      rfile<<(*procs)[i]->Name()<<" ";
    }
    rfile<<endl;
    

    rfile.precision(6);
    rfile<<setw(10)<<ycut<<" ";

    for (int i=0; i<procs->Size();++i) {
      double xstot = (*procs)[i]->Total()*rpa.Picobarn();
      double xserr = (*procs)[i]->TotalError()*rpa.Picobarn();
      double njet  = (*procs)[i]->Nout();
      rfile<<setw(10)<<xstot<<" "<<setw(10)<<xserr<<" ";
    }
    rfile<<endl;
    rfile.close();
  }
  return true;
}

bool Matrix_Element_Handler::PrepareXSecTables() { return true; };
bool Matrix_Element_Handler::LookUpXSec(double,bool,std::string) { return true; };


bool Matrix_Element_Handler::GenerateOneEvent() 
{
  if (m_eventmode) return UnweightedEvent();
  m_weight = WeightedEvent() * rpa.Picobarn();
  return (m_weight>0.);
}

bool Matrix_Element_Handler::GenerateOneEvent(ATOOLS::Decay_Channel * _dc,double _mass) 
{
  switch (m_mode) {
  case 1: return p_amegic->GetAllDecays()->UnweightedEvent(_dc,_mass);
  }
  return false;
}

bool Matrix_Element_Handler::GenerateSameEvent() 
{
  if (m_eventmode) {
    switch (m_mode) {
    case 1: return p_amegic->SameEvent();
    case 2: return p_simplexs->OneEvent();
    }
  }
  else {
    switch (m_mode) {
    case 1: return p_amegic->SameWeightedEvent();
    case 2: return p_simplexs->WeightedEvent();
    }
  }
  return false;
}

bool Matrix_Element_Handler::UnweightedEvent() 
{
  switch (m_mode) {
  case 1: return p_amegic->UnweightedEvent();
  case 2: return p_simplexs->UnweightedEvent();
  }
  return 0;
}



double Matrix_Element_Handler::WeightedEvent() 
{
  switch (m_mode) {
  case 1: return p_amegic->WeightedEvent();
  case 2: return p_simplexs->WeightedEvent();
  }
  return 0;
}

unsigned int Matrix_Element_Handler::MaxJets() {
  switch (m_mode) {
  case 1: return p_amegic->MaxJets();
  case 2: return p_simplexs->MaxJets();
  }
  return 0;
}


std::string Matrix_Element_Handler::Name() { return m_name; }

unsigned int Matrix_Element_Handler::NIn() {
  switch (m_mode) {
  case 1: return p_amegic->Nin();
  case 2: return p_simplexs->Nin();
  }
  return 0;
}

unsigned int Matrix_Element_Handler::NOut() {
  switch (m_mode) {
  case 1: return p_amegic->Nout();
  case 2: return p_simplexs->Nout();
  }
  return 0;
}

unsigned int Matrix_Element_Handler::NDecOut() {
  switch (m_mode) {
  case 1: return p_amegic->GetAllDecays()->Nout();
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

ATOOLS::Vec4D * Matrix_Element_Handler::Momenta() {
  switch (m_mode) {
  case 1: return p_amegic->Momenta();
  case 2: return p_simplexs->Momenta();
  }
  return NULL;
}

ATOOLS::Vec4D * Matrix_Element_Handler::DecMomenta() {
  switch (m_mode) {
  case 1: return p_amegic->Momenta();
  }
  return NULL;
}

ATOOLS::Flavour * Matrix_Element_Handler::Flavs() {
  switch (m_mode) {
  case 1: return p_amegic->Flavs();
  case 2: return p_simplexs->Flavs();
  }
  return NULL;
}

int Matrix_Element_Handler::InSwaped() {
  switch (m_mode) {
  case 1: return p_amegic->InSwaped();
  case 2: return p_simplexs->InSwaped();
  }
  return 0;
}

unsigned int Matrix_Element_Handler::NumberOfDiagrams() 
{
  if (m_mode==1) return p_amegic->NumberOfDiagrams();
  msg.Error()<<"Error in Matrix_Element_Handler::NumberOfDiagrams()."<<endl
	     <<"   Wrong mode for "<<m_signalgenerator<<", abort."<<endl;
  abort();
}

AMEGIC::Amegic * Matrix_Element_Handler::GetAmegic() 
{
  if (m_mode==1) return p_amegic;
  return NULL;
}

AMEGIC::Point * Matrix_Element_Handler::GetDiagram(int _diag) 
{
  if (m_mode==1) return p_amegic->Diagram(_diag);
  msg.Error()<<"Error in Matrix_Element_Handler::GetDiagram("<<_diag<<")."<<endl
	     <<"   Wrong mode for "<<m_signalgenerator<<", abort."<<endl;
  abort();
}


EXTRAXS::XS_Base * Matrix_Element_Handler::GetXS() 
{
  if (m_mode==2) return p_simplexs->Selected();
  msg.Error()<<"Error in Matrix_Element_Handler::GetXS()."<<endl
	     <<"   Wrong mode for "<<m_signalgenerator<<", abort."<<endl;
  std::cout<<p_simplexs<<std::endl;
  abort();
}

double  Matrix_Element_Handler::Weight() 
{
  return m_weight;
}

void Matrix_Element_Handler::SetAmegic(AMEGIC::Amegic *_p_amegic)
{
  if ((p_amegic!=NULL)||(p_simplexs!=NULL)) return;
  p_amegic=_p_amegic;
  p_isr=p_amegic->Processes()->ISR();
  m_name=string("Amegic");
  m_mode=1;
}

void Matrix_Element_Handler::SetXS(EXTRAXS::SimpleXSecs *_p_simplexs)
{
  if ((p_amegic!=NULL)||(p_simplexs!=NULL)) return;
  p_simplexs=_p_simplexs;
  p_isr=p_simplexs->ISR();
  m_name=string("SimpleXS");
  m_mode=2;
}

