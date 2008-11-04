#include "Matrix_Element_Handler.H"
#include "Message.H"
#include "Amegic.H"
#include "Simple_XS.H"
#include "Random.H"
#include "Exception.H"
#include "Spin_Structure.H"
#include "Decay_Table.H"
#include <iomanip>

using namespace SHERPA;
using namespace MODEL;
using namespace BEAM;
using namespace PDF;
using namespace ATOOLS;
using namespace std;

// ============================================================
//                Particle_Map
// ============================================================

Particle_Map::Particle_Map():
  m_flip_anti(false), m_change_flavs(false)
{
}

bool Particle_Map::Unity() 
{
  return !(m_flip_anti || m_change_flavs);
}

void Particle_Map::SetFlipAnti() 
{
  m_flip_anti=true;
}

void Particle_Map::Add(const ATOOLS::Flavour & a,const ATOOLS::Flavour & b)
{
  m_flmap[a]=b;
}

bool Particle_Map::Apply(int n, ATOOLS::Flavour * flavs, ATOOLS::Vec4D * moms)
{
  for (int i=0; i<n; ++i) {
    Flavour_Map::iterator it = m_flmap.find(flavs[i]);
    if (it!=m_flmap.end()) {
      flavs[i] = it->second;
    }
  }

  if (m_flip_anti) {
    for (int i=0; i<n; ++i) {
      flavs[i] = flavs[i].Bar();
      moms[i]  = Vec4D(moms[i][0],-1.*Vec3D(moms[i]));
      // check that momenta in lab system !!
      
    }
    Flavour help = flavs[0];
    flavs[0] = flavs[1];
    flavs[1] = help;
    Vec4D mom = moms[0];
    moms[0] = moms[1];
    moms[1] = mom;
  }
  return m_flip_anti;
}



// ============================================================
//               Matrix_Element_Handler
// ============================================================

Matrix_Element_Handler::Matrix_Element_Handler() :
  m_dir("./"), m_file(""), p_amegic(NULL), p_simplexs(NULL),
  p_isr(NULL), p_model(NULL), m_mode(0), m_weight(1.), m_ntrial(1), m_xsecntrial(0),
  m_sntrial(0), m_name(""), m_eventmode(1), 
  m_sudakovon(0), m_apply_hhmf(0), m_ini_swaped(0), p_reader(NULL), p_flavs(NULL), p_moms(NULL) {}

Matrix_Element_Handler::Matrix_Element_Handler(std::string _dir,std::string _file,
					       MODEL::Model_Base * _model,
					       Matrix_Element_Handler * _me) :
  m_dir(_dir), m_file(_file), p_amegic(NULL), p_simplexs(NULL),
  p_isr(NULL), p_model(_model), m_mode(0), m_weight(1.), m_ntrial(1), m_xsecntrial(0), 
  m_sntrial(0), m_name(""), m_eventmode(1),
  m_sudakovon(0), m_apply_hhmf(0), m_ini_swaped(0), p_reader(NULL), p_flavs(NULL), p_moms(NULL) 
{
  if (_me) p_amegic = _me->GetAmegic(); 
  m_mode      = InitializeAmegic(_model,NULL,NULL);
  if (m_mode>0) return;
  THROW(normal_exit,"Failed to initialize "+m_signalgenerator
	+" for hard interactions.");
}

Matrix_Element_Handler::Matrix_Element_Handler(std::string _dir,std::string _file,
					       MODEL::Model_Base * _model,
					       BEAM::Beam_Spectra_Handler * _beam,
					       PDF::ISR_Handler * _isr,
					       Matrix_Element_Handler * _me) :
  m_dir(_dir), m_file(_file), p_amegic(NULL), p_simplexs(NULL),
  p_isr(_isr), p_model(_model), m_mode(0), m_weight(1.), m_ntrial(1), m_xsecntrial(0),
  m_sntrial(0), m_sudakovon(0), m_apply_hhmf(0),
  m_ini_swaped(0), p_flavs(NULL), p_moms(NULL)
{
  p_reader = new Data_Reader(" ",";","!","=");
  p_reader->SetInputPath(m_dir);
  p_reader->SetInputFile(m_file);
  if (!p_reader->ReadFromFile(m_signalgenerator,"ME_SIGNAL_GENERATOR")) m_signalgenerator="Amegic";
  if (!p_reader->ReadFromFile(m_sudakovon,"SUDAKOV_WEIGHT")) m_sudakovon=1;
  if (!p_reader->ReadFromFile(m_apply_hhmf,"TEVATRON_WpWm")) m_apply_hhmf=0;
  if (m_signalgenerator==string("None")) {
    m_mode=0;
    m_name="None";
    return;
  }
  if (m_signalgenerator==string("Amegic")) {
    if (_me) p_amegic = _me->GetAmegic(); 
    m_mode = InitializeAmegic(_model,_beam,_isr);
  }
  if (m_signalgenerator==string("Internal")) m_mode = InitializeSimpleXS(_model,_beam,_isr);

  string evtm;
  if (!p_reader->ReadFromFile(evtm,"EVENT_GENERATION_MODE")) evtm="Unweighted";
  if (evtm==string("Unweighted"))
    m_eventmode=1;
  else {
    if (evtm==string("WeightedNS"))
      m_eventmode=-31;
    else m_eventmode=0;
  }
  m_nmoms = MaxJets();
  p_flavs = new Flavour[MaxJets()+2];
  p_moms  = new Vec4D[MaxJets()+2];
  if (m_apply_hhmf) SetupHHMF();
  if (m_mode>0) return;
  THROW(normal_exit,"Failed to initialize "+m_signalgenerator
	+" for hard interactions.");
}

Matrix_Element_Handler::~Matrix_Element_Handler()
{
  if (p_moms)  delete [] p_moms;
  if (p_flavs) delete [] p_flavs;
  if (p_reader) delete p_reader;
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
  p_simplexs = new EXTRAXS::Simple_XS(m_dir,m_file,_model);
  if (p_simplexs->InitializeProcesses(_beam,_isr)) return 2;
  return 0;
}

bool Matrix_Element_Handler::AddToDecays(const ATOOLS::Flavour & _flav) 
{
  switch (m_mode) {
    //case 1 : return p_amegic->GetAllDecays()->AddToDecays(_flav);
  }
  msg_Error()<<"Error in Matrix_Element_Handler::AddToDecays("<<_flav<<") : "<<endl
	     <<"   m_mode = "<<m_mode<<" Abort."<<endl;
  abort();
}

bool Matrix_Element_Handler::AddToDecays(ATOOLS::Decay_Channel * _dec)
{
  switch (m_mode) {
    //case 1 : return p_amegic->GetAllDecays()->AddToDecays(_dec);
  }
  msg_Error()<<"Error in Matrix_Element_Handler::AddToDecays(";_dec->Output();
  msg_Error()<<"   m_mode = "<<m_mode<<" Abort."<<endl;
  abort();

}


bool Matrix_Element_Handler::InitializeDecayTables()
{
  switch (m_mode) {
    //  case 1 : 
    //return p_amegic->GetAllDecays()->InitializeDecayTables();
  }
  msg_Error()<<"Error in Matrix_Element_Handler::InitializeDecayTables() : "<<endl
	     <<"   m_mode = "<<m_mode<<" Abort."<<endl;
  abort();
}

bool Matrix_Element_Handler::CalculateWidths() 
{
  switch (m_mode) {
    //case 1: 
    //return p_amegic->GetAllDecays()->CalculateWidths();
  }
  msg_Error()<<"Error in Matrix_Element_Handler::CalculateWidths() : "<<endl
	     <<"   m_mode = "<<m_mode<<" Abort."<<endl;
  abort();
}

bool Matrix_Element_Handler::FillDecayTable(ATOOLS::Decay_Table * _dt,bool _ow) 
{
  switch (m_mode) {
    //case 1: 
    //AMEGIC::Full_Decay_Table * fdt;
    //fdt = p_amegic->GetAllDecays()->GetFullDecayTable(_dt->Flav());
    //for (int i=0;i<fdt->NumberOfChannels();i++) _dt->AddDecayChannel(fdt->GetChannel(i));
    //if (_ow) _dt->Flav().SetWidth(fdt->Width());
    //return 1;
  }
  msg_Error()<<"Error in Matrix_Element_Handler::FillDecayTable() : "<<endl
	     <<"   m_mode = "<<m_mode<<" Abort."<<endl;
  abort();
}


bool Matrix_Element_Handler::CalculateTotalXSecs(int scalechoice) 
{
  if (!p_reader->ReadFromFile(m_readin,"RESULT_DIRECTORY")) m_readin="./Results";
  switch (m_mode) { 
  case 0:
    return 1;
    break;
  case 1: 
    if (p_amegic->CalculateTotalXSec(m_readin,m_eventmode<0?-1:1)) {
      PrintTotalXSec();
      return 1;
    }
    msg_Error()<<"Error in Matrix_Element_Handler::CalculateTotalXSecs()."<<endl
	       <<"   Failed to Calculate total XSec through Amegic. Abort."<<endl;
    abort();
  case 2:
    if (p_simplexs->CalculateTotalXSec(m_readin)) return 1;
    msg_Error()<<"Error in Matrix_Element_Handler::CalculateTotalXSecs()."<<endl
	       <<"   Failed to Calculate total XSec through SimpleXS. Abort."<<endl;
    abort();
  }
  msg_Error()<<"Error in Matrix_Element_Handler::CalculateTotalXSecs()."<<endl
	     <<"   Failed to Calculate total XSec. m_mode = "<<m_mode<<" Abort."<<endl;
  abort();
}

double Matrix_Element_Handler::ExpectedEvents() const
{
  switch (m_mode) { 
  case 1: 
    return p_amegic->Processes()->ExpectedEvents();
  case 2:
    return p_simplexs->ExpectedEvents();
  }
  return 0.0;
}

bool Matrix_Element_Handler::PrintTotalXSec() 
{
  // processes not rescaled at the moment only status printed
  AMEGIC::Process_Base * procs = p_amegic->Processes();

  double errsum=0;
  for (size_t i=0; i<procs->Size();++i) {
    errsum+= (*procs)[i]->TotalError();
  }
  if (errsum!=0.) {
    MyStrStream sstr;
    int ecms = int(rpa.gen.Ecms()*10.);
    sstr<<"xsections_"<<ecms<<".dat"<<endl;
    std::string filename;
    sstr>>filename;
    std::ofstream  rfile(filename.c_str(),std::ios::app);
    rfile<<"# ";
    for (size_t i=0; i<procs->Size();++i) {
      rfile<<(*procs)[i]->Name()<<" ";
    }
    rfile<<endl;
    

    rfile.precision(6);
    rfile<<setw(30)<<rpa.gen.Variable("Y_CUT")<<" ";

    for (size_t i=0; i<procs->Size();++i) {
      double xstot = (*procs)[i]->TotalXS()/((*procs)[i]->EnhanceFactor())*rpa.Picobarn();
      double xserr = (*procs)[i]->TotalError()*rpa.Picobarn();
      //      double njet  = (*procs)[i]->Nout();
      rfile<<setw(10)<<xstot<<" "<<setw(10)<<xserr<<" ";
    }
    rfile<<endl;
    rfile.close();
  }
  return true;
}

bool Matrix_Element_Handler::PrepareXSecTables() { return true; }
bool Matrix_Element_Handler::LookUpXSec(double,bool,std::string) { return true; }


bool Matrix_Element_Handler::GenerateOneEvent() 
{
  if (m_eventmode==1) {
    Blob_Data_Base * message = UnweightedEvent();
    if (message) {
      PHASIC::Weight_Info winfo = message->Get<PHASIC::Weight_Info>();
      m_xf1=winfo.xf1;
      m_xf2=winfo.xf2;
      m_weight =  winfo.weight;
      m_procweight = winfo.procweight;
      m_xsecweight = winfo.xsecweight;
      m_xsecntrial = winfo.xsecntrial;
      m_ntrial =  winfo.ntrial;
      msg_Debugging()<<"MEH::GOE: "<<m_procweight<<" "
		     <<m_weight<<" "<<m_ntrial<<"\n";
      delete message;
    }
    else THROW(fatal_error,"No weight information.");
    GetMomentaNFlavours();
    ApplyHHMF();
    return m_weight>0.0;
  }
  Blob_Data_Base * message = WeightedEvent();
  if (GetMomentaNFlavours()) {
    ApplyHHMF();
  }
  if (message) {
    PHASIC::Weight_Info winfo = message->Get<PHASIC::Weight_Info>();
    m_xf1=winfo.xf1;
    m_xf2=winfo.xf2;
    m_weight =  winfo.weight * rpa.Picobarn();
    m_procweight = winfo.procweight;
    m_xsecweight = winfo.xsecweight;
    m_xsecntrial = winfo.xsecntrial;
    m_ntrial =  winfo.ntrial;
    //PRINT_INFO(m_procweight<<" "<<m_weight<<" "<<m_ntrial);
    delete message;
  }
  else {
    m_weight=0.;
    if (rpa.gen.NumberOfDicedEvents()==rpa.gen.NumberOfEvents()) return true;
  }
  return (m_weight>0.);
}

bool Matrix_Element_Handler::GenerateOneEvent(ATOOLS::Decay_Channel * _dc,double _mass) 
{
  switch (m_mode) {
    //case 1: return p_amegic->GetAllDecays()->UnweightedEvent(_dc,_mass);
  }
  return false;
}

double Matrix_Element_Handler::FactorisationScale()
{
  switch (m_mode) {
  case 1: return p_amegic->GetProcess()->Scale(PHASIC::stp::fac);
  case 2: return p_simplexs->Selected()->Scale(PHASIC::stp::fac);
  }
  return 0.;
}

bool Matrix_Element_Handler::GenerateSameEvent() 
{
  if (m_eventmode==1) {
    Blob_Data_Base * message = UnweightedSameEvent();
    if (message) {
      PHASIC::Weight_Info winfo = message->Get<PHASIC::Weight_Info>();
      m_xf1=winfo.xf1;
      m_xf2=winfo.xf2;
      m_weight     =  winfo.weight;
      m_procweight = winfo.procweight;
      m_xsecweight = winfo.xsecweight;
      m_xsecntrial = winfo.xsecntrial;
      m_ntrial     =  winfo.ntrial;
      delete message;
    }
    else THROW(fatal_error,"No weight information.");
    GetMomentaNFlavours();
    ApplyHHMF();
    return m_weight>0.0;
  }
  Blob_Data_Base * message = WeightedSameEvent();
  if (GetMomentaNFlavours()) {
    ApplyHHMF();
  }
  if (message) {
    PHASIC::Weight_Info winfo = message->Get<PHASIC::Weight_Info>();
    m_xf1=winfo.xf1;
    m_xf2=winfo.xf2;
    m_weight =  winfo.weight * rpa.Picobarn();
    m_ntrial =  winfo.ntrial;
    delete message;
  }
  else {
    m_weight=0.;
    if (rpa.gen.NumberOfDicedEvents()==rpa.gen.NumberOfEvents()) return true;
  }
  return (m_weight>0.);
}

bool Matrix_Element_Handler::GetMomentaNFlavours() {    
  switch (m_mode) {
  case 1: 
    m_nmoms = p_amegic->NIn()+p_amegic->NOut();
    if (p_amegic->Flavours() && p_amegic->Momenta()) {
      for (int i=0;i<m_nmoms;++i) p_flavs[i] = p_amegic->Flavours()[i];
      for (int i=0;i<m_nmoms;++i) p_moms[i]  = p_amegic->Momenta()[i];
      return true;
    }
    break;
  case 2: 
    int nmoms(m_nmoms);
    m_nmoms = p_simplexs->NIn()+p_simplexs->NOut();
    if (m_nmoms>nmoms) {
      delete [] p_flavs;
      delete [] p_moms;
      p_flavs = new Flavour[m_nmoms];
      p_moms = new Vec4D[m_nmoms];
    }
    if (p_simplexs->Flavours() && p_simplexs->Momenta()) {
      for (int i=0;i<m_nmoms;++i) p_flavs[i] = p_simplexs->Flavours()[i];
      for (int i=0;i<m_nmoms;++i) p_moms[i]  = p_simplexs->Momenta()[i];
      return true;
    }
  }
  if (rpa.gen.NumberOfDicedEvents()==rpa.gen.NumberOfEvents()) return false;
  msg_Error()<<"Warning in Matrix_Element_Handler::SetMomenta()"<<endl
	     <<"   No ME generator available to get momenta from."<<endl
	     <<"   Continue run and hope for the best."<<endl;
  return false;
}


ATOOLS::Blob_Data_Base *Matrix_Element_Handler::UnweightedSameEvent() 
{
  switch (m_mode) {
  case 1: return p_amegic->SameEvent();
  case 2: return p_simplexs->OneEvent();
  }
  return 0;
}

  
Blob_Data_Base * Matrix_Element_Handler::WeightedSameEvent() 
{
  switch (m_mode) {
  case 1: return p_amegic->SameWeightedEvent();
  case 2: return p_simplexs->WeightedEvent();
  }
  return 0;
}

Blob_Data_Base *Matrix_Element_Handler::UnweightedEvent() 
{
  switch (m_mode) {
  case 1: return p_amegic->UnweightedEvent();
  case 2: return p_simplexs->OneEvent();
  }
  return 0;
}



Blob_Data_Base * Matrix_Element_Handler::WeightedEvent() 
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

unsigned int Matrix_Element_Handler::MinQCDJets() {
  switch (m_mode) {
  case 1: return p_amegic->MinQCDJets();
  case 2: return p_simplexs->MinQCDJets();
  }
  return 0;
}

unsigned int Matrix_Element_Handler::MaxQCDJets() {
  switch (m_mode) {
  case 1: return p_amegic->MaxQCDJets();
  case 2: return p_simplexs->MaxQCDJets();
  }
  return 0;
}


std::string Matrix_Element_Handler::Name() { return m_name; }

unsigned int Matrix_Element_Handler::NIn() {
  switch (m_mode) {
  case 1: return p_amegic->NIn();
  case 2: return p_simplexs->NIn();
  }
  return 0;
}

unsigned int Matrix_Element_Handler::NOut() {
  switch (m_mode) {
  case 1: return p_amegic->NOut();
  case 2: return p_simplexs->NOut();
  }
  return 0;
}

unsigned int Matrix_Element_Handler::NDecOut() {
  switch (m_mode) {
    //case 1: return p_amegic->GetAllDecays()->NOut();
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

const ATOOLS::Vec4D * Matrix_Element_Handler::Momenta() {
  switch (m_mode) {
  case 1: return p_amegic->Momenta();
  case 2: return p_simplexs->Momenta();
  }
  return NULL;
}

const ATOOLS::Vec4D * Matrix_Element_Handler::DecMomenta() {
  switch (m_mode) {
    //case 1: return p_amegic->GetAllDecays()->Momenta();
  }
  return NULL;
}

const ATOOLS::Flavour * Matrix_Element_Handler::Flavours() {
  switch (m_mode) {
  case 1: return p_amegic->Flavours();
  case 2: return p_simplexs->Flavours();
  }
  return NULL;
}

void  Matrix_Element_Handler::SetupHHMF() 
{
  // Tevatron :
  //   starting with   W- -> e- neb
  //   generate  W- -> e- nueb
  //             W+ -> e+ nue
  //             W- -> mu- nmub
  //             W+ -> mu+ nmu

  Particle_Map  unity;
  m_particle_maps.push_back(unity);

  Particle_Map  chargeconj;
  chargeconj.SetFlipAnti();
  m_particle_maps.push_back(chargeconj);

  Particle_Map  muon;
  muon.Add(Flavour(kf_e),Flavour(kf_mu));
  muon.Add(Flavour(kf_nue).Bar(),Flavour(kf_numu).Bar());
  m_particle_maps.push_back(muon);

  muon.SetFlipAnti();
  m_particle_maps.push_back(muon);
}

bool  Matrix_Element_Handler::ApplyHHMF()
{
  if (!m_apply_hhmf || m_particle_maps.size()==0) return 1;

  size_t no=(size_t)(m_particle_maps.size()*ATOOLS::ran.Get());
  if (no==m_particle_maps.size()) no=0;

  m_ini_swaped=m_particle_maps[no].Apply(m_nmoms,p_flavs,p_moms);

  return 1;
}


int Matrix_Element_Handler::InSwaped() {
  switch (m_mode) {
  case 1: return m_ini_swaped^p_amegic->InSwaped();
  case 2: return m_ini_swaped^p_simplexs->InSwaped();
  }
  return 0;
}

unsigned int Matrix_Element_Handler::NumberOfDiagrams() 
{
  if (m_mode==1) return p_amegic->NumberOfDiagrams();
  msg_Error()<<"Error in Matrix_Element_Handler::NumberOfDiagrams()."<<endl
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
  msg_Error()<<"Error in Matrix_Element_Handler::GetDiagram("<<_diag<<")."<<endl
	     <<"   Wrong mode for "<<m_signalgenerator<<", abort."<<endl;
  abort();
}

// ATOOLS::Spin_Correlation_Tensor * Matrix_Element_Handler::GetSpinCorrelations() 
// {
//   if (Spin_Correlation_Tensor::Mode()==scmode::None) return NULL;
//   if (m_mode==1) return p_amegic->GetProcess()->GetSpinCorrelations();
//   msg_Error()<<"Error in Matrix_Element_Handler::GetSpinCorrelations()."<<endl
// 	     <<"   Wrong mode for ME generator "<<m_signalgenerator<<", abort."<<endl;
//   abort();
// }

void Matrix_Element_Handler::FillAmplitudes(HELICITIES::Amplitude_Tensor* atensor) 
{
  p_amegic->GetProcess()->FillAmplitudes(atensor);
}

bool Matrix_Element_Handler::SpinCorrelations()
{
  return m_spincorrelations;
}

void Matrix_Element_Handler::SetSpinCorrelations(bool sc)
{
  m_spincorrelations = sc;
}

EXTRAXS::XS_Base * Matrix_Element_Handler::GetXS(const int mode) 
{
  if (m_mode==2) {
    if (p_simplexs->Selected()!=NULL) {
      return static_cast<EXTRAXS::XS_Base*>(p_simplexs->Selected());
    }
    else {
      return p_simplexs;
    }
  }
  if (mode>0) return NULL;
  msg_Error()<<"Error in Matrix_Element_Handler::GetXS("<<mode<<")."<<endl
	     <<"   Wrong mode for "<<m_signalgenerator<<", abort."<<endl;
  abort();
}

double  Matrix_Element_Handler::Weight() 
{
  return m_weight;
}

double  Matrix_Element_Handler::ProcessWeight() 
{
  return m_procweight;
}

double  Matrix_Element_Handler::XSecWeight() 
{
  return m_xsecweight;
}

unsigned long Matrix_Element_Handler::NumberOfTrials()
{
  return m_ntrial+m_sntrial;
}

unsigned long Matrix_Element_Handler::NumberOfXSecTrials()
{
  return m_xsecntrial;
}

int Matrix_Element_Handler::OrderStrong()
{
  if (m_mode==1) return p_amegic->OrderStrong();
  msg_Error()<<"Error in Matrix_Element_Handler::OrderStrong()."<<endl
	     <<"   Wrong mode for "<<m_signalgenerator<<", abort."<<endl;
  abort();
}

int Matrix_Element_Handler::OrderEWeak()
{
  if (m_mode==1) return p_amegic->OrderEWeak();
  msg_Error()<<"Error in Matrix_Element_Handler::OrderStrong()."<<endl
	     <<"   Wrong mode for "<<m_signalgenerator<<", abort."<<endl;
  abort();
}

void Matrix_Element_Handler::SetAmegic(AMEGIC::Amegic *amegic)
{
  if (p_amegic!=NULL || p_simplexs!=NULL || amegic==NULL) return;
  p_amegic=amegic;
  p_isr=p_amegic->Processes()->ISR();
  m_name=string("Amegic");
  m_mode=1;
}

void Matrix_Element_Handler::SetXS(EXTRAXS::Simple_XS *simplexs)
{
  if (p_amegic!=NULL || p_simplexs!=NULL || simplexs==NULL) return;
  p_simplexs=simplexs;
  p_isr=p_simplexs->ISR();
  m_name=string("SimpleXS");
  m_mode=2;
}

