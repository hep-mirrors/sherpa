#include "Matrix_Element_Handler.H"
#include "Data_Read.H"
#include "Message.H"
#include "Amegic.H"
#include "SimpleXSecs.H"
#include "Random.H"
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
  p_isr(NULL), m_mode(0), m_weight(1.), m_ntrial(1), m_name(""), m_eventmode(1), 
  m_sudakovon(0), m_apply_hhmf(0), m_ini_swaped(0), p_dataread(NULL), p_flavs(NULL), p_moms(NULL) {}

Matrix_Element_Handler::Matrix_Element_Handler(std::string _dir,std::string _file,
					       MODEL::Model_Base * _model,
					       Matrix_Element_Handler * _me) :
  m_dir(_dir), m_file(_file), p_amegic(NULL), p_simplexs(NULL),
  p_isr(NULL), m_mode(0), m_weight(1.), m_ntrial(1), m_name(""), m_eventmode(1),
  m_sudakovon(0), m_apply_hhmf(0), m_ini_swaped(0), p_dataread(NULL), p_flavs(NULL), p_moms(NULL) 
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
  p_isr(_isr), m_mode(0), m_weight(1.), m_ntrial(1), m_sudakovon(0), m_apply_hhmf(0),
  m_ini_swaped(0), p_flavs(NULL), p_moms(NULL)
{
  p_dataread        = new Data_Read(m_dir+m_file);
  m_signalgenerator = p_dataread->GetValue<string>("ME_SIGNAL_GENERATOR",std::string("Amegic"));
  m_sudakovon       = p_dataread->GetValue<int>("SUDAKOV WEIGHT",0);
  m_apply_hhmf      = p_dataread->GetValue<int>("TEVATRON_WpWm",0);
  if (m_signalgenerator==string("Amegic")) {
    if (_me) p_amegic = _me->GetAmegic(); 
    m_mode = InitializeAmegic(_model,_beam,_isr);
  }
  if (m_signalgenerator==string("Internal")) m_mode = InitializeSimpleXS(_model,_beam,_isr);

  if (p_dataread->GetValue<string>("EVENT_GENERATION_MODE",std::string("Unweighted"))==string("Unweighted"))
    m_eventmode=1;
  else
    m_eventmode=0;

  p_flavs = new Flavour[MaxJets()+2];
  p_moms  = new Vec4D[MaxJets()+2];
  if (m_apply_hhmf) SetupHHMF();

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
  if (p_moms)  delete [] p_moms;
  if (p_flavs) delete [] p_flavs;
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

bool Matrix_Element_Handler::AddToDecays(const ATOOLS::Flavour & _flav) 
{
  switch (m_mode) {
  case 1 : return p_amegic->GetAllDecays()->AddToDecays(_flav);
  }
  msg.Error()<<"Error in Matrix_Element_Handler::AddToDecays("<<_flav<<") : "<<endl
	     <<"   m_mode = "<<m_mode<<" Abort."<<endl;
  abort();
}

bool Matrix_Element_Handler::AddToDecays(ATOOLS::Decay_Channel * _dec)
{
  switch (m_mode) {
  case 1 : return p_amegic->GetAllDecays()->AddToDecays(_dec);
  }
  msg.Error()<<"Error in Matrix_Element_Handler::AddToDecays(";_dec->Output();
  msg.Error()<<"   m_mode = "<<m_mode<<" Abort."<<endl;
  abort();

}


bool Matrix_Element_Handler::InitializeDecayTables()
{
  switch (m_mode) {
  case 1 : 
    return p_amegic->GetAllDecays()->InitializeDecayTables();
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
  switch (m_mode) { 
  case 1: 
    m_readin = p_dataread->GetValue<string>("RESULT DIRECTORY",string("./Results"));
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


void RescaleProcesses(AMEGIC::Process_Base * procs, double fac ) {
  if (fac==1.) return;
  if (!procs) return;
  if ((*procs)[0]==procs) {
    double xs=procs->Total();
    procs->SetTotal(xs*fac);
    std::cout<<" changing xs from "<<xs<<" to "<<procs->Total()<<std::endl;
  }
  else {
    double xs=procs->Total();
    for (int i=0; i<procs->Size();++i) {
      RescaleProcesses((*procs)[i],fac);
    }
    procs->SetTotal(xs*fac);
    std::cout<<" changing xs from "<<xs<<" to "<<procs->Total()<<std::endl;
  }
}

bool Matrix_Element_Handler::RescaleJetrates() 
{
  // processes not rescaled in the moment only status printed
  AMEGIC::Process_Base * procs = p_amegic->Processes();

  double errsum=0;
  for (int i=0; i<procs->Size();++i) {
    errsum+= (*procs)[i]->TotalError();
  }

  cout<<" rescale Jetrates : "<<endl;
  //vs.facs[10] = { 1., 1., 1. , 0.1, 1., 1.,1., 1., 1., 1.};
  double facs[10] = { 1., 1., 1. , 1., 1., 1.,1., 1., 1., 1.};
  for (int i=0; i<procs->Size();++i) {
    //    double xstot = (*procs)[i]->Total()*rpa.Picobarn();
    //    double xserr = (*procs)[i]->TotalError()*rpa.Picobarn();
    int njet  = (*procs)[i]->Nout();
    RescaleProcesses((*procs)[i],facs[njet]);
  }
  procs->SetMax(0.);

  if (errsum!=0.) {
    MyStrStream sstr;
    int ecms = int(rpa.gen.Ecms()*10.);
    double ycut=log(rpa.gen.Ycut())/log(10.);
    sstr<<"xsections_"<<ecms<<".dat"<<endl;
    std::string filename;
    sstr>>filename;
    msg.Debugging()<<" looking for "<<filename<<endl;
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
      //      double njet  = (*procs)[i]->Nout();
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
  if (m_eventmode) {
    bool stat = UnweightedEvent();
    if (!stat) return stat;
    GetMomentaNFlavours();
    ApplyHHMF();
    return stat;
  }
  Blob_Data_Base * message = WeightedEvent();
  GetMomentaNFlavours();
  ApplyHHMF();
  if (message) {
    PHASIC::Weight_Info winfo = message->Get<PHASIC::Weight_Info>();
    m_weight =  winfo.weight * rpa.Picobarn();
    m_ntrial =  winfo.ntrial;
    delete message;
  }
  else {
    m_weight=0.;
  }
  return (m_weight>0.);
}

bool Matrix_Element_Handler::GenerateOneEvent(ATOOLS::Decay_Channel * _dc,double _mass) 
{
  switch (m_mode) {
  case 1: return p_amegic->GetAllDecays()->UnweightedEvent(_dc,_mass);
  }
  return false;
}

double Matrix_Element_Handler::FactorisationScale()
{
  if (m_mode==1) return p_amegic->GetProcess()->FactorisationScale();
  msg.Out()<<" Warning: Matrix_Element_Handler::FactorisationScale() called without AMEGIC!"<<std::endl;
  return 0.;
}

bool Matrix_Element_Handler::GenerateSameEvent() 
{
  if (m_eventmode) {
    bool stat = UnweightedSameEvent();
    if (!stat) return stat;
    GetMomentaNFlavours();
    ApplyHHMF();
    return stat;
  }
  Blob_Data_Base * message = WeightedSameEvent();
  GetMomentaNFlavours();
  ApplyHHMF();
  if (message) {
    PHASIC::Weight_Info winfo = message->Get<PHASIC::Weight_Info>();
    m_weight =  winfo.weight * rpa.Picobarn();
    m_ntrial =  winfo.ntrial;
    delete message;
  }
  else {
    m_weight=0.;
  }
  return (m_weight>0.);
}

void Matrix_Element_Handler::GetMomentaNFlavours() {    
  switch (m_mode) {
  case 1: 
    m_nmoms = p_amegic->Nin()+p_amegic->Nout();
    if (p_amegic->Flavs() && p_amegic->Momenta()) {
      for (int i=0;i<m_nmoms;++i) p_flavs[i] = p_amegic->Flavs()[i];
      for (int i=0;i<m_nmoms;++i) p_moms[i]  = p_amegic->Momenta()[i];
      return;
    }
    break;
  case 2: 
    m_nmoms = p_simplexs->Nin()+p_simplexs->Nout();
    if (p_simplexs->Flavs() && p_simplexs->Momenta()) {
      for (int i=0;i<m_nmoms;++i) p_flavs[i] = p_simplexs->Flavs()[i];
      for (int i=0;i<m_nmoms;++i) p_moms[i]  = p_simplexs->Momenta()[i];
      return;
    }
  }
  msg.Error()<<"Warning in Matrix_Element_Handler::SetMomenta()"<<endl
	     <<"   No ME generator available to get momenta from."<<endl
	     <<"   Continue run and hope for the best."<<endl;
}


bool Matrix_Element_Handler::UnweightedSameEvent() 
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

bool Matrix_Element_Handler::UnweightedEvent() 
{
  switch (m_mode) {
  case 1: return p_amegic->UnweightedEvent();
  case 2: return p_simplexs->UnweightedEvent();
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

ATOOLS::Vec4D * Matrix_Element_Handler::DecMomenta() {
  switch (m_mode) {
  case 1: return p_amegic->GetAllDecays()->Momenta();
  }
  return NULL;
}


ATOOLS::Vec4D * Matrix_Element_Handler::Momenta() {
  return p_moms;
}

ATOOLS::Flavour * Matrix_Element_Handler::Flavs() {
  return p_flavs;
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
  muon.Add(Flavour(kf::e),Flavour(kf::mu));
  muon.Add(Flavour(kf::nue).Bar(),Flavour(kf::numu).Bar());
  m_particle_maps.push_back(muon);

  muon.SetFlipAnti();
  m_particle_maps.push_back(muon);
}

bool  Matrix_Element_Handler::ApplyHHMF()
{
  if (!m_apply_hhmf || m_particle_maps.size()==0) return 1;

  size_t no=(size_t)(m_particle_maps.size()*ran.Get());
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
  msg.Debugging()<<p_simplexs<<std::endl;
  abort();
}

double  Matrix_Element_Handler::Weight() 
{
  return m_weight;
}

unsigned long Matrix_Element_Handler::NumberOfTrials()
{
  return m_ntrial;
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

