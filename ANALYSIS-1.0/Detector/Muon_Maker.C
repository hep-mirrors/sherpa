#include "Muon_Maker.H"
#include "Primitive_Analysis.H"
#include "Detector.H"
#include "Message.H"
#include "MyStrStream.H"
#include "Exception.H"


using namespace ANALYSIS;
using namespace ATOOLS;

DECLARE_GETTER(Muon_Maker_Getter,"Muon_Maker",
	       Analysis_Object,Argument_Matrix);

Analysis_Object *
Muon_Maker_Getter::operator()(const Argument_Matrix &parameters) const
{			
  if (parameters.size()<1) return NULL;
  //if (parameters.size()==1) abort(); // For read-in of, like 'ATLAS'


  std::string mode("ET_UP");
  Muon_Maker * maker = new Muon_Maker(parameters(),mode);
  double Estart,cut_Ebyp,RHad,relHad,cutHad,REM,totEM,cutEM;

  for (size_t i=0;i<parameters.size();++i) {
    const std::vector<std::string> &cur=parameters[i];
    if (cur.size()<2) continue;
    else if (cur[0]=="Track") {
      Estart   = ATOOLS::ToType<double>(cur[1]);
      cut_Ebyp = ATOOLS::ToType<double>(cur[2]);
      maker->SetTracking(Estart);
    }
    else if (cur[0]=="HadIso") {
      RHad   = ATOOLS::ToType<double>(cur[1]);
      relHad = ATOOLS::ToType<double>(cur[2]);
      cutHad = ATOOLS::ToType<double>(cur[3]);
      maker->SetHadIso(RHad,relHad,cutHad);
    }
    else if (cur[0]=="EMIso") {
      REM   = ATOOLS::ToType<double>(cur[1]);
      totEM = ATOOLS::ToType<double>(cur[2]);
      cutEM = ATOOLS::ToType<double>(cur[3]);
      maker->SetEMIso(REM,totEM,cutEM);
    }
    else if (cur[0]=="Ordering") {
      mode = cur[1];
      maker->SetOrdering(mode);
    }
    else if (cur[0]=="ECorrection") {
      maker->SetECorrection(ATOOLS::ToType<double>(cur[1]));
    }
  }
  return maker;
}									

void Muon_Maker_Getter::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"Keyword=Acceptance   Parameters=etamin, etamax"<<std::endl; 
}

Muon_Maker::Muon_Maker(Primitive_Analysis * ana,const std::string mode) :
  Object_Definition_Base(ana,"MuonMaker",mode), 
  m_dim(3),m_Estart(3.),m_R2track(sqr(0.2)),m_cut_Ebyp(2.),
  m_R2hadiso(sqr(0.2)),m_relhad(0.01),m_minhad(1.),
  m_R2EMiso(sqr(0.2)),m_totEM(5.0),m_minEM(1.)
{ 
  m_kfcode=kf::mu;
  GetElements();
}

Muon_Maker::~Muon_Maker() {}

void Muon_Maker::SetTracking(const double E) {
  m_Estart = E; 
}

void Muon_Maker::SetHadIso(const double R,const double rel,const double cut) {
  m_R2hadiso = R*R; m_relhad = rel; m_minhad = cut;
}

void Muon_Maker::SetEMIso(const double R,const double tot,const double cut) {
  m_R2EMiso = R*R; m_totEM = tot; m_minEM = cut;
}

void Muon_Maker::SetECorrection(const double inv) {
  m_inv = inv;
}

void Muon_Maker::ReconstructObjects(Particle_List * plist,ATOOLS::Vec4D & METvector) {
  //std::cout<<METHOD<<std::endl;
  //std::cout<<METHOD<<" for "<<p_ECal->GetHitCells()->size()<<" hit cells in ECal and "
  //	   <<p_HCal->GetHitCells()->size()<<" hit cells in HCal."<<std::endl;

  m_objects.clear();
  GetTracks();
  IsolateTracks();
  CorrectEnergies();

  DropUsedTracks();
  Particle * part;
  while (!m_objects.empty()) {
    part = m_objects.front()->CreateParticle();
    delete m_objects.front();
    m_objects.pop_front();
    plist->push_back(part);
    METvector -= part->Momentum(); 
  }
}

void Muon_Maker::GetTracks() {
  std::list<Track *> tracks;
  p_chambers->GetTracks(tracks);
  if (tracks.size()==0) return;
  std::list<Track *>::iterator trit;
  for (trit=tracks.begin(); trit!=tracks.end(); trit++) {
    if (!(*trit)->used &&
	((*trit)->flav==Flavour(kf::mu)||(*trit)->flav==Flavour(kf::mu).Bar()) &&
	(*trit)->mom[0]>m_Estart) {
      Reconstructed_Object * object = new Reconstructed_Object((*trit));
      object->SetIncludeTracks(true);
      m_objects.push_back(object);
      (*trit)->used = true;
    }
  }
}

void Muon_Maker::IsolateTracks() {
  if (m_objects.size()==0) return;
  double E_HCal, E_ECal;
  std::list<Cell *> * ECal_cells = p_ECal->GetHitCells();
  std::list<Cell *> * HCal_cells = p_HCal->GetHitCells();
  std::list<Cell *>::iterator cit;
  if (!ECal_cells || ECal_cells->size()==0) return;

  double E,eta,phi;
  bool   veto;
  for (ObjectListIterator olit=m_objects.begin();olit!=m_objects.end();) {
    //std::cout<<"--------------------------------------------------"<<std::endl;
    E      = (*olit)->E();
    eta    = (*olit)->Eta();
    phi    = (*olit)->Phi();
    veto   = false;
    E_HCal = E_ECal = 0.;    
    for (cit=HCal_cells->begin();cit!=HCal_cells->end();cit++) {
      if ((*cit)->R2(eta,phi)<m_R2hadiso && !(*olit)->IsIncluded((*cit))) {
	E_HCal += (*cit)->TotalDeposit();
      }
      if (E_HCal>m_minhad && E_HCal/E>m_relhad) {
	veto = true;
	break;
      }
    }
    for (cit=ECal_cells->begin();cit!=ECal_cells->end();cit++) {
      if ((*cit)->R2(eta,phi)<m_R2EMiso && !(*olit)->IsIncluded((*cit))) {
	E_ECal += (*cit)->TotalDeposit();
      }
      if (E_ECal>m_totEM) {
	veto = true;
	break;
      }
    }
    if (veto) {
      (*(*olit)->GetTracks().begin())->used = false;
      olit = m_objects.erase(olit);
    }
    else olit++;     
  }
}

void Muon_Maker::CorrectEnergies() {
  for (ObjectListIterator olit=m_objects.begin();olit!=m_objects.end();olit++) {
    (*olit)->CorrectE(m_inv);
  }  
}
