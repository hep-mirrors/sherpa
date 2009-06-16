#include "AddOns/Analysis/Detector/Photon_Maker.H"
#include "AddOns/Analysis/Main/Primitive_Analysis.H"
#include "AddOns/Analysis/Detector/Detector.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"


using namespace ANALYSIS;
using namespace ATOOLS;

DECLARE_GETTER(Photon_Maker_Getter,"Photon_Maker",
	       Analysis_Object,Argument_Matrix);

Analysis_Object *
Photon_Maker_Getter::operator()(const Argument_Matrix &parameters) const
{			
  if (parameters.size()<1) return NULL;

  std::string mode("ET_UP");
  Photon_Maker * maker = new Photon_Maker(parameters(),mode);
  double Estart,Rtrack,trackcut,RHad,relHad,cutHad,REM,totEM,cutEM;
  int    size;

  for (size_t i=0;i<parameters.size();++i) {
    const std::vector<std::string> &cur=parameters[i];
    if (cur.size()<2) continue;
    else if (cur[0]=="Cluster") {
      Estart = ATOOLS::ToType<double>(cur[1]);
      size   = ATOOLS::ToType<int>(cur[2]);
      maker->SetClustering(Estart,size);
    }
    else if (cur[0]=="Tracking") {
      Rtrack   = ATOOLS::ToType<double>(cur[1]);
      trackcut = ATOOLS::ToType<double>(cur[2]);
      maker->SetTracking(Rtrack,trackcut);
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

void Photon_Maker_Getter::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"Keyword=Acceptance   Parameters=etamin, etamax"<<std::endl; 
}

Photon_Maker::Photon_Maker(Primitive_Analysis * ana,const std::string mode) :
  Object_Definition_Base(ana,"PhotonMaker",mode), 
  m_dim(3),m_Estart(3.),m_R2track(sqr(0.2)),m_trackcut(2.),
  m_R2hadiso(sqr(0.2)),m_relhad(0.01),m_minhad(1.),
  m_R2EMiso(sqr(0.2)),m_totEM(5.0),m_minEM(1.)
{ 
  m_kfcode=kf_photon;
  GetElements();
}

Photon_Maker::~Photon_Maker() {}

void Photon_Maker::SetClustering(const double E,const int dim) {
  m_Estart = E; m_dim = dim;
}

void Photon_Maker::SetTracking(const double R,const double cut) {
  m_R2track = R*R; m_trackcut = cut;
}

void Photon_Maker::SetHadIso(const double R,const double rel,const double cut) {
  m_R2hadiso = R*R; m_relhad = rel; m_minhad = cut;
}

void Photon_Maker::SetEMIso(const double R,const double tot,const double cut) {
  m_R2EMiso = R*R; m_totEM = tot; m_minEM = cut;
}

void Photon_Maker::SetECorrection(const double inv) {
  m_inv = inv;
}

void Photon_Maker::ReconstructObjects(Particle_List * plist,ATOOLS::Vec4D & METvector) {
  m_objects.clear();
  BuildMatchedClusters();
  IsolateClusters();
  CorrectEnergies();

  msg_Debugging()<<METHOD<<":"<<std::endl;
  Particle * part;
  while (!m_objects.empty()) {
    part = m_objects.front()->CreateParticle();
    plist->push_back(part);
    msg_Debugging()<<"   Found photon : "<<m_objects.front()->Mom()
    	     <<"/"<<part->Momentum()<<" with "<<m_objects.front()->GetCells().size()
    	     <<"/"<<m_objects.front()->GetTracks().size()<<std::endl;
    delete m_objects.front();
    m_objects.pop_front();
    METvector -= part->Momentum(); 
  }
  DropUsedCells();
}

void Photon_Maker::BuildMatchedClusters() {
  msg_Debugging()<<METHOD<<":"<<std::endl;
  std::list<Cell *> * cells = p_ECal->GetHitCells();
  if (!cells || cells->size()==0) return;
  Cell * cell;
  std::list<Track *> tracks;
  std::list<Track *>::iterator trit;
  std::vector<Cell *> cluster;
  double E,eta,phi;
  bool vetoit(false);
  for (std::list<Cell *>::iterator cit=cells->begin();cit!=cells->end();cit++) {
    cell = (*cit);
    if (!cell->Used() && cell->TotalDeposit()>m_Estart) {
      msg_Debugging()<<"   Found potential gamma-seed: "
	       <<cell->ParticleEntries()->begin()->first->Flav()<<" with "
	       <<cell->ParticleEntries()->begin()->first->Momentum()<<" --> ";

      cell->Centroid(eta,phi);
      cluster.clear();
      tracks.clear();
      p_ECal->BuildCluster(cell,cluster,m_dim,E,eta,phi);
      p_tracker->GetTracks(tracks,eta,phi,m_R2track,kf_none);
      if (tracks.size()>0) {
	vetoit = false;
	for (trit=tracks.begin(); trit!=tracks.end(); trit++) {
	  if ((*trit)->mom.PPerp()>m_trackcut) { vetoit=true; break; }
	} 
      }
      if (!vetoit) {
	msg_Debugging()<<" not track-matched, take it."<<std::endl;
	Reconstructed_Object * object = new Reconstructed_Object(m_kfcode,E,eta,phi);
	object->SetCells(cluster);
	m_objects.push_back(object);
      }
      else {
	msg_Debugging()<<" but track-matched, veto it."<<std::endl;
	for (std::vector<Cell *>::iterator cur=cluster.begin();cur!=cluster.end();cur++)
	  (*cur)->SetUsed(false);
      }
    }
  }
}

void Photon_Maker::IsolateClusters() {
  msg_Debugging()<<METHOD<<":"<<std::endl;
  if (m_objects.size()==0) return;
  double E_HCal, E_ECal;
  std::list<Cell *> * ECal_cells = p_ECal->GetHitCells();
  std::list<Cell *> * HCal_cells = p_HCal->GetHitCells();
  std::list<Cell *>::iterator cit;
  if (!ECal_cells || ECal_cells->size()==0) return;

  double E,eta,phi;
  bool   veto;
  for (ObjectListIterator olit=m_objects.begin();olit!=m_objects.end();) {
    E      = (*olit)->E();
    eta    = (*olit)->Eta();
    phi    = (*olit)->Phi();
    veto   = false;
    E_HCal = E_ECal = 0.;    
    for (cit=HCal_cells->begin();cit!=HCal_cells->end();cit++) {
      if ((*cit)->R2(eta,phi)<m_R2hadiso && !(*olit)->IsIncluded((*cit))) {
	E_HCal += (*cit)->TotalDeposit();
      }
      if (E_HCal>m_minhad && E_HCal/E>m_relhad) { veto = true; break; }
    }
    for (cit=ECal_cells->begin();cit!=ECal_cells->end();cit++) {
      if ((*cit)->R2(eta,phi)<m_R2EMiso && !(*olit)->IsIncluded((*cit))) {
	E_ECal += (*cit)->TotalDeposit();
      }
      if (E_ECal>m_totEM) { veto = true; break; }
    }
    msg_Debugging()<<"    Check for isolation of "<<(*olit)->Flav()<<" --> ";
    if (veto) { 
      msg_Debugging()<<" not isolated (HCal = "<<E_HCal<<", ECAL = "<<E_ECal<<"), veto it."<<std::endl;
      delete (*olit); 
      olit = m_objects.erase(olit); 
    }
    else { 
      msg_Debugging()<<" keep it."<<std::endl;
      olit++;
    }
  }
}

void Photon_Maker::CorrectEnergies() {
  for (ObjectListIterator olit=m_objects.begin();olit!=m_objects.end();olit++) {
    (*olit)->CorrectE(1./m_inv);
  }  
}
