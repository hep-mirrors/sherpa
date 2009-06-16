#include "AddOns/Analysis/Detector/Electron_Maker.H"
#include "AddOns/Analysis/Main/Primitive_Analysis.H"
#include "AddOns/Analysis/Detector/Detector.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"


using namespace ANALYSIS;
using namespace ATOOLS;

DECLARE_GETTER(Electron_Maker_Getter,"Electron_Maker",
	       Analysis_Object,Argument_Matrix);

Analysis_Object *
Electron_Maker_Getter::operator()(const Argument_Matrix &parameters) const
{			
  if (parameters.size()<1) return NULL;

  std::string mode("ET_UP");
  Electron_Maker * maker = new Electron_Maker(parameters(),mode);
  double Estart,Rtrack,cut_Ebyp,RHad,relHad,cutHad,REM,totEM,cutEM;
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
      Rtrack   = ATOOLS::ToType<double>(cur[2]);
      cut_Ebyp = ATOOLS::ToType<double>(cur[3]);
      maker->SetTracking(cur[1],Rtrack,cut_Ebyp);
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

void Electron_Maker_Getter::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"Keyword=Acceptance   Parameters=etamin, etamax"<<std::endl; 
}

Electron_Maker::Electron_Maker(Primitive_Analysis * ana,const std::string mode) :
  Object_Definition_Base(ana,"ElectronMaker",mode), 
  m_dim(3),m_Estart(3.),m_R2track(sqr(0.2)),m_cut_Ebyp(2.),
  m_R2hadiso(sqr(0.2)),m_relhad(0.01),m_minhad(1.),
  m_R2EMiso(sqr(0.2)),m_totEM(5.0),m_minEM(1.),
  m_trackmode(std::string("exact"))
{ 
  m_kfcode=kf_e;
  GetElements();
}

Electron_Maker::~Electron_Maker() {}

void Electron_Maker::SetClustering(const double E,const int dim) {
  m_Estart = E; m_dim = dim;
}

void Electron_Maker::SetTracking(const std::string mode,const double R,const double cut) {
  m_trackmode = mode; m_R2track = R*R; m_cut_Ebyp = cut;
}

void Electron_Maker::SetHadIso(const double R,const double rel,const double cut) {
  m_R2hadiso = R*R; m_relhad = rel; m_minhad = cut;
}

void Electron_Maker::SetEMIso(const double R,const double tot,const double cut) {
  m_R2EMiso = R*R; m_totEM = tot; m_minEM = cut;
}

void Electron_Maker::SetECorrection(const double inv) {
  m_inv = inv;
}

void Electron_Maker::ReconstructObjects(Particle_List * plist,ATOOLS::Vec4D & METvector) {
  m_objects.clear();

  BuildMatchedClusters();
  IsolateClusters();
  CorrectEnergies();

  Particle * part;
  msg_Debugging()<<METHOD<<":"<<std::endl;
  while (!m_objects.empty()) {
    part = m_objects.front()->CreateParticle();
    msg_Debugging()<<"   Found electron : "<<m_objects.front()->Mom()
    	     <<"/"<<part->Momentum()<<" with "<<m_objects.front()->GetCells().size()
    	     <<"/"<<m_objects.front()->GetTracks().size()<<std::endl;
    delete m_objects.front();
    m_objects.pop_front();
    plist->push_back(part);
    METvector -= part->Momentum(); 
  }
  DropUsedCells();
}

void Electron_Maker::BuildMatchedClusters() {
  std::list<Cell *> * cells(p_ECal->GetHitCells());
  if (!cells || cells->size()==0) return;
  Cell  * cell(NULL);
  std::list<Track *> tracks;
  std::list<Track *>::iterator trit;
  std::vector<Cell *> cluster;
  Track * track(NULL);

  msg_Debugging()<<METHOD<<":"<<std::endl;
  double E,eta,phi;
  bool matched(true);
  for (std::list<Cell *>::iterator cit=cells->begin();cit!=cells->end();cit++) {
    cell = (*cit);
    cell->Centroid(eta,phi);
    if (!cell->Used() && cell->TotalDeposit()>m_Estart) {
      msg_Debugging()<<"   Found potential e-seed: "
	       <<cell->ParticleEntries()->begin()->first->Flav()<<" with "
	       <<cell->ParticleEntries()->begin()->first->Momentum()<<" --> ";

      cluster.clear();
      tracks.clear();
      p_ECal->BuildCluster(cell,cluster,m_dim,E,eta,phi);
      p_tracker->GetTracks(tracks,eta,phi,m_R2track,m_kfcode); 
      matched = false;
 
      if (tracks.size()>0) {
	if (m_trackmode=="exact") {
	  if (tracks.size()==1 && !(*tracks.begin())->used &&
	      (*tracks.begin())->flav.Kfcode()==kf_e) {
	    matched = true;
	    track   = (*tracks.begin());
	  }
	}
	else if (m_trackmode=="just_one") {
	  msg_Debugging()<<"(too many tracks : "<<tracks.size()<<", "
		   <<tracks.front()->flav<<")";
	  if (tracks.size()==1 && !tracks.front()->used) {
	    matched = true;
	    track   = tracks.front();
	  }
	}
	else if (m_trackmode=="any_electron") {
	  for (trit=tracks.begin(); trit!=tracks.end(); trit++) {
	    if (!(*trit)->used &&
		((*trit)->flav==Flavour(kf_e)||(*trit)->flav==Flavour(kf_e).Bar())) {
	      matched = true; 
	      track   = (*trit);
	      break;
	    }
	  }
	}
	else matched = true;
      }
      if (matched) {
	msg_Debugging()<<" matched, take it."<<std::endl;
	Flavour flav = (track->flav.Charge()>0)?Flavour(m_kfcode).Bar():Flavour(m_kfcode);
	Reconstructed_Object * object(new Reconstructed_Object(flav,E,eta,phi));
	object->SetCells(cluster);
	object->AddTrack(track);
	track->used = true;
	m_objects.push_back(object);
      }
      else {
	msg_Debugging()<<" not matched, discard it."<<std::endl;
	for (std::vector<Cell *>::iterator cur=cluster.begin();cur!=cluster.end();cur++)
	  (*cur)->SetUsed(false);
      }
    }
  }
  msg_Debugging()<<METHOD<<": "<<m_objects.size()
	   <<" track-matched candidates waiting for isolation."<<std::endl;
}

void Electron_Maker::IsolateClusters() {
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

void Electron_Maker::CorrectEnergies() {
  for (ObjectListIterator olit=m_objects.begin();olit!=m_objects.end();olit++) {
    (*olit)->CorrectE(1./m_inv);
  }  
}
