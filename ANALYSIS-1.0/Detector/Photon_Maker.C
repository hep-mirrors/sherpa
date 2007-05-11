#include "Photon_Maker.H"
#include "Primitive_Analysis.H"
#include "Detector.H"
#include "Message.H"
#include "MyStrStream.H"
#include "Exception.H"


using namespace ANALYSIS;
using namespace ATOOLS;

DECLARE_GETTER(Photon_Maker_Getter,"Photon_Maker",
	       Analysis_Object,Argument_Matrix);

Analysis_Object *
Photon_Maker_Getter::operator()(const Argument_Matrix &parameters) const
{			
  if (parameters.size()<1) return NULL;
  //if (parameters.size()==1) abort(); // For read-in of, like 'ATLAS'


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
  m_kfcode=kf::photon;
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

void Photon_Maker::ReconstructObjects(Particle_List * plist) {
  //std::cout<<METHOD<<std::endl;
  m_objects.clear();
  BuildMatchedClusters();
  IsolateClusters();
  CorrectEnergies();
  DropUsedCells();

  Particle * part;
  //std::cout<<METHOD<<" : "<<m_objects.size()<<std::endl;
  while (!m_objects.empty()) {
    part = m_objects.front()->CreateParticle();
    plist->push_back(part);
    delete m_objects.front();
    m_objects.pop_front();
  }
  //std::cout<<METHOD<<" --> "<<plist->size()<<", therefore "
  //	   <<p_ECal->GetHitCells()->size()<<"/"<<p_HCal->GetHitCells()->size()<<" total = "
  //	   <<(p_ECal->GetHitCells()->size()+p_HCal->GetHitCells()->size())<<std::endl;
}

void Photon_Maker::BuildMatchedClusters() {
  //std::cout<<"---------------------------------------------------"<<std::endl
  //	   <<"---------------------------------------------------"<<std::endl
  //	   <<METHOD<<" : check ecal = "<<p_ECal<<"."<<std::endl;
  std::list<Cell *> * cells = p_ECal->GetHitCells();
  if (!cells || cells->size()==0) return;
  //std::cout<<METHOD<<" for "<<cells->size()<<" hit cells in ECal."<<std::endl;
  Cell * cell;
  std::list<Track *> * tracks;
  std::list<Track *>::iterator trit;
  double E,eta,phi;
  bool vetoit;
  for (std::list<Cell *>::iterator cit=cells->begin();cit!=cells->end();cit++) {
    cell = (*cit);
    //std::cout<<METHOD<<" : E = "<<cell->TotalDeposit()
    //	     <<" for "<<cell->ParticleEntries()->begin()->first->Flav()<<":"<<std::endl;
    if (cell->TotalDeposit()>m_Estart) {
      cell->Centroid(eta,phi);
      //std::cout<<"    ====> seed found in ("<<eta<<","<<phi<<"), "
      //	       <<" cluster with size "<<m_dim<<"."<<std::endl;
      p_cluster = p_ECal->BuildCluster(cell,m_dim,E,eta,phi);
      //std::cout<<"Built cluster "<<p_cluster<<"("<<p_cluster->size()<<") : "<<E<<std::endl;
      tracks    = p_tracker->GetTracks(eta,phi,m_R2track,kf::none);
      //std::cout<<"   Tracker "<<p_tracker<<" gives "<<tracks<<" "<<tracks->size()<<" "
      //	       <<tracks->front()<<"."<<std::endl;
      if (tracks && tracks->size()>0) {
	vetoit = false;
	for (trit=tracks->begin(); trit!=tracks->end(); trit++) {
	  if ((*trit)->mom.PPerp()>m_trackcut) {
	    vetoit=true; break;
	  }
	} 
      }
      //std::cout<<"   vetoit = "<<vetoit<<std::endl;
      if (vetoit) { 
	//std::cout<<"   ... delete testcluster : "<<p_cluster<<std::endl;
	delete p_cluster; p_cluster = NULL;
      } 
      else {
	Reconstructed_Object * object = new Reconstructed_Object(m_kfcode,E,eta,phi);
	object->SetCells(p_cluster);
	m_objects.push_back(object);
	//std::cout<<"   new photon ("<<m_objects.size()<<") : "
	//	 <<"E = "<<E<<" at ("<<eta<<", "<<phi<<") : object = "
	//	 <<object<<" cluster = "<<p_cluster<<" ("<<p_cluster->size()<<")"<<std::endl;
	//std::cout<<"Check this: "
	//	       <<cell->ParticleEntries()->begin()->first->Flav()
	//	       <<"("<<cell->ParticleEntries()->size()<<") vs. "
	//	       <<object->Flav()<<" "<<object->E()
	//	       <<" "<<object->Eta()<<" "<<object->Phi()<<std::endl;
      }
      //std::cout<<"   ... delete testtracks : "<<tracks<<std::endl;      
      delete tracks; tracks = NULL; 
    }
  }
  //std::cout<<" ............ out of "<<METHOD<<std::endl;
}

void Photon_Maker::IsolateClusters() {
  //std::cout<<METHOD<<std::endl;
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
	//std::cout<<"Distance of "<<(*cit)->Direction()
	//	 <<"("<<(*cit)->Direction().Eta()<<","<<(*cit)->Direction().Phi()<<") = "
	//	 <<(*cit)->R2(eta,phi)<<"("<<m_R2hadiso<<") from "<<eta<<"/"<<phi
	//	 <<" -> E = "<<E_HCal<<" ===> Veto !!!"<<std::endl;
	veto = true;
	break;
      }
    }
    for (cit=ECal_cells->begin();cit!=ECal_cells->end();cit++) {
      if ((*cit)->R2(eta,phi)<m_R2EMiso && !(*olit)->IsIncluded((*cit))) {
	E_ECal += (*cit)->TotalDeposit();
      }
      if (E_ECal>m_totEM) {
	//std::cout<<"Distance of "<<(*cit)->Direction()
	//	 <<"("<<(*cit)->Direction().Eta()<<","<<(*cit)->Direction().Phi()<<") = "
	//	 <<(*cit)->R2(eta,phi)<<"("<<m_R2EMiso<<") from "<<eta<<"/"<<phi
	//	 <<" -> E = "<<E_ECal<<" ===> Veto !!!"<<std::endl;
	veto = true;
	break;
      }
    }
    if (veto) {
      delete (*olit);
      olit = m_objects.erase(olit);
    }
    else olit++;     
  }
  //std::cout<<"--------------------------------------------------"<<std::endl
  //	   <<"--------------------------------------------------"<<std::endl;
}

void Photon_Maker::CorrectEnergies() {
  for (ObjectListIterator olit=m_objects.begin();olit!=m_objects.end();olit++) {
    (*olit)->CorrectE(m_inv);
  }  
}
