#include "Electron_Maker.H"
#include "Primitive_Analysis.H"
#include "Detector.H"
#include "Message.H"
#include "MyStrStream.H"
#include "Exception.H"


using namespace ANALYSIS;
using namespace ATOOLS;

DECLARE_GETTER(Electron_Maker_Getter,"Electron_Maker",
	       Analysis_Object,Argument_Matrix);

Analysis_Object *
Electron_Maker_Getter::operator()(const Argument_Matrix &parameters) const
{			
  if (parameters.size()<1) return NULL;
  //if (parameters.size()==1) abort(); // For read-in of, like 'ATLAS'


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
  m_kfcode=kf::e;
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

void Electron_Maker::ReconstructObjects(Particle_List * plist) {
  std::cout<<METHOD<<std::endl;
  m_objects.clear();
  if (!m_elements) GetElements();
  BuildMatchedClusters();
  IsolateClusters();
  CorrectEnergies();

  DropUsedCells();
  Particle * part;
  while (!m_objects.empty()) {
    part = m_objects.front()->CreateParticle();
    delete m_objects.front();
    m_objects.pop_front();
    plist->push_back(part);
  }
  std::cout<<METHOD<<" --> "<<plist->size()<<", therefore "
	   <<p_ECal->GetHitCells()->size()<<"/"<<p_HCal->GetHitCells()->size()<<" total = "
	   <<(p_ECal->GetHitCells()->size()+p_HCal->GetHitCells()->size())<<std::endl;
}

void Electron_Maker::BuildMatchedClusters() {
  //std::cout<<"---------------------------------------------------"<<std::endl
  //	   <<"---------------------------------------------------"<<std::endl;
  std::list<Cell *> * cells(p_ECal->GetHitCells());
  if (!cells || cells->size()==0) return;
  std::cout<<METHOD<<" for "<<cells->size()<<" hit cells in ECal."<<std::endl;
  Cell  * cell(NULL);
  std::list<Track *> * tracks(NULL);
  std::list<Track *>::iterator trit;
  Track * track(NULL);

  double E,eta,phi;
  bool matched(true);
  for (std::list<Cell *>::iterator cit=cells->begin();cit!=cells->end();cit++) {
    cell = (*cit);
    cell->Centroid(eta,phi);
    //std::cout<<METHOD<<" : E = "<<cell->TotalDeposit()
    //	     <<" for "<<cell->ParticleEntries()->begin()->first->Flav()
    //	     <<" at ("<<eta<<", "<<phi<<")."<<std::endl;
    if (cell->TotalDeposit()>m_Estart) {
      //std::cout<<"    ====> seed found in ("<<eta<<","<<phi<<"), "
      //	       <<" cluster with size "<<m_dim<<"."<<std::endl;
      p_cluster = p_ECal->BuildCluster(cell,m_dim,E,eta,phi);
      std::cout<<" built cluster "<<p_cluster<<"("<<p_cluster->size()<<") : "<<E<<std::endl;
      tracks    = p_tracker->GetTracks(eta,phi,m_R2track,m_kfcode); 
      //std::cout<<" out of gettracks : "<<matched<<"/"<<tracks<<"/"
      //	       <<cell->ParticleEntries()->begin()->first->Flav()<<std::endl;
      matched = false;
      std::cout<<tracks<<" "<<tracks->size()<<" "
	       <<tracks->front()<<" "<<tracks->front()->flav<<" "<<m_trackmode<<std::endl;

      if (tracks && tracks->size()>0) {
	if (m_trackmode=="exact") {
	  if (tracks->size()==1 && (*tracks->begin())->flav.Kfcode()==kf::e) {
	    matched = true;
	    track   = (*tracks->begin());
	  }
	}
	else if (m_trackmode=="just_one") {
	  if (tracks->size()==1) {
	    matched = true;
	    track   = tracks->front();
	  }
	}
	else if (m_trackmode=="any_electron") {
	  for (trit=tracks->begin(); trit!=tracks->end(); trit++) {
	    if ((*trit)->flav==Flavour(kf::e)||(*trit)->flav==Flavour(kf::e).Bar()) {
	      matched = true; 
	      track   = (*trit);
	      break;
	    }
	  }
	}
	else matched = true;
      }
      //std::cout<<"   matched = "<<matched<<std::endl;
      if (!matched) delete p_cluster;
      else {
	Reconstructed_Object * object = new Reconstructed_Object((*trit)->flav,E,eta,phi);
	std::cout<<"new electron : E = "<<E<<" at ("<<eta<<", "<<phi<<") :"
		 <<object<<" cluster = "<<p_cluster<<" ("<<p_cluster->size()<<")"<<std::endl;
	object->SetCells(p_cluster);
	object->SetTracks(new std::vector<Track *>);
	object->AddTrack(track);
	m_objects.push_back(object);
	//std::cout<<"Check this: "
	//	 <<(*trit)->flav<<" "<<(*trit)->mom[0]
	//	 <<" "<<(*trit)->mom.Eta()<<" "<<(*trit)->mom.Phi()<<" vs. "
	//	 <<cell->ParticleEntries()->begin()->first->Flav()
	//	 <<"("<<cell->ParticleEntries()->size()<<") vs. "
	//	 <<object->Flav()<<" "<<object->E()
	//	 <<" "<<object->Eta()<<" "<<object->Phi()<<std::endl;
      }
      delete tracks; tracks = NULL; 
    }
  }
}

void Electron_Maker::IsolateClusters() {
  //std::cout<<"--------------------------------------------------"<<std::endl
  //	   <<"--------------------------------------------------"<<std::endl;
  std::cout<<METHOD<<" for "<<m_objects.size()<<std::endl;
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
	//std::cout<<"Distance of "<<(*cit)->Direction()
	//	 <<"("<<(*cit)->Direction().Eta()<<","<<(*cit)->Direction().Phi()<<") = "
	//	 <<(*cit)->R2(eta,phi)<<"("<<m_R2hadiso<<") from "<<eta<<"/"<<phi
	//	 <<" -> E = "<<E_HCal;
      }
      if (E_HCal>m_minhad && E_HCal/E>m_relhad) {
	//std::cout<<" ===> Veto !!!"<<std::endl;
	veto = true;
	break;
      }
      //std::cout<<std::endl;
    }
    for (cit=ECal_cells->begin();cit!=ECal_cells->end();cit++) {
      if ((*cit)->R2(eta,phi)<m_R2EMiso && !(*olit)->IsIncluded((*cit))) {
	E_ECal += (*cit)->TotalDeposit();
	//std::cout<<"Distance of "<<(*cit)->Direction()
	//	 <<"("<<(*cit)->Direction().Eta()<<","<<(*cit)->Direction().Phi()<<") = "
	//	 <<(*cit)->R2(eta,phi)<<"("<<m_R2EMiso<<") from "<<eta<<"/"<<phi
	//	 <<" -> E = "<<E_ECal;
      }
      if (E_ECal>m_totEM) {
	//std::cout<<" ===> Veto !!!"<<std::endl;
	veto = true;
	break;
      }
      //std::cout<<std::endl;
    }
    if (veto) olit = m_objects.erase(olit);
    else olit++;     
  }
  //std::cout<<"--------------------------------------------------"<<std::endl
  //	   <<"--------------------------------------------------"<<std::endl;
}

void Electron_Maker::CorrectEnergies() {
  for (ObjectListIterator olit=m_objects.begin();olit!=m_objects.end();olit++) {
    (*olit)->CorrectE(m_inv);
  }  
}
