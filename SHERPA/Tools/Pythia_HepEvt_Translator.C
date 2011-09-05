#include "SHERPA/Tools/Pythia_HepEvt_Translator.H"
#include "SHERPA/Tools/HepEvt_Interface.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

Pythia_HepEvt_Translator::Pythia_HepEvt_Translator(HepEvt_Interface * interface) :
  p_interface(interface), p_blobs(NULL)
{}

bool Pythia_HepEvt_Translator::ConstructBlobs(ATOOLS::Blob_List * const blobs)
{
  p_blobs = blobs;

  m_bunchints.clear();
  m_signalints.clear();
  CopyHepEvtBlock();

  HepEvt2Particles();
  ReconstructBeamsAndBunches();
  if (ReconstructSignalBlob()==false)        return false;
  if (ReconstructShowerBlob()==false)        return false;
  if (ReconstructFragmentationBlob()==false) return false;
  CleanUp();
  //cout<<(*p_blobs)<<endl<<endl;
  return true;
}

void Pythia_HepEvt_Translator::CopyHepEvtBlock() {
  m_nhep   = p_interface->m_nhep;
  p_isthep = p_interface->p_isthep;
  p_idhep  = p_interface->p_idhep;
  p_jmohep = p_interface->p_jmohep;
  p_jdahep = p_interface->p_jdahep;
  p_phep   = p_interface->p_phep;
  p_vhep   = p_interface->p_vhep;
}

void Pythia_HepEvt_Translator::HepEvt2Particles()
{
  if (!m_convertH2S.empty()) {
    for (Translation_Map::iterator m_piter=m_convertH2S.begin();
	 m_piter!=m_convertH2S.end();m_piter++) {
      if (m_piter->second.second) {
	delete (m_piter->second.first); m_piter->second.first=NULL; 
      }
    }
    m_convertH2S.clear();
  }

  Flavour flav;
  for (int pos=0;pos<m_nhep;pos++) {
    if (abs(p_idhep[pos])==9902210) return;
    //cout<<pos<<": stat,id "<<p_isthep[pos]<<","<<p_idhep[pos]
    //	     <<"; mos : "<<p_jmohep[2*pos]<<","<<p_jmohep[2*pos+1]
    //	     <<"; das : "<<p_jdahep[2*pos]<<","<<p_jdahep[2*pos+1]<<endl
    //	     <<"    mom: "<<p_phep[5*pos+3]<<" "<<p_phep[5*pos+0]<<" "<<p_phep[5*pos+1]
    //	     <<" "<<p_phep[5*pos+2]<<" "<<p_phep[5*pos+4]<<endl;

    flav.FromHepEvt(p_idhep[pos]);
    Vec4D momentum     = Vec4D(p_phep[3+pos*5],p_phep[0+pos*5],p_phep[1+pos*5],p_phep[2+pos*5]);
    Particle * newpart = new Particle(pos+1,flav,momentum);
    newpart->SetFinalMass(p_phep[4+pos*5]);
    newpart->SetStatus(part_status::code(p_isthep[pos]));
    m_convertH2S[pos]=pair<Particle*,bool>(newpart,true);
    //cout<<"   "<<(*newpart)<<endl;
  }
}

bool Pythia_HepEvt_Translator::ReconstructBeamsAndBunches()
{
  int        pos;
  Particle * part, * mother;
  Blob     * blob;
  for (m_piter=m_convertH2S.begin();m_piter!=m_convertH2S.end();m_piter++) {
    if (!m_piter->second.second) continue;
    pos  = m_piter->first;
    part = m_piter->second.first;
    if (part->Status()==part_status::documentation) {
      if ((pos==0 && part->Flav()==rpa->gen.Beam1()) ||
	  (pos==1 && part->Flav()==rpa->gen.Beam2())) {
	blob = new Blob();
	blob->SetId();
	blob->SetType(btp::Beam);
	blob->SetStatus(blob_status::inactive);
	blob->AddToOutParticles(part);
	part->SetStatus(part_status::active);
	p_blobs->push_back(blob);
	m_piter->second.second = false;
      }
    }
    if (pos>1) {
      if ((p_jmohep[2*pos]-1==0 && p_jmohep[2*pos+1]-1==-1) ||
	  (p_jmohep[2*pos]-1==1 && p_jmohep[2*pos+1]-1==-1)) {
	m_piter1 = m_convertH2S.find(p_jmohep[2*pos]-1); 
	if (m_piter1==m_convertH2S.end()) {
	  msg_Error()<<"WARNING : Error in "<<METHOD<<" : "<<endl
		     <<"    Potential bunch particle does not have beam particle."<<endl
		     <<"    Will return .false. and hope that event is discarded."<<endl;
	  return false;
	}
	else mother = m_piter1->second.first;
	blob = mother->DecayBlob(); 
	if (blob==NULL) {
	  blob = new Blob();
	  blob->SetId();
	  blob->SetType(btp::Bunch);
	  blob->SetStatus(blob_status::inactive);
	  blob->AddToInParticles(mother);
	  part->SetStatus(part_status::decayed);
	  p_blobs->push_back(blob);
	}
	blob->AddToOutParticles(part);
	part->SetStatus(part_status::active);
	if (part->Flav().IsHadron() && !part->Flav().IsDiQuark()) {
	  int dau1 = p_jdahep[2*pos]-1, dau2 = p_jdahep[2*pos+1]-1;
	  if (dau1>pos) {
	    Blob * blob = new Blob();
	    blob->SetId();
	    blob->SetStatus(blob_status::inactive);
	    blob->SetType(btp::Hadron_Decay);
	    blob->AddToInParticles(part);
	    p_blobs->push_back(blob);
	    for (int i=dau1;i<=dau2;i++) ReconstructDecayChain(blob,i);
	  }
	}   
	else m_bunchints.insert(pos);
	m_piter->second.second = false;
      }
    }
  }
  return true;
}

bool Pythia_HepEvt_Translator::ReconstructSignalBlob()
{
  int        pos;
  Particle * part, * mother;

  p_signal=NULL; 

  for (m_piter=m_convertH2S.begin();m_piter!=m_convertH2S.end();m_piter++) {
    if (!m_piter->second.second) continue;
    pos  = m_piter->first;
    part = m_piter->second.first;
    //cout<<"   Particle : "<<pos<<" : "<<part->Flav()<<" "<<part->Status()<<endl;
    if (part->Status()==part_status::documentation) {
      //cout<<"   Belongs to signal."<<endl;
      if (IsZero(part->Momentum()[0])) {
	msg_Error()<<"WARNING : Error in "<<METHOD<<" : "<<endl
		   <<"    Signal particles with zero energy: Looks like a nonsensical event."<<endl
		   <<"    Will return .false. and hope that event is discarded."<<endl;
	return false;
      }
      m_signalints.insert(pos);
    }
  }
  if (m_signalints.size()<=0) {
    msg_Error()<<"WARNING : Error in "<<METHOD<<" : "<<endl
	       <<"    Zero signal particles: Looks like a nonsensical event."<<endl
	       <<"    Will return .false. and hope that event is discarded."<<endl;
    return false;
  }

  p_signal = new Blob;
  p_signal->SetId();
  p_signal->SetStatus(blob_status::inactive);
  p_signal->ClearAllData();
  p_signal->AddData("ME_Weight",new Blob_Data<double>(1.));
  p_signal->SetType(btp::Signal_Process);
  p_blobs->push_back(p_signal);

  for (m_spiter=m_signalints.begin();m_spiter!=m_signalints.end();m_spiter++) {
    m_piter = m_convertH2S.find(*m_spiter);
    pos   = p_jmohep[2*(*m_spiter)]-1;
    part  = m_piter->second.first;
    if (pos==-1) {
      p_signal->AddToOutParticles(part);
      part->SetStatus(part_status::active);
      m_piter->second.second = false;
    }
    else {
      m_piter1 = m_convertH2S.find(pos); 
      mother = m_piter1->second.first;
      if (mother->ProductionBlob() && 
	  mother->ProductionBlob()->Type()==btp::Bunch) {
	p_signal->AddToInParticles(part);
	part->SetStatus(part_status::decayed);
	m_piter->second.second = false;
      }
      else if ((mother->DecayBlob() && 
		mother->DecayBlob()->Type()==btp::Signal_Process) ||
	       (mother->ProductionBlob() && 
		mother->ProductionBlob()->Type()==btp::Signal_Process)) {
	p_signal->AddToOutParticles(part);
	part->SetStatus(part_status::active);
	m_piter->second.second = false;
      }
    }
  }
  for (m_spiter=m_signalints.begin();m_spiter!=m_signalints.end();m_spiter++) {
    pos    = p_jmohep[2*(*m_spiter)]-1;
    m_piter  = m_convertH2S.find(pos); 
    mother = m_piter->second.first;
    if (mother->ProductionBlob() && 
	mother->ProductionBlob()->Type()==btp::Signal_Process) {
      mother = p_signal->RemoveOutParticle(mother,true);
      m_piter->second.second = true;
    }
  }

  //cout<<(*p_signal)<<endl;
  return true;
}

bool Pythia_HepEvt_Translator::ReconstructShowerBlob() {
  set<int>   partints;
  for (m_piter=m_convertH2S.begin();m_piter!=m_convertH2S.end();m_piter++) {
    if (!m_piter->second.second) continue;
    if (m_piter->second.first->Status()==part_status::documentation) {
      msg_Error()<<"WARNING : Error in "<<METHOD<<" : "<<endl
		 <<"    Shower/hadronization particle with documentation tag?"<<endl
		 <<"    "<<(*m_piter->second.first)<<endl
		 <<"    Will return .false. and hope that event is discarded."<<endl;
      return false;
    }
    if (!m_piter->second.first->Flav().IsHadron()) partints.insert(m_piter->first);
  }


  if (partints.size()<=0) {
    msg_Error()<<"WARNING : Error in "<<METHOD<<" : "<<endl
	       <<"    Zero shower particles: Will return .true. and hope for the best."<<endl;
    return true;
  }

  int        pos;
  Particle * part, * mother;

  Blob * blob = new Blob;
  blob->SetId();
  blob->SetStatus(blob_status::inactive);
  blob->SetType(btp::Shower);
  p_blobs->push_back(blob);

  // Naively: 
  // Add signals outparticles as inparticles of shower (FS initialisers), and
  // add signals inparticles as outparticles of shower (IS products).
  for (int i=0;i<p_signal->NInP();i++)  {
    part = p_signal->InParticle(i);
    part->SetStatus(part_status::decayed);
    blob->AddToOutParticles(part);
  }
  for (int i=0;i<p_signal->NOutP();i++) {
    part = p_signal->OutParticle(i);
    part->SetStatus(part_status::decayed);
    blob->AddToInParticles(part);
  }
  for (m_spiter=m_bunchints.begin();m_spiter!=m_bunchints.end();m_spiter++) {
    m_piter = m_convertH2S.find((*m_spiter));
    part    = m_piter->second.first;
    part->SetStatus(part_status::decayed);
    blob->AddToInParticles(part);
  }

  bool flag;
  set<int>::iterator sit;
  int dau1,dau2;
  
  for (m_spiter=partints.begin();m_spiter!=partints.end();m_spiter++) {
    m_piter = m_convertH2S.find(*m_spiter);
    if (!m_piter->second.second) continue;
    part    = m_piter->second.first;
    if (!part->Flav().IsHadron() && 
	part->Flav()!=Flavour(kf_cluster) &&
	part->Flav()!=Flavour(kf_string)) {
      pos   = p_jmohep[2*(*m_spiter)]-1;
      flag  = false;
      //cout<<"Check for "<<(*m_spiter)<<", mother = "<<pos<<" ("
      //	  <<(partints.find(pos)!=partints.end())<<") "
      //	  <<part->Flav()<<" "<<m_piter->second.second<<endl;

      // Not so obvious:  Check if mother is a shower particle - 
      // happens if, e.g., low-mass clusters neccessitate momenta
      // shuffling of final state particles in shower.
      if ((*m_spiter)>pos && partints.find(pos)!=partints.end()) {
	//cout<<"Found a match (FS->FS): "<<pos<<" -> "<<(*m_spiter)
	//   <<" : "<<part->Flav()<<endl;
	m_piter1 = m_convertH2S.find(pos);
	blob->RemoveOutParticle(m_piter1->second.first);
	m_piter1->second.second = true;
	flag = true;
      }
      // Obvious Checks, for "normal" shower:
      // FS Shower:  if mother in signal, but particle not in signal   .or.     
      // IS Shower:  if same mother as a signal particle, which itself is not a signal particle .or.
      //             if mother in bunch.
      if (!flag) {
	for (set<int>::iterator auntit=m_signalints.begin();
	     auntit!=m_signalints.end();auntit++) {
	  if ((pos==(*auntit) && m_signalints.find((*m_spiter))==m_signalints.end()) ||
	      (p_jmohep[2*(*auntit)]-1==pos &&
	       m_signalints.find(p_jmohep[2*(*auntit)]-1)==m_signalints.end()) ||
	      (m_bunchints.find(pos)!=m_bunchints.end())) {
	    //cout<<"Found a match (IS/FS): "<<pos<<" -> "<<(*m_spiter)
	    //	<<" : "<<part->Flav()<<endl;
	    flag = true;
	    break;
	  }
	}
      }
      // Backup Checks:
      // If not identifed as shower particle so far:
      // if mother has signal as production or decay blob
      if (!flag) {
	m_piter1 = m_convertH2S.find(pos);
	if (m_piter1!=m_convertH2S.end()) mother = m_piter1->second.first;
	                             else mother = NULL;
	if (m_signalints.find(*m_spiter)==m_signalints.end() &&
	    (mother && (mother->DecayBlob()==p_signal || 
			mother->ProductionBlob()==p_signal))) {
	  //cout<<"Found a match (Alternative): "<<pos<<" -> "<<(*m_spiter)
	  //   <<" : "<<part->Flav()<<endl;
	  flag = true;
	}
      }
      if (flag) {
	blob->AddToOutParticles(part);
	part->SetStatus(part_status::active);
	m_piter->second.second = false;
	dau1 = p_jdahep[2*(*m_spiter)]-1;
	dau2 = p_jdahep[2*(*m_spiter)+1]-1;
	m_piter1 = m_convertH2S.find(dau1);
        // If decay chain
	if (dau1!=dau2 && 
	    ((partints.find(dau1)==partints.end() || partints.find(dau2)==partints.end()) ||
	     part->Flav().Kfcode()==kf_tau)) {
	  //cout<<"Decay this : "<<(*m_spiter)<<" -> "
	  //   <<dau1<<" / "<<dau2<<" "<<m_piter1->second.first->Flav()<<"."<<endl;
	  Blob * blob = new Blob();
	  blob->SetId();
	  blob->SetStatus(blob_status::inactive);
	  blob->SetType(btp::Hadron_Decay);
	  blob->AddToInParticles(part);
	  p_blobs->push_back(blob);
	  part->SetStatus(part_status::fragmented);
	  for (int i=dau1;i<=dau2;i++) ReconstructDecayChain(blob,i);
	}
	else if (dau1!=dau2 && 
		 partints.find(dau1)!=partints.end() && partints.find(dau2)!=partints.end()) {
	  blob->RemoveOutParticle(part);
	  m_piter->second.second = true;
	}
      }
    }
  }
  return true;
}

bool Pythia_HepEvt_Translator::ReconstructFragmentationBlob() {
  Blob * blob = new Blob();
  blob->SetId();
  blob->SetStatus(blob_status::inactive);
  blob->SetType(btp::Fragmentation);
  p_blobs->push_back(blob);

  Particle * part,* extra;
  int pos,frag,dau1,dau2;
  for (m_piter=m_convertH2S.begin();m_piter!=m_convertH2S.end();m_piter++) {
    pos      = m_piter->first;
    part     = m_piter->second.first;
    frag     = p_jdahep[2*pos]-1;
    m_piter1 = m_convertH2S.find(frag);
    if (m_piter1->second.first->Flav()==Flavour(kf_string) ||
	m_piter1->second.first->Flav()==Flavour(kf_cluster)) {
      // Make sure, particle does not show up in 2 blobs - may happen to bunchparticles.
      if (part->DecayBlob()) {
	//cout<<METHOD<<" EXTRA!!!"<<endl;
	extra = new Particle((*part));
	extra->SetNumber(0);
	extra->SetFinalMass(part->FinalMass());
	part->DecayBlob()->AddToOutParticles(extra);
	blob->AddToInParticles(extra);
      }
      else blob->AddToInParticles(part);
      m_piter->second.first->SetStatus(part_status::fragmented);
      dau1 = p_jdahep[2*frag]-1;
      dau2 = p_jdahep[2*frag+1]-1;
      for (int i=dau1;i<=dau2;i++) ReconstructDecayChain(blob,i);
    }
  }
  return true;
}

void Pythia_HepEvt_Translator::ReconstructDecayChain(ATOOLS::Blob * prod,const int pos)
{
  ATOOLS::Translation_Map::iterator pit = m_convertH2S.find(pos);
  Particle * part = pit->second.first;
  if (!pit->second.second || part->DecayBlob()!=NULL) return;
  prod->AddToOutParticles(part);
  prod->SetPosition(Vec4D(p_vhep[3+pos*4],p_vhep[0+pos*4],p_vhep[1+pos*4],p_vhep[2+pos*4]));
  pit->second.second = false;
  int dau1 = p_jdahep[2*pos]-1, dau2 = p_jdahep[2*pos+1]-1;
  if (dau1>pos) {
    Blob * blob = new Blob();
    blob->SetId();
    blob->SetStatus(blob_status::inactive);
    blob->SetType(btp::Hadron_Decay);
    blob->AddToInParticles(part);
    p_blobs->push_back(blob);
    for (int i=dau1;i<=dau2;i++) ReconstructDecayChain(blob,i);
  }
}

void Pythia_HepEvt_Translator::CleanUp() {
  for (m_piter=m_convertH2S.begin();m_piter!=m_convertH2S.end();) {
    if (m_piter->second.second) {
      //cout<<"Active particle: "<<m_piter->first<<" "
      //	  <<m_piter->second.first->Flav()<<" : "
      //	  <<m_piter->second.first->ProductionBlob();
      //if (m_piter->second.first->ProductionBlob()) 
      // 	cout<<"("<<int(m_piter->second.first->ProductionBlob()->Type())<<")";
      //cout<<" , "<<m_piter->second.first->DecayBlob();
      //if (m_piter->second.first->DecayBlob()) 
      // 	cout<<"("<<int(m_piter->second.first->DecayBlob()->Type())<<")";
      //cout<<endl;
      m_piter1 = m_piter; m_piter1++; 
      delete m_piter->second.first;
      m_convertH2S.erase(m_piter);
      m_piter = m_piter1;
    }
    else {
      //cout<<"Passive particle: "<<m_piter->first<<" "
      //	  <<m_piter->second.first->Flav()<<" : "
      //  <<m_piter->second.first->ProductionBlob();
      //if (m_piter->second.first->ProductionBlob()) 
      //	cout<<"("<<int(m_piter->second.first->ProductionBlob()->Type())<<")";
      //cout<<" , "<<m_piter->second.first->DecayBlob();
      //if (m_piter->second.first->DecayBlob()) 
      //	cout<<"("<<int(m_piter->second.first->DecayBlob()->Type())<<")";
      //cout<<endl;
      m_piter++;
    }
  }
}
