#include "Beam_Remnant_Handler.H"

#include "Hadron_Remnant.H"
#include "Electron_Remnant.H"
#include "Photon_Remnant.H"
#include "No_Remnant.H"
#include "Data_Read.H"

#ifdef PROFILE__all
#define PROFILE__Beam_Remnant_Handler
#endif
#ifdef PROFILE__Beam_Remnant_Handler
#include "prof.hh" 
#else
#define PROFILE_HERE
#endif

using namespace SHERPA;

//#define DEBUG__Beam_Remnant_Handler

#ifndef DEBUG__Beam_Remnant_Handler
inline void SumMomenta(ATOOLS::Blob *bl) {}
#else
std::set<ATOOLS::Blob*> checked_blobs;

void SumMomenta(ATOOLS::Blob *bl,ATOOLS::Vec4D &inisum,
		ATOOLS::Vec4D & finsum,bool iterate=true) 
{
  // check if caught in loop
  std::set<ATOOLS::Blob*>::const_iterator bit=checked_blobs.find(bl);
  if (bit!=checked_blobs.end()) return;
  checked_blobs.insert(bl);
  for (int i=0;i<bl->NInP();++i) {
    ATOOLS::Particle * p =bl->InParticle(i);
    if (p->ProductionBlob()==NULL) inisum+=p->Momentum(); 
    else if (iterate) SumMomenta(p->ProductionBlob(),inisum,finsum);
    else inisum+=p->Momentum();
  }
  for (int i=0;i<bl->NOutP();++i) {
    ATOOLS::Particle * p =bl->OutParticle(i);
    if (p->DecayBlob()==NULL) finsum+=p->Momentum(); 
    else if (iterate) SumMomenta(p->DecayBlob(),inisum,finsum);
    else finsum+=p->Momentum();
  }
}

bool SumMomenta(ATOOLS::Blob *bl)
{
  ATOOLS::Vec4D inisum,finsum;
  checked_blobs.clear();
  SumMomenta(bl,inisum,finsum);
  bool test=inisum==finsum;
  if (!test) {
    ATOOLS::msg.Error()<<"SumMomenta(..): Summation does not agree."<<std::endl
		       <<"initial = "<<inisum<<" vs. final = "<<finsum<<std::endl;
  }
  return test;
}
#endif

Beam_Remnant_Handler::Beam_Remnant_Handler(std::string _m_path,std::string _m_file,
					   PDF::ISR_Handler * _p_isr,BEAM::Beam_Spectra_Handler * _p_beam,
					   double _m_scale) :
  p_isr(_p_isr), 
  p_beam(_p_beam),
  m_path(_m_path), 
  m_file(_m_file),
  m_maxtrials(1),
  m_fill(true)
{
  for (size_t i=0;i<2;++i) {
    if (p_isr->Flav(i).IsHadron()) p_beampart[i] = new Hadron_Remnant(p_isr,i,_m_scale);
    else if (p_isr->Flav(i).IsLepton()) p_beampart[i] = new Electron_Remnant(p_isr,i,_m_scale);
    else if (p_isr->Flav(i).IsPhoton()) p_beampart[i] = new Photon_Remnant(i);
    else p_beampart[i] = new No_Remnant(i);
  }
  for (size_t i=0;i<2;++i) p_beampart[i]->SetPartner(p_beampart[1-i]);
  p_kperp = new Primordial_KPerp(_m_path,_m_file);
}

Beam_Remnant_Handler::~Beam_Remnant_Handler() 
{  
  for (size_t i=0;i<2;++i) delete p_beampart[i];
  delete p_kperp;
}


bool Beam_Remnant_Handler::FillBunchBlobs(ATOOLS::Blob_List * bloblist,
					  ATOOLS::Particle_List * particlelist)
{
  p_bloblist=bloblist;
  p_particlelist=particlelist;
  ATOOLS::Blob_List::iterator endblob=bloblist->end(); 
  ATOOLS::Blob *blob;
  bool flag=false;
  ATOOLS::Particle *p;
  for (short unsigned int i=0;i<2;i++) {
    for (ATOOLS::Blob_List::iterator biter=bloblist->begin();biter!=endblob;++biter) {
      if ((*biter)->Status()==1 && (*biter)->Beam()==i && 
	  ((*biter)->Type()==ATOOLS::btp::Beam || 
	   (*biter)->Type()==ATOOLS::btp::IS_Shower)) {
	(*biter)->SetStatus(2);
	blob = new ATOOLS::Blob();
	bloblist->insert(bloblist->begin(),blob);
	blob->SetType(ATOOLS::btp::Bunch);
	blob->SetBeam(i);
	blob->SetId();
	blob->AddToOutParticles((*biter)->InParticle(0));
	if ((*biter)->InParticle(0)->Flav()==p_beam->GetBeam(i)->Beam() &&
	    ATOOLS::IsEqual((*biter)->InParticle(0)->E(),
			    p_beam->GetBeam(i)->InMomentum()[0])) {
	  p = new ATOOLS::Particle(*(*biter)->InParticle(0));
	  if (particlelist!=NULL) p->SetNumber(-particlelist->size());
	  else p->SetNumber(0);
	  blob->AddToInParticles(p);
	}
	else {
	  p = new ATOOLS::Particle(-1,p_beam->GetBeam(i)->Beam(),
				   p_beam->GetBeam(i)->InMomentum());
	  if (particlelist!=NULL) p->SetNumber(-particlelist->size());
	  else p->SetNumber(0);
	  p->SetStatus(2);
	  blob->AddToInParticles(p);
	  ATOOLS::Particle *p = new 
	    ATOOLS::Particle(-1,p_beam->GetBeam(i)->Remnant(),
			     p_beam->GetBeam(i)->InMomentum()+
			     (-1.)*(*biter)->InParticle(0)->Momentum());
	  if (particlelist!=NULL) p->SetNumber(-particlelist->size());
	  else p->SetNumber(0);
	  blob->AddToOutParticles(p);
	}
	flag=true;
      }
    }
  }
  return flag;
}

void Beam_Remnant_Handler::RemoveFromList(ATOOLS::Particle *cur)
{
  if (p_particlelist!=NULL) {
    for (ATOOLS::Particle_List::iterator pit=p_particlelist->begin();
	 pit!=p_particlelist->end();++pit) {
      if (*pit==cur) {
	p_particlelist->erase(pit--);
	return;
      }
    }
  }
}

void Beam_Remnant_Handler::RemoveFromList(ATOOLS::Blob *cur)
{
  for (ATOOLS::Blob_List::iterator bit=p_bloblist->begin();
       bit!=p_bloblist->end();++bit) {
    if (*bit==cur) {
      p_bloblist->erase(bit--);
      return;
    }
  }
}

bool Beam_Remnant_Handler::RemoveConnected(ATOOLS::Blob *blob)
{
  if (blob==NULL) return false;
  if (m_removed.find(blob)!=m_removed.end() ||
      blob->Type()==ATOOLS::btp::Beam) return false;
  m_removed.insert(blob);
  while (blob->NOutP()>0) {
    RemoveFromList(blob->OutParticle(0));
    if (!RemoveConnected(blob->OutParticle(0)->DecayBlob())) {
      blob->DeleteOutParticle(blob->OutParticle(0));
    }
  }
  while (blob->NInP()>0) {
    RemoveFromList(blob->InParticle(0));
    if (!RemoveConnected(blob->InParticle(0)->ProductionBlob())) {
      blob->DeleteInParticle(blob->InParticle(0));
    }
  }
  RemoveFromList(blob);
  delete blob;
  return true;
}

bool Beam_Remnant_Handler::RemoveLast()
{
  PROFILE_HERE;
  m_removed.clear();
  for (ATOOLS::Blob_List::iterator bit=p_bloblist->end()-1;
       bit!=p_bloblist->begin();--bit) {
    if ((*bit)->Type()==ATOOLS::btp::Hard_Collision) {
      RemoveConnected(*bit);
      return true;
    }
  }
  return false;
}

bool Beam_Remnant_Handler::FillBeamBlobs(ATOOLS::Blob_List *bloblist,
					 ATOOLS::Particle_List *particlelist)
{ 
  PROFILE_HERE;
  p_bloblist=bloblist;
  p_particlelist=particlelist;
  if (!m_fill) return false;
  ATOOLS::Blob_Iterator endblob = bloblist->end(); 
  ATOOLS::Blob * blob;
  bool okay=false, success=true, treat[2];
  for (size_t trials=0;trials<m_maxtrials;++trials) {
    success=false;
    for (short unsigned int i=0;i<2;++i) {
      p_beampart[i]->Clear();
      p_beampart[i]->ClearErrors();
      okay=false;
      treat[i]=false;
      for (ATOOLS::Blob_List::iterator biter=bloblist->begin();biter!=endblob;++biter) {
	if ((*biter)->Beam()==i && (*biter)->Type()==ATOOLS::btp::IS_Shower) { 
	  if (p_beampart[i]->Type()&rtp::qcd_remnant) {
	    if (trials==0 && !okay) {
	      blob = new ATOOLS::Blob();
	      bloblist->insert(bloblist->begin(),blob);
	      blob->SetId();
	      blob->SetType(ATOOLS::btp::Beam);
	      blob->SetBeam(i);
	      blob->SetStatus(1);
	      ATOOLS::Particle *p = new 
		ATOOLS::Particle(-1,p_isr->Flav(i),p_beam->GetBeam(i)->OutMomentum());
	      p->SetStatus(2);
	      blob->AddToInParticles(p);
	      p_beamblob[i]=blob;
	      okay=true;
	    }
	    p_beampart[i]->Extract((*biter)->InParticle(0));
	    if (trials>0) p_beamblob[i]->RemoveOutParticle((*biter)->InParticle(0));
	    (*biter)->SetStatus(2);
	    treat[i]=true;
	  }
	  else if (trials==0 && !ATOOLS::IsEqual((*biter)->InParticle(0)->E(),
						 p_beam->GetBeam(i)->OutMomentum()[0])) {
	    (*biter)->SetStatus(2);
	    blob = new ATOOLS::Blob();
	    blob->SetType(ATOOLS::btp::Beam);
	    blob->SetBeam(i);
	    blob->SetStatus(1);
	    p_beampart[i]->Extract((*biter)->InParticle(0));
	    if (trials>0) p_beamblob[i]->RemoveOutParticle((*biter)->InParticle(0));
	    ATOOLS::Particle *p = new 
	      ATOOLS::Particle(-1,p_isr->Flav(i),p_beam->GetBeam(i)->OutMomentum());
	    if (particlelist!=NULL) p->SetNumber(-particlelist->size());
	    else p->SetNumber(0);
	    p->SetStatus(2);
	    blob->AddToInParticles(p);
	    bloblist->insert(bloblist->begin(),blob);
	    blob->SetId();
	    p_beamblob[i]=blob;
	    okay=true;
	    treat[i]=true;
	  }
	}
      }
    }
    success=true;
    for (short unsigned int i=0;i<2;++i) {
      if (treat[i]) if (!p_beampart[i]->
			FillBlob(p_beamblob[i],particlelist)) success=false;
    }
    if (p_beampart[0]->Type()==rtp::hadron || 
	p_beampart[1]->Type()==rtp::hadron) {
      p_kperp->CreateKPerp(p_beamblob[0],p_beamblob[1]);
      for (short unsigned int i=0;i<2;++i) p_kperp->FillKPerp(p_beamblob[i]);
    }
    bool empty=false;
    for (short unsigned int i=0;i<2;++i) {
      if (!(success=success&&p_beampart[i]->AdjustKinematics())) {
	for (short unsigned int i=0; i<2;++i) p_beampart[i]->UnDo();
	if (!RemoveLast()) empty=true;
	break;
      }
    }
    if (empty || success) break;
  }
  SumMomenta(bloblist->front());
  if (!success) {
    msg_Tracking()<<"Beam_Remnant_Handler::FillBeamBlobs(..): "
		  <<"Not enough energy to preserve remnant masses."<<std::endl;
    for (size_t i=0;i<2;++i) p_beampart[i]->Clear();
    while (bloblist->size()>0) {
      delete *bloblist->begin();
      bloblist->erase(bloblist->begin());
    }
    return false;
  }
  return okay;
}

