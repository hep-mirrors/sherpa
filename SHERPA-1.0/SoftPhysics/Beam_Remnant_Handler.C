#include "Beam_Remnant_Handler.H"

#include "Hadron_Remnant.H"
#include "Electron_Remnant.H"
#include "Photon_Remnant.H"
#include "No_Remnant.H"
#include "Data_Reader.H"
#include "Run_Parameter.H"
#include "Matrix_Element_Handler.H"
#include "Exception.H"

#ifdef PROFILE__all
#define PROFILE__Beam_Remnant_Handler
#endif
#ifdef PROFILE__Beam_Remnant_Handler
#include "prof.hh" 
#else
#define PROFILE_HERE
#endif

using namespace SHERPA;

Beam_Remnant_Handler::
Beam_Remnant_Handler(const std::string path,const std::string file,
		     PDF::ISR_Handler *const isr,
		     BEAM::Beam_Spectra_Handler *const beam):
  p_isr(isr), p_beam(beam), 
  p_mehandler(NULL), m_path(path), m_file(file), m_fill(true)
{
  for (size_t i=0;i<2;++i) {
    if (p_isr->Flav(i).IsHadron()) {
      Hadron_Remnant *remnant = new Hadron_Remnant(p_isr,i);
      ATOOLS::Data_Reader reader;
      reader.SetInputPath(m_path);
      reader.SetInputFile(m_file);
      double helpd;
      if (!reader.ReadFromFile(helpd,"REMNANT_STRING_TENSION")) helpd=1.0;
      remnant->SetStringDrawing(helpd,0);
      if (!reader.ReadFromFile(helpd,"REMNANT_RANDOM_STRINGS")) helpd=0.0;
      remnant->SetStringDrawing(helpd,1);
      p_beampart[i]=remnant;
    }
    else if (p_isr->Flav(i).IsLepton()) 
      p_beampart[i] = new Electron_Remnant(p_isr,i);
    else if (p_isr->Flav(i).IsPhoton()) 
      p_beampart[i] = new Photon_Remnant(i);
    else p_beampart[i] = new No_Remnant(i);
    p_beampart[i]->SetBeamEnergy(beam->GetBeam(i)->Energy());
  }
  for (size_t i=0;i<2;++i) p_beampart[i]->SetPartner(p_beampart[1-i]);
  p_kperp = new Primordial_KPerp(path,file);
  for (size_t i=0;i<2;++i) p_kperp->SetRemnant(p_beampart[i],i);
}

Beam_Remnant_Handler::~Beam_Remnant_Handler() 
{  
  for (size_t i=0;i<2;++i) delete p_beampart[i];
  delete p_kperp;
}

bool Beam_Remnant_Handler::
FillBunchBlobs(ATOOLS::Blob_List *const  bloblist,
	       ATOOLS::Particle_List *const particlelist)
{
  p_bloblist=bloblist;
  p_particlelist=particlelist;
  ATOOLS::Blob_List::iterator endblob=bloblist->end(); 
  ATOOLS::Blob *blob;
  bool flag=false;
  ATOOLS::Particle *p;
  for (short unsigned int i=0;i<2;i++) {
    for (ATOOLS::Blob_List::iterator biter=bloblist->begin();
	 biter!=endblob;++biter) {
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
	if ((*biter)->InParticle(0)->Number()<0) (*biter)->InParticle(0)->SetNumber(0);
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

bool Beam_Remnant_Handler::
FillBeamBlobs(ATOOLS::Blob_List *const bloblist,
	      ATOOLS::Particle_List *const particlelist)
{ 
  PROFILE_HERE;
  p_bloblist=bloblist;
  p_particlelist=particlelist;
  if (p_mehandler==NULL) {
    p_mehandler=GET_OBJECT(Matrix_Element_Handler,"ME_Handler");
    if (p_mehandler==NULL) THROW(fatal_error,"No matrix element handler found.");
  }
  if (!m_fill) return false;
  bool adjusted=false;
  bool okay=false;
  while (!adjusted) {
    ATOOLS::Blob_List::iterator endblob = bloblist->end(); 
    ATOOLS::Blob * blob;
    bool treat[2];
    for (short unsigned int i=0;i<2;++i) {
      p_beampart[i]->Clear();
      p_beampart[i]->ClearErrors();
      okay=false;
      treat[i]=false;
      p_beamblob[i]=NULL;
      for (ATOOLS::Blob_List::iterator biter=bloblist->begin();
	   biter!=endblob;++biter) {
	if ((*biter)->Beam()==i && (*biter)->Type()==ATOOLS::btp::IS_Shower) { 
	  if (p_beampart[i]->Type()&rtp::qcd_remnant) {
	    if (!okay) {
	      blob = new ATOOLS::Blob();
	      bloblist->insert(bloblist->begin(),blob);
	      blob->SetId();
	      blob->SetType(ATOOLS::btp::Beam);
	      blob->SetBeam(i);
	      blob->SetStatus(1);
	      ATOOLS::Particle *p = 
		new ATOOLS::Particle(-1,p_isr->Flav(i),
				     p_beam->GetBeam(i)->OutMomentum());
	      p->SetStatus(2);
	      blob->AddToInParticles(p);
	      p_beamblob[i]=blob;
	      okay=true;
	    }
	    if (!p_beampart[i]->Extract((*biter)->InParticle(0))) {
	      ATOOLS::msg.Error()<<"Beam_Remnant_Handler::FillBeamBlobs(..): "
				 <<"Extract parton failed.\n   Retry event "
				 <<ATOOLS::rpa.gen.NumberOfDicedEvents()
				 <<"."<<std::endl;
	      msg_Tracking()<<*bloblist<<std::endl;
	      bloblist->Clear();
	      if (p_mehandler->Weight()!=1.) p_mehandler->SaveNumberOfTrials();
	      return false;
	    }
	    (*biter)->SetStatus(2);
	    treat[i]=true;
	  }
	  else if (!ATOOLS::IsEqual((*biter)->InParticle(0)->E(),
				    p_beam->GetBeam(i)->OutMomentum()[0])) {
	    (*biter)->SetStatus(2);
	    blob = new ATOOLS::Blob();
	    blob->SetType(ATOOLS::btp::Beam);
	    blob->SetBeam(i);
	    blob->SetStatus(1);
	    p_beampart[i]->Extract((*biter)->InParticle(0));
	    ATOOLS::Particle *p = new 
	      ATOOLS::Particle(-1,p_isr->Flav(i),
			       p_beam->GetBeam(i)->OutMomentum());
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
    if (p_beamblob[0]==NULL || p_beamblob[1]==NULL) return false;
    for (short unsigned int i=0;i<2;++i) 
      if (p_beamblob[i]) if (!p_beampart[i]->FillBlob(p_beamblob[i],particlelist)) {
	ATOOLS::msg.Error()<<*bloblist<<std::endl;
	if (i==0) p_beampart[1]->FillBlob(p_beamblob[i],particlelist);
	bloblist->Clear();
	if (p_mehandler->Weight()!=1.) p_mehandler->SaveNumberOfTrials();
	return false;
      }
    if (p_beampart[0]->Type()==rtp::hadron || 
	p_beampart[1]->Type()==rtp::hadron) {
      p_kperp->CreateKPerp(p_beamblob[0],p_beamblob[1]);
      for (short unsigned int i=0;i<2;++i) p_kperp->FillKPerp(p_beamblob[i]);
    }
    adjusted=true;
    for (short unsigned int i=0;i<2;++i) 
      if (!p_beampart[i]->AdjustKinematics()) adjusted=false;
    if (!adjusted) {
      ATOOLS::Blob *lastmi(bloblist->FindLast(ATOOLS::btp::Hard_Collision));
      if (lastmi!=NULL) {
	for (short unsigned int i=0;i<2;++i) bloblist->Delete(p_beamblob[i]);
	bloblist->DeleteConnected(lastmi);
      }
      else break;
    }
  }
  if (!bloblist->FourMomentumConservation()) {
    msg_Info()<<"Beam_Remnant_Handler::FillBeamBlobs(..): Retry event "
	      <<ATOOLS::rpa.gen.NumberOfDicedEvents()<<".\n"
	      <<*bloblist<<std::endl;
    bloblist->Clear();
  }
  p_mehandler->ResetNumberOfTrials();
  return okay;
}

void Beam_Remnant_Handler::SetScale(const double scale)
{
  for (short unsigned int i=0;i<2;++i) p_beampart[i]->SetScale(scale);
}
