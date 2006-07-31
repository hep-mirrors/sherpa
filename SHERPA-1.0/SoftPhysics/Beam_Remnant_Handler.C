#include "Beam_Remnant_Handler.H"

#include "Hadron_Remnant.H"
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
using namespace ATOOLS;

Beam_Remnant_Handler::
Beam_Remnant_Handler(const std::string path,const std::string file,
		     PDF::ISR_Handler *const isr,
		     BEAM::Beam_Spectra_Handler *const beam):
  p_isr(isr), p_beam(beam), 
  p_mehandler(NULL), m_path(path), m_file(file), m_fill(true)
{
  p_kperp = new Primordial_KPerp(path,file);
  for (size_t i=0;i<2;++i) {
    p_beampart[i]=p_isr->GetRemnant(i);
    p_beampart[i]->SetBeam(beam->GetBeam(i));
    p_kperp->SetRemnant(p_beampart[i],i);
  }
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
  PROFILE_HERE;
  p_bloblist=bloblist;
  p_particlelist=particlelist;
  bool flag=false;
  for (short unsigned int i=0;i<2;i++) {
    for (Blob_List::iterator bit=bloblist->begin();
	 bit!=bloblist->end();++bit) {
      if ((*bit)->Status()==1 && (*bit)->Beam()==i && 
	  ((*bit)->Type()==btp::Beam || (*bit)->Type()==btp::IS_Shower)) {
	(*bit)->SetStatus(2);
	Blob *blob = new Blob();
	bloblist->push_front(blob);
	blob->SetType(btp::Bunch);
	blob->SetBeam(i);
	blob->SetId();
	blob->AddToOutParticles((*bit)->InParticle(0));
	if ((*bit)->InParticle(0)->Flav()==p_beam->GetBeam(i)->Beam() &&
	    IsEqual((*bit)->InParticle(0)->E(),
		    p_beam->GetBeam(i)->InMomentum()[0])) {
	  Particle *p = new Particle(*(*bit)->InParticle(0));
	  p->SetNumber(0);
	  blob->AddToInParticles(p);
	}
	else {
	  Particle *p = new Particle(-1,p_beam->GetBeam(i)->Beam(),
				     p_beam->GetBeam(i)->InMomentum());
	  p->SetNumber(0);
	  p->SetStatus(part_status::decayed);
	  p->SetFinalMass();
	  blob->AddToInParticles(p);
	  p = new Particle(-1,p_beam->GetBeam(i)->Remnant(),
			   p_beam->GetBeam(i)->InMomentum()-
			   (*bit)->InParticle(0)->Momentum());
	  p->SetNumber(0);
	  p->SetStatus(part_status::active);
	  p->SetFinalMass();
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
    if (p_mehandler==NULL) 
      THROW(fatal_error,"No matrix element handler found.");
  }
  if (!m_fill) return false;
  p_beamblob[1]=p_beamblob[0]=NULL;
  Blob_List::iterator endblob=bloblist->end(); 
  for (short unsigned int i=0;i<2;++i) {
    p_beampart[i]->Clear();
    p_beampart[i]->ClearErrors();
    for (Blob_List::iterator bit=bloblist->begin();
	 bit!=bloblist->end();++bit) {
      if ((*bit)->Beam()==i && (*bit)->Type()==btp::IS_Shower &&
	  (*bit)->TypeSpec()!="ADICIC++0.0") {
	ATOOLS::Particle *in=(*bit)->InParticle(0);
	if (in->Flav().Strong() &&
	    in->GetFlow(1)==0 && in->GetFlow(2)==0) continue;
	if (p_beamblob[i]==NULL) {
	    p_beamblob[i] = new Blob();
	    p_beamblob[i]->SetType(btp::Beam);
	    bloblist->push_front(p_beamblob[i]);
	    p_beamblob[i]->SetId();
	    p_beamblob[i]->SetBeam(i);
	    p_beamblob[i]->SetStatus(1);
	    Particle *p = new Particle(-1,p_isr->Flav(i),
				       p_beam->GetBeam(i)->OutMomentum());
	    p->SetNumber(0);
	    p->SetStatus(part_status::decayed);
	    p->SetFinalMass();
	    p_beamblob[i]->AddToInParticles(p);
	}
	if (!p_beampart[i]->Extract((*bit)->InParticle(0))) {
	  msg.Error()<<"Beam_Remnant_Handler::FillBeamBlobs(..): "
		     <<"Extract parton failed for\n   "
		     <<*(*bit)->InParticle(0)<<".\n   Retry event "
		     <<rpa.gen.NumberOfDicedEvents()<<"."<<std::endl;
	  msg_Tracking()<<*bloblist<<std::endl;
	  bloblist->Clear();
	  if (p_mehandler->Weight()!=1.0) p_mehandler->SaveNumberOfTrials();
	  return false;
	}
	(*bit)->SetStatus(2);
      }
      if((*bit)->Beam()==0 &&
	 (*bit)->Type()==btp::IS_Shower &&
	 (*bit)->TypeSpec()=="ADICIC++0.0") {
	//Special treatment for ADICIC IS Blobs!
	Particle* in=(*bit)->InParticle(i);
	//if(in->Flav().Strong() &&
	//   in->GetFlow(1)==0 && in->GetFlow(2)==0) continue;
	if(p_beamblob[i]==NULL) {
	  p_beamblob[i]=new Blob();
	  p_beamblob[i]->SetType(btp::Beam);
	  bloblist->push_front(p_beamblob[i]);
	  p_beamblob[i]->SetId();
	  p_beamblob[i]->SetBeam(i);
	  p_beamblob[i]->SetStatus(1);
	  Particle* p=new Particle(-1,p_isr->Flav(i),
				   p_beam->GetBeam(i)->OutMomentum());
	  p->SetNumber(0);
	  p->SetFinalMass();
	  p->SetStatus(part_status::decayed);
	  p_beamblob[i]->AddToInParticles(p);
	}
	if(!p_beampart[i]->Extract(in)) {
	  msg.Error()<<"Beam_Remnant_Handler::FillBeamBlobs(..): "
		     <<"Extract parton failed for\n   "
		     <<*(*bit)->InParticle(i)<<".\n   Retry event "
		     <<rpa.gen.NumberOfDicedEvents()<<"."<<std::endl;
	  msg_Tracking()<<*bloblist<<std::endl;
	  bloblist->Clear();
	  if(p_mehandler->Weight()!=1.0) p_mehandler->SaveNumberOfTrials();
	  return false;
	}
	(*bit)->SetStatus(2);
      }
    }
  }
  if (p_beamblob[0]==NULL || p_beamblob[1]==NULL) {
    if (bloblist->FourMomentumConservation()) {
      p_mehandler->ResetNumberOfTrials();
      return true;
    }
  }
  for (short unsigned int i=0;i<2;++i) 
    if (p_beamblob[i]) 
      if (!p_beampart[i]->FillBlob(p_beamblob[i],particlelist)) {
	msg_Tracking()<<*bloblist<<std::endl;
	if (i==0) p_beampart[1]->FillBlob(p_beamblob[i],particlelist);
	bloblist->Clear();
	if (p_mehandler->Weight()!=1.0) p_mehandler->SaveNumberOfTrials();
	return false;
      }
  if (p_beampart[0]->Type()==PDF::rtp::hadron || 
      p_beampart[1]->Type()==PDF::rtp::hadron) {
    p_kperp->CreateKPerp(p_beamblob[0],p_beamblob[1]);
    for (short unsigned int i=0;i<2;++i) p_kperp->FillKPerp(p_beamblob[i]);
  }
  bool adjusted=true;
  for (short unsigned int i=0;i<2;++i) {
    if (!p_beampart[i]->AdjustKinematics()) adjusted=false;
    if (!p_beampart[i]->AdjustColors()) adjusted=false;
  }
  if (!adjusted) {
    bloblist->Clear();
    if (p_mehandler->Weight()!=1.0) p_mehandler->SaveNumberOfTrials();
    return false;
  }
  if (bloblist->FourMomentumConservation() && bloblist->ColorConservation()) {
    p_mehandler->ResetNumberOfTrials();
    return true;
  }
  ATOOLS::msg.Error()<<"Beam_Remnant_Handler::FillBeamBlobs(..): Retry event "
		     <<rpa.gen.NumberOfDicedEvents()<<"."<<std::endl;
  bloblist->Clear();
  return false;
}

void Beam_Remnant_Handler::SetScale(const double scale)
{
  for (short unsigned int i=0;i<2;++i) p_beampart[i]->SetScale(scale);
}

