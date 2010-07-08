#include "SHERPA/SoftPhysics/Beam_Remnant_Handler.H"

#include "PDF/Remnant/Hadron_Remnant.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"

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
  p_signal(NULL), p_shower(NULL),
  m_path(path), m_file(file), m_fill(true)
{
  Data_Reader read(" ",";","!","=");
  read.AddComment("#");
  read.SetInputPath(m_path);
  read.SetInputFile(m_file);
  if (!read.ReadFromFile(m_vmode,"BRH_VMODE")) m_vmode=0;
  else msg_Info()<<METHOD<<"(): Set check mode "<<m_vmode<<"."<<std::endl;
  p_kperp = new Primordial_KPerp(path,file);
  for (size_t i=0;i<2;++i) {
    p_beampart[i] = p_isr->GetRemnant(i);
    p_beampart[i]->SetBeam(beam->GetBeam(i));
    p_kperp->SetRemnant(p_beampart[i],i);
  }
}

Beam_Remnant_Handler::~Beam_Remnant_Handler() 
{  
  delete p_kperp;
}


Return_Value::code 
Beam_Remnant_Handler::FillBeamAndBunchBlobs(Blob_List *const bloblist)
{
  StorePerturbativeBlobs(bloblist);
  Return_Value::code fbc(FillBeamBlobs(bloblist));
  switch(fbc) {
  case Return_Value::Success:
    if (p_signal) delete p_signal; p_signal=NULL;
    if (p_shower) delete p_shower; p_shower=NULL;
    return FillBunchBlobs(bloblist);
  case Return_Value::Retry_Phase:
    bloblist->Clear();
    bloblist->push_back(p_signal);
    p_signal=NULL;
    bloblist->push_back(p_shower);
    p_shower=NULL;
    DEBUG_VAR(*bloblist);
    return Return_Value::Retry_ShowerEvent;
  default:
    if (p_signal) delete p_signal; p_signal=NULL;
    if (p_shower) delete p_shower; p_shower=NULL;
    return fbc;
  }
}

void Beam_Remnant_Handler::StorePerturbativeBlobs(ATOOLS::Blob_List* blobs)
{
  DEBUG_FUNC(*blobs);
  Blob* signal = blobs->FindFirst(btp::Signal_Process);
  if (signal) {
    p_signal = new Blob(signal);
  }

  Blob* shower = blobs->FindFirst(btp::Shower);
  if (shower) {
    p_shower = new Blob(shower);
    // delete copied particles because we have to take some of them from SP
    p_shower->DeleteInParticles();
    p_shower->DeleteOutParticles();

    Particle * part(NULL);
    for (int i=0;i<shower->NInP();i++) {
      if (shower->ConstInParticle(i)->ProductionBlob() &&
          shower->ConstInParticle(i)->ProductionBlob()->Type()==btp::Signal_Process) continue;
      part = new Particle((*shower->ConstInParticle(i)));
      p_shower->AddToInParticles(part);
    }
    for (int i=0;i<shower->NOutP();i++) {
      if (shower->ConstOutParticle(i)->DecayBlob() &&
          shower->ConstOutParticle(i)->DecayBlob()->Type()==btp::Signal_Process) continue;
      part = new Particle((*shower->ConstOutParticle(i)));
      part->SetDecayBlob(NULL);
      part->SetStatus(part_status::active);
      p_shower->AddToOutParticles(part);
    }

    p_shower->SetStatus(blob_status::needs_harddecays|
                        blob_status::needs_beams|
                        blob_status::needs_hadronization|
                        blob_status::needs_extraQED);
  }

  if (p_signal && p_shower) {
    for (int i=0; i<p_signal->NInP(); ++i) {
      p_shower->AddToOutParticles(p_signal->InParticle(i));
    }
    for (int i=0; i<p_signal->NOutP(); ++i) {
      p_shower->AddToInParticles(p_signal->OutParticle(i));
    }
  }
  DEBUG_VAR(*p_signal);
  DEBUG_VAR(*p_shower);
}

Return_Value::code Beam_Remnant_Handler::
FillBunchBlobs(Blob_List *const  bloblist,
	       Particle_List *const particlelist)
{
  PROFILE_HERE;
  for (Blob_List::iterator bit=bloblist->begin();
	 bit!=bloblist->end();++bit) {
    if ((*bit)->Type()==btp::Bunch) return Return_Value::Nothing;
  }
  bool flag(false);
  m_beam = 0;
  for (Blob_List::iterator bit=bloblist->begin();
       bit!=bloblist->end();++bit) {
    if ((*bit)->Has(blob_status::needs_beams) && 
	((*bit)->Type()==btp::Beam || (*bit)->Type()==btp::Shower)) {
      (*bit)->UnsetStatus(blob_status::needs_beams);
      bloblist->push_front
	(FillBunchBlob((*bit)->Beam(),(*bit)->InParticle(0)));
      if (m_beam>2) {
	msg_Error()<<"ERROR in "<<METHOD<<" : "<<std::endl
		   <<"   Too many bunch blobs required, "
		   <<"return 'Error' and hope for the best."<<std::endl;
	return Return_Value::Error;
      }
      flag=true;
    }
  }
  return (flag?Return_Value::Success:Return_Value::Nothing);
}

Return_Value::code Beam_Remnant_Handler::
FillBeamBlobs(Blob_List *const bloblist,
	      Particle_List *const particlelist)
{ 
  if (!m_fill) return Return_Value::Nothing;
  for (Blob_List::iterator bit=bloblist->begin();
	 bit!=bloblist->end();++bit) {
    if ((*bit)->Type()==btp::Beam) return Return_Value::Nothing;
  }
  for (short unsigned int i=0;i<2;++i) {
    p_beampart[i]->Clear();
    p_beampart[i]->ClearErrors();
    InitBeamBlob(i);
  }
  for (Blob_List::iterator bit=bloblist->begin();
       bit!=bloblist->end();++bit) {
    if ((*bit)->Has(blob_status::needs_beams) && 
	(*bit)->Type()==btp::Shower) {
      for (int i(0);i<(*bit)->NInP();++i) {
	Particle *isr_init((*bit)->InParticle(i));
	if (isr_init->ProductionBlob()!=NULL) continue;
	int beam((*bit)->Beam());
	if (beam<0) beam=isr_init->Momentum()[3]>0.0?0:1;
	if (isr_init->Flav().Strong() && 
	    isr_init->GetFlow(1)==0 && isr_init->GetFlow(2)==0) {
	  delete p_beamblob[beam]; 
	  p_beamblob[beam]=NULL; 
	  continue;
	}
	else {
	  (*bit)->AddStatus(blob_status::internal_flag);
	  if (!p_beampart[beam]->Extract(isr_init)) {
	    msg_Debugging()<<METHOD<<"(): Cannot extract "
			   <<*p_beamblob[beam]->InParticle(0)<<"\n  from"
			   <<*isr_init<<"\n  retry event "<<std::endl;
	    msg_Debugging()<<*bloblist<<std::endl;
	    for (short unsigned int i(0);i<2;++i) 
	      bloblist->push_front(p_beamblob[i]);
	    return Return_Value::Retry_Event;
	  }
	}
      }
      (*bit)->UnsetStatus(blob_status::needs_beams);
    }
  }
  for (short unsigned int i=0;i<2;++i) 
    if (p_beamblob[i]!=NULL) bloblist->push_front(p_beamblob[i]);
  if (p_beamblob[0]==NULL || p_beamblob[1]==NULL) {
    return Return_Value::Success;
  }
  for (short unsigned int i=0;i<2;++i) {
    if (!p_beampart[i]->FillBlob(p_beamblob[i],NULL)) {
      return Return_Value::Retry_Event;
    }
  }
  if (p_beampart[0]->Type()==PDF::rtp::hadron || 
      p_beampart[1]->Type()==PDF::rtp::hadron) {
    p_kperp->CreateKPerp(p_beamblob[0],p_beamblob[1]);
    for (short unsigned int i=0;i<2;++i) p_kperp->FillKPerp(p_beamblob[i]);
  }
  
  for (short unsigned int i=0;i<2;++i) {
    if (!p_beampart[i]->AdjustKinematics() || 
	!p_beampart[i]->AdjustColors()) { 
      return Return_Value::Retry_Phase;
    }
  }
  // set status for all partons entering the shower to decayed
  for (short unsigned int i=0;i<2;++i) {
    for (int j=0;j<p_beamblob[i]->NOutP();++j) {
      if (p_beamblob[i]->OutParticle(j)->DecayBlob()) {
	if (p_beamblob[i]->OutParticle(j)->DecayBlob()->Type() == btp::Shower) {
	  p_beamblob[i]->OutParticle(j)->SetStatus(part_status::decayed);
	}
      }
    }
  }
  if (bloblist->FourMomentumConservation() && bloblist->ColorConservation()) {
    for (Blob_List::iterator bit=bloblist->begin();
	 bit!=bloblist->end();++bit) {
      if ((*bit)->Has(blob_status::internal_flag) && 
	  (*bit)->Type()==btp::Shower) { 
	(*bit)->UnsetStatus(blob_status::internal_flag);
      }
    }
    return Return_Value::Success;
  }
  if (m_vmode) abort();
  return Return_Value::Retry_Event;
}



Blob * Beam_Remnant_Handler::FillBunchBlob(const int beam,Particle * particle) 
{
  Blob *blob = new Blob();
  blob->SetType(btp::Bunch);
  blob->SetBeam(beam);
  blob->SetId();
  blob->SetStatus(blob_status::needs_beams &
		  blob_status::needs_softUE &
		  blob_status::needs_hadronization);
  blob->AddToOutParticles(particle);
  if (particle->Flav()==p_beam->GetBeam(beam)->Beam() &&
      IsEqual(particle->E(),p_beam->GetBeam(beam)->InMomentum()[0])) {
    Particle *p = new Particle(*particle);
    p->SetNumber(0);
    blob->AddToInParticles(p);
  }
  else {
    Particle *p = new Particle(-1,p_beam->GetBeam(beam)->Beam(),
			       p_beam->GetBeam(beam)->InMomentum());
    p->SetNumber(0);
    p->SetStatus(part_status::decayed);
    p->SetFinalMass();
    blob->AddToInParticles(p);
    p = new Particle(-1,p_beam->GetBeam(beam)->Remnant(),
		     p_beam->GetBeam(beam)->InMomentum()-particle->Momentum());
    p->SetNumber(0);
    p->SetStatus(part_status::active);
    p->SetFinalMass();
    blob->AddToOutParticles(p);
  }
  m_beam++;
  return blob;
}


void Beam_Remnant_Handler::InitBeamBlob(const int beam) 
{
  p_beamblob[beam] = new Blob();
  p_beamblob[beam]->SetType(btp::Beam);
  p_beamblob[beam]->SetId();
  p_beamblob[beam]->SetBeam(beam);
  p_beamblob[beam]->SetStatus(blob_status::needs_beams |
			      blob_status::needs_softUE |
			      blob_status::needs_hadronization);
  Particle * beampart = new Particle(-1,p_isr->Flav(beam),
				     p_beam->GetBeam(beam)->OutMomentum());
  beampart->SetNumber(0);
  beampart->SetStatus(part_status::decayed);
  beampart->SetFinalMass();
  p_beamblob[beam]->AddToInParticles(beampart);
}

void Beam_Remnant_Handler::SetScale(const double scale)
{
  for (short unsigned int i=0;i<2;++i) p_beampart[i]->SetScale(scale);
}

void Beam_Remnant_Handler::CleanUp()
{
  for (short unsigned int i=0;i<2;++i) p_beampart[i]->Clear();
}
