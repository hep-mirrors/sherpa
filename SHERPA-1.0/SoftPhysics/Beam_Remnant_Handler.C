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
    p_beampart[i] = p_isr->GetRemnant(i);
    p_beampart[i]->SetBeam(beam->GetBeam(i));
    p_kperp->SetRemnant(p_beampart[i],i);
  }
}

Beam_Remnant_Handler::~Beam_Remnant_Handler() 
{  
  for (size_t i=0;i<2;++i) delete p_beampart[i];
  delete p_kperp;
}


Return_Value::code 
Beam_Remnant_Handler::FillBeamAndBunchBlobs(Blob_List *const bloblist)
{
  switch (int(FillBeamBlobs(bloblist))) {
    case int(Return_Value::Nothing) : return Return_Value::Nothing;
    case int(Return_Value::Success) : break;
    case int(Return_Value::Error)   : return Return_Value::Error;
    default :
      msg.Error()<<"ERROR in "<<METHOD<<":"<<std::endl
		 <<"   FillBeamBlobs yields unknown outcome,"
		 <<" return 'Error' and hope for the best."<<std::endl;
      return Return_Value::Error;
  }
  switch (int(FillBunchBlobs(bloblist))) {
    case int(Return_Value::Nothing)     : return Return_Value::Nothing;
    case int(Return_Value::Success)     : break;
    case int(Return_Value::Retry_Event) : return Return_Value::Retry_Event;
    case int(Return_Value::Error)       : return Return_Value::Error;
    default :
      msg.Error()<<"ERROR in "<<METHOD<<":"<<std::endl
		 <<"   FillBunchBlobs yields unknown outcome,"
		 <<" return 'Error' and hope for the best."<<std::endl;
      return Return_Value::Error;
  }
  return Return_Value::Success;
}

Return_Value::code Beam_Remnant_Handler::
FillBunchBlobs(Blob_List *const  bloblist,
	       Particle_List *const particlelist)
{
  int bunchblobcount(0);
  for (Blob_List::iterator bit=bloblist->begin();
	 bit!=bloblist->end();++bit) {
    if ((*bit)->Type()==btp::Bunch) bunchblobcount++;
  }
  if (bunchblobcount>1) return Return_Value::Nothing;
  PROFILE_HERE;
  bool flag(false);
  m_beam = 0;
  for (Blob_List::iterator bit=bloblist->begin();
       bit!=bloblist->end();++bit) {
    if ((*bit)->Has(blob_status::needs_beams) && 
	((*bit)->Type()==btp::Beam || (*bit)->Type()==btp::IS_Shower)) {
      (*bit)->UnsetStatus(blob_status::needs_beams);
      bloblist->push_front(FillBunchBlob((*bit)->Beam(),(*bit)->InParticle(0)));
      if (m_beam>2) {
	msg.Error()<<"ERROR in "<<METHOD<<" : "<<std::endl
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
  PROFILE_HERE;
  if (!m_fill) return Return_Value::Nothing;
  if (p_mehandler==NULL) {
    p_mehandler=GET_OBJECT(Matrix_Element_Handler,"ME_Handler");
    if (p_mehandler==NULL) 
      THROW(fatal_error,"No matrix element handler found.");
  }
  int beamblobcount(0);
  for (Blob_List::iterator bit=bloblist->begin();
	 bit!=bloblist->end();++bit) {
    if ((*bit)->Type()==btp::Beam) beamblobcount++;
  }
  if (beamblobcount>1) return Return_Value::Nothing;
	
  p_beamblob[1]=p_beamblob[0]=NULL;
  Blob_List::iterator endblob=bloblist->end(); 
  for (short unsigned int i=0;i<2;++i) {
    p_beampart[i]->Clear();
    p_beampart[i]->ClearErrors();
  }

  Particle * isr_init;
  short unsigned int beam(0);
  for (Blob_List::iterator bit=bloblist->begin();
       bit!=bloblist->end();++bit) {
    if ((*bit)->Has(blob_status::needs_beams) && (*bit)->Type()==btp::IS_Shower) { 
      isr_init = (*bit)->InParticle(0);
      beam     = (*bit)->Beam();
      if (isr_init->Flav().Strong() && isr_init->GetFlow(1)==0 && isr_init->GetFlow(2)==0) {
	delete p_beamblob[beam]; p_beamblob[beam]=NULL; continue;
      }
      else {
	(*bit)->UnsetStatus(blob_status::needs_beams);
	(*bit)->AddStatus(blob_status::internal_flag);
	if (!FillBeamBlob(beam,isr_init)) {
	  msg.Error()<<"ERROR in "<<METHOD<<":"<<std::endl
		     <<"   FillBeamBlob failed for\n   "
		     <<p_beamblob[beam]->InParticle(0)<<" --> \n"
		     <<isr_init<<".\n   Retry event "
		     <<rpa.gen.NumberOfDicedEvents()<<"."<<std::endl;
	  msg_Tracking()<<*bloblist<<std::endl;
	  Reset(bloblist);
	  return Return_Value::Retry_Event;
	}
	bloblist->push_front(p_beamblob[beam]);
      }
    }
  }
  if (p_beamblob[0]==NULL || p_beamblob[1]==NULL) {
    Reset(bloblist);
    return Return_Value::Success;
  }

  if (p_beampart[0]->Type()==PDF::rtp::hadron || 
      p_beampart[1]->Type()==PDF::rtp::hadron) {
    p_kperp->CreateKPerp(p_beamblob[0],p_beamblob[1]);
    for (short unsigned int i=0;i<2;++i) p_kperp->FillKPerp(p_beamblob[i]);
  }

  for (short unsigned int i=0;i<2;++i) {
    if (!p_beampart[i]->AdjustKinematics() || !p_beampart[i]->AdjustColors()) { 
      Reset(bloblist); 
      return Return_Value::Retry_Event; 
    }
  }
    
  if (bloblist->FourMomentumConservation() && bloblist->ColorConservation()) {
    p_mehandler->ResetNumberOfTrials();
    for (Blob_List::iterator bit=bloblist->begin();bit!=bloblist->end();++bit) {
      if ((*bit)->Has(blob_status::internal_flag) && (*bit)->Type()==btp::IS_Shower) { 
	(*bit)->UnsetStatus(blob_status::internal_flag);
      }
    }
    return Return_Value::Success;
  }
  Reset(bloblist);
  return Return_Value::Retry_Event;
}



Blob * Beam_Remnant_Handler::FillBunchBlob(const int beam,Particle * particle) {
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


bool Beam_Remnant_Handler::FillBeamBlob(const int beam,Particle * particle) {
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
  //std::cout<<"Extract "<<isr_init->Flav()<<" out of "
  //	 <<p_beamblob[beam]->InParticle(0)->Flav()<<std::endl;
  if (!p_beampart[beam]->Extract(particle) ||
      !p_beampart[beam]->FillBlob(p_beamblob[beam],NULL)) return false;
  return true;
}

void Beam_Remnant_Handler::Reset(Blob_List * const bloblist) {
  while (bloblist->front()==p_beamblob[0] || bloblist->front()==p_beamblob[1]) {
    bloblist->pop_front();
  }
  for (short unsigned int beam=0;beam<2;beam++) {
    if (p_mehandler->Weight()!=1.0) p_mehandler->SaveNumberOfTrials();
    if (p_beamblob[beam]) { delete p_beamblob[beam]; p_beamblob[beam]=NULL; }
  }
  for (Blob_List::iterator bit=bloblist->begin();
       bit!=bloblist->end();++bit) {
    if ((*bit)->Has(blob_status::internal_flag) && (*bit)->Type()==btp::IS_Shower) { 
	(*bit)->SetStatus(blob_status::needs_beams);
	(*bit)->UnsetStatus(blob_status::internal_flag);
    }
  }
}

void Beam_Remnant_Handler::SetScale(const double scale)
{
  for (short unsigned int i=0;i<2;++i) p_beampart[i]->SetScale(scale);
}

