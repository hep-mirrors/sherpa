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
using namespace ATOOLS;

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
      Data_Reader reader;
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
	  p->SetStatus(2);
	  blob->AddToInParticles(p);
	  p = new Particle(-1,p_beam->GetBeam(i)->Remnant(),
			   p_beam->GetBeam(i)->InMomentum()-
			   (*bit)->InParticle(0)->Momentum());
	  p->SetNumber(0);
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
      if ((*bit)->Beam()==i && (*bit)->Type()==btp::IS_Shower) { 
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
	    p->SetStatus(2);
	    p_beamblob[i]->AddToInParticles(p);
	}
	if (p_beampart[i]->Type()&rtp::qcd_remnant ||
	    !IsEqual((*bit)->InParticle(0)->E(),
		     p_beam->GetBeam(i)->OutMomentum()[0])) {
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
	}
	(*bit)->SetStatus(2);
      }
    }
  }
  if (p_beamblob[0]==NULL || p_beamblob[1]==NULL) return false;
  for (short unsigned int i=0;i<2;++i) 
    if (p_beamblob[i]) 
      if (!p_beampart[i]->FillBlob(p_beamblob[i],particlelist)) {
	msg.Error()<<*bloblist<<std::endl;
	if (i==0) p_beampart[1]->FillBlob(p_beamblob[i],particlelist);
	bloblist->Clear();
	if (p_mehandler->Weight()!=1.0) p_mehandler->SaveNumberOfTrials();
	return false;
      }
  if (p_beampart[0]->Type()==rtp::hadron || 
      p_beampart[1]->Type()==rtp::hadron) {
    p_kperp->CreateKPerp(p_beamblob[0],p_beamblob[1]);
    for (short unsigned int i=0;i<2;++i) p_kperp->FillKPerp(p_beamblob[i]);
  }
  bool adjusted=true;
  for (short unsigned int i=0;i<2;++i) 
    if (!p_beampart[i]->AdjustKinematics()) adjusted=false;
  if (!adjusted) {
    bloblist->Clear();
    if (p_mehandler->Weight()!=1.0) p_mehandler->SaveNumberOfTrials();
    return false;
  }
  if (!bloblist->FourMomentumConservation()) {
    msg_Info()<<"Beam_Remnant_Handler::FillBeamBlobs(..): Retry event "
	      <<rpa.gen.NumberOfDicedEvents()<<".\n"
	      <<*bloblist<<std::endl;
    bloblist->Clear();
    return false;
  }
  p_mehandler->ResetNumberOfTrials();
  return true;
}

void Beam_Remnant_Handler::SetScale(const double scale)
{
  for (short unsigned int i=0;i<2;++i) p_beampart[i]->SetScale(scale);
}
