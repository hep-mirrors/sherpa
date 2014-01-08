#include "SHERPA/SoftPhysics/Parametrised_Beam_Remnants.H"

#include "PDF/Remnant/Hadron_Remnant.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"

#ifdef PROFILE__all
#define PROFILE__Parametrised_Beam_Remnants
#endif
#ifdef PROFILE__Parametrised_Beam_Remnants
#include "prof.hh" 
#else
#define PROFILE_HERE
#endif

using namespace SHERPA;
using namespace ATOOLS;

Parametrised_Beam_Remnants::
Parametrised_Beam_Remnants(const std::string path,const std::string file,
			   PDF::ISR_Handler * const isr,
			   BEAM::Beam_Spectra_Handler *const beam):
  p_isr(isr), m_path(path), m_file(file)
{
  p_kperp = new Primordial_KPerp(path,file);
  for (size_t i=0;i<2;++i) {
    p_beampart[i] = p_isr->GetRemnant(i);
    p_beampart[i]->SetBeam(beam->GetBeam(i));
    p_kperp->SetRemnant(p_beampart[i],i);
  }
}

Parametrised_Beam_Remnants::~Parametrised_Beam_Remnants() 
{  
  delete p_kperp;
}


Return_Value::code Parametrised_Beam_Remnants::
FillBeamBlobs(Blob_List *const bloblist,
	      Particle_List *const particlelist)
{ 
  for (Blob_List::iterator bit=bloblist->begin();
       bit!=bloblist->end();++bit) {
    if ((*bit)->Has(blob_status::needs_beams) && 
	(*bit)->Type()==btp::Shower) {
      for (int i(0);i<(*bit)->NInP();++i) {
	Particle *isr_init((*bit)->InParticle(i));
	if (isr_init->ProductionBlob()!=NULL) continue;
	int beam((*bit)->Beam());
	if (beam<0) beam=isr_init->Beam();
	if (isr_init->Flav().Strong() && 
	    isr_init->GetFlow(1)==0 && isr_init->GetFlow(2)==0) {
	  delete p_beamblob[beam]; 
	  p_beamblob[beam]=NULL; 
	  continue;
	}
	else {
	  (*bit)->AddStatus(blob_status::internal_flag);
	  if (!p_beampart[beam]->Extract(isr_init)) {
	    msg_Debugging()<<METHOD<<"(): Cannot extract particle:\n"
			   <<*isr_init<<"\n  from \n"
			   <<*p_beamblob[beam]->InParticle(0)
			   <<"\n  retry event\n"<<*bloblist<<std::endl;
	    for (short unsigned int i(0);i<2;++i) 
	      bloblist->push_front(p_beamblob[i]);
	    return Return_Value::Retry_Event;
	  }
	}
      }
      (*bit)->UnsetStatus(blob_status::needs_beams);
    }
  }
  for (short unsigned int i=0;i<2;++i) {
    if (p_beamblob[i]!=NULL) bloblist->push_front(p_beamblob[i]);
  }
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
    for (short unsigned int i=0;i<2;++i) {
      //msg_Out()<<METHOD<<"("<<i<<"):\n"<<(*p_beamblob[i])<<"\n";
      p_kperp->FillKPerp(p_beamblob[i]);
    }
  }
  
  for (short unsigned int i=0;i<2;++i) {
    if (!p_beampart[i]->AdjustKinematics() || 
	!p_beampart[i]->AdjustColors()) { 
      //msg_Out()<<(*bloblist);
      return Return_Value::Retry_Event; 
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
  msg_Tracking()<<METHOD<<" found four momentum conservation error.\n";
  return Return_Value::New_Event;
}



