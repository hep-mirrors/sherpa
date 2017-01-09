#include "SHERPA/SoftPhysics/Beam_Remnant_Handler.H"
#include "ATOOLS/Org/Default_Reader.H"

using namespace SHERPA;
using namespace ATOOLS;

Beam_Remnant_Handler::
Beam_Remnant_Handler(const std::string path,const std::string file,
		     BEAM::Beam_Spectra_Handler *const beam,
		     PDF::ISR_Handler *const isr,
		     Soft_Collision_Handler *const softcollisions):
  p_parametrised(NULL), 
  p_shrimps(softcollisions?softcollisions->GetShrimps():NULL),
  p_beam(beam), m_fill(1)
{
  Default_Reader reader;
  reader.SetInputPath(path);
  reader.SetInputFile(file);
  m_fill  = reader.Get("BEAM_REMNANTS", 1, "remnants", METHOD);
  m_vmode = reader.Get("BRH_VMODE", 0, "check mode", METHOD);
  if (p_shrimps==NULL) {
    p_parametrised = new Parametrised_Beam_Remnants(path,file,isr,p_beam);
    p_parametrised->SetScale(4.0);
  }
}

Beam_Remnant_Handler::~Beam_Remnant_Handler() 
{  
  if (p_parametrised) delete p_parametrised;
}


Return_Value::code 
Beam_Remnant_Handler::FillBeamAndBunchBlobs(Blob_List *const bloblist)
{
  if (!m_fill) {
    bool set(false);
    for (Blob_List::iterator bit=bloblist->begin();
	 bit!=bloblist->end();++bit) {
      if ((*bit)->Has(blob_status::needs_beams)) {
	(*bit)->UnsetStatus(blob_status::needs_beams);
	(*bit)->UnsetStatus(blob_status::internal_flag);
	set=true;
      }
    }
    if (!set) return Return_Value::Nothing;
    if (bloblist->FourMomentumConservation())
      return Return_Value::Success;
    msg_Tracking()<<METHOD<<" found four momentum conservation error.\n";
    if (m_vmode) THROW(fatal_error,"Four Momentum not conserved.");
    return Return_Value::New_Event;
  }
  for (Blob_List::iterator bit=bloblist->begin();
	 bit!=bloblist->end();++bit) {
    if ((*bit)->Type()==btp::Beam) {
      return Return_Value::Nothing;
    }
  }
  for (short unsigned int i=0;i<2;++i) {
    p_beamblobs[i] = InitBeamBlob(i);
    if (p_shrimps)      p_shrimps->SetBeamBlob(p_beamblobs[i],i);
    if (p_parametrised) p_parametrised->SetBeamBlob(p_beamblobs[i],i);
  }
  Return_Value::code fbc(Return_Value::Error);
  if (p_shrimps) {
    fbc =  p_shrimps->FillBeamBlobs(bloblist);
  }
  else if (p_parametrised) {
    fbc = p_parametrised->FillBeamBlobs(bloblist);
    if (fbc==Return_Value::New_Event && m_vmode)
      THROW(fatal_error,"Four Momentum not conserved.");
  } else {
    THROW(fatal_error, "No physics model initialised.");
  }
  if (fbc!=Return_Value::Success) return fbc;
  fbc = FillBunchBlobs(bloblist);
  return fbc;
}

Return_Value::code Beam_Remnant_Handler::
FillBunchBlobs(Blob_List *const  bloblist,
	       Particle_List *const particlelist)
{
  for (Blob_List::iterator bit=bloblist->begin();
	 bit!=bloblist->end();++bit) {
    if ((*bit)->Type()==btp::Bunch) return Return_Value::Nothing;
  }
  bool flag(false);
  m_beam = 0;
  Blob * bunch;
  for (Blob_List::iterator bit=bloblist->begin();
       bit!=bloblist->end();++bit) {
    if ((*bit)->Has(blob_status::needs_beams) && 
	((*bit)->Type()==btp::Beam || (*bit)->Type()==btp::Shower)) {
      (*bit)->UnsetStatus(blob_status::needs_beams);
      bunch = FillBunchBlob((*bit)->Beam(),(*bit)->InParticle(0));
      bloblist->push_front(bunch);
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


ATOOLS::Blob * Beam_Remnant_Handler::InitBeamBlob(const int beam) 
{
  ATOOLS::Blob * blob = new Blob();
  blob->SetType(btp::Beam);
  blob->SetId();
  blob->SetBeam(beam);
  blob->SetStatus(blob_status::needs_beams |
		  blob_status::needs_softUE |
		  blob_status::needs_hadronization);
  Particle * beampart = new Particle(-1,p_beam->GetBeam(beam)->Beam(),
				     p_beam->GetBeam(beam)->OutMomentum());
  beampart->SetNumber(0);
  beampart->SetBeam(beam);
  beampart->SetStatus(part_status::decayed);
  beampart->SetFinalMass();
  blob->AddToInParticles(beampart);
  return blob;
}

void Beam_Remnant_Handler::CleanUp(const size_t & mode)
{
  if (p_shrimps) {
    p_shrimps->CleanUp(mode);
  }
  else p_parametrised->CleanUp();
}
