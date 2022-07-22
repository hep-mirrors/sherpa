#include "SHERPA/SoftPhysics/Beam_Remnant_Handler.H"
#include "BEAM/Main/Beam_Spectra_Handler.H"
#include "PDF/Main/ISR_Handler.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Org/Exception.H"

using namespace SHERPA;
using namespace ATOOLS;

Beam_Remnant_Handler::
Beam_Remnant_Handler(BEAM::Beam_Spectra_Handler *const beam,
		     REMNANTS::Remnant_Handler *const remnants,
		     Soft_Collision_Handler *const softcollisions):
  p_remnants(remnants),
  p_shrimps(softcollisions?softcollisions->GetShrimps():NULL),
  p_beam(beam), m_fill(true), m_vmode(false)
{
  /*
    if (p_shrimps) {
    m_name = std::string("Shrimps");
    return;
    }
  */
  p_shrimps = NULL;
  Settings& s = Settings::GetMainSettings();
  m_fill  = s["BEAM_REMNANTS"].SetDefault(true).Get<bool>();
  m_vmode = s["BRH_VMODE"].SetDefault(false).Get<bool>();
  p_remnants->SetScale2(sqr(4.0));
  m_name = std::string("Parametrised");
}

Beam_Remnant_Handler::~Beam_Remnant_Handler() {}

Return_Value::code
Beam_Remnant_Handler::FillBeamAndBunchBlobs(Blob_List *const bloblist,
					    const bool & onlyBunch)
{
  if (!m_fill) return TreatNoFill(bloblist);
  Return_Value::code fbc(Return_Value::Nothing);
  for (Blob_List::iterator bit=bloblist->begin();
       bit!=bloblist->end();++bit) {
    if ((*bit)->Type()==btp::Beam) return fbc;
  }
  if (!onlyBunch) {
    if (p_shrimps) fbc = p_shrimps->MakeBeamBlobs(bloblist);
    else           fbc = p_remnants->MakeBeamBlobs(bloblist);
  }
  if (fbc==Return_Value::New_Event && m_vmode)
    THROW(fatal_error,"Four Momentum not conserved.");
  if (fbc!=Return_Value::Success) return fbc;
  return FillBunchBlobs(bloblist);
}


Return_Value::code 
Beam_Remnant_Handler::TreatNoFill(Blob_List *const bloblist)
{
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
  if (bloblist->FourMomentumConservation()) return Return_Value::Success;
  msg_Tracking()<<METHOD<<" found four momentum conservation error.\n";
  if (m_vmode) THROW(fatal_error,"Four Momentum not conserved.");
  return Return_Value::New_Event;
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
	msg_Error()<<"ERROR in "<<METHOD<<": Too many bunch blobs required, "
		   <<"return 'Error' and hope for the best.\n";
	return Return_Value::Error;
      }
      flag=true;
    }
    else if ((*bit)->Has(blob_status::needs_beams) ||
	     (*bit)->Type()==btp::Elastic_Collision ||
	     (*bit)->Type()==btp::Soft_Diffractive_Collision ||
	     (*bit)->Type()==btp::Quasi_Elastic_Collision) {
      (*bit)->UnsetStatus(blob_status::needs_beams);
      for (size_t i=0;i<(*bit)->NInP();i++) {
	bunch = FillBunchBlob((*bit)->InParticle(i)->Beam(),(*bit)->InParticle(i));
	bloblist->push_front(bunch);
      }
    }
  }
  return (flag?Return_Value::Success:Return_Value::Nothing);
}

Blob * Beam_Remnant_Handler::FillBunchBlob(int beam,Particle * particle) 
{
  Blob *blob = new Blob();
  blob->SetType(btp::Bunch);
  blob->SetBeam(beam);
  blob->SetId();
  blob->SetStatus(blob_status::needs_beams &
		  blob_status::needs_reconnections &
		  blob_status::needs_softUE &
		  blob_status::needs_hadronization);
  blob->AddToOutParticles(particle);
  blob->SetPosition(p_remnants->GetRemnant(beam)->Position());
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

void Beam_Remnant_Handler::CleanUp(const size_t & mode)
{
  p_remnants->Reset();
}
