#include "SHERPA/SoftPhysics/Beam_Remnant_Handler.H"
#include "BEAM/Main/Beam_Spectra_Handler.H"
#include "PDF/Main/ISR_Handler.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Org/Exception.H"

using namespace SHERPA;
using namespace ATOOLS;

Beam_Remnant_Handler::Beam_Remnant_Handler(
        BEAM::Beam_Spectra_Handler* const beam,
        REMNANTS::Remnant_Handler* const  remnants,
        Soft_Collision_Handler* const     softcollisions)
    : p_remnants(remnants), p_bunchremnants(nullptr),
      p_shrimps(softcollisions ? softcollisions->GetShrimps() : nullptr),
      p_beam(beam), m_bunchrescatter(false), m_fill(true), m_vmode(false)
{
  Settings& s = Settings::GetMainSettings();
  m_fill  = s["BEAM_REMNANTS"].Get<bool>();
  m_vmode = s["BRH_VMODE"].SetDefault(false).Get<bool>();
  p_remnants->SetScale2(sqr(4.0));
  m_name = std::string("Parametrised");
}

Beam_Remnant_Handler::~Beam_Remnant_Handler() = default;

Return_Value::code Beam_Remnant_Handler::FillBeamAndBunchBlobs(
        Blob_List* const bloblist, const bool& onlyBunch)
{
  if (!m_fill) return TreatNoFill(bloblist);
  Return_Value::code fbc(Return_Value::Nothing);
  for (auto* bit : *bloblist) {
    if (bit->Type() == btp::Beam) return fbc;
  }
  if (!onlyBunch) {
    if (p_shrimps) fbc = p_shrimps->MakeBeamBlobs(bloblist);
    else           fbc = p_remnants->MakeBeamBlobs(bloblist);
    if (fbc==Return_Value::New_Event && m_vmode)
      THROW(fatal_error,"Four Momentum not conserved.");
    if (fbc!=Return_Value::Success) return fbc;
  }
  return FillBunchBlobs(bloblist);
}

Return_Value::code
Beam_Remnant_Handler::FillRescatterBeamBlobs(Blob_List *const bloblist) {
  Return_Value::code fbc = p_bunchremnants->MakeBeamBlobs(bloblist, true);
  return fbc;
}

Return_Value::code
Beam_Remnant_Handler::TreatNoFill(Blob_List *const bloblist)
{
  bool set(false);
  for (auto* bit : *bloblist) {
    if (bit->Has(blob_status::needs_beams)) {
      bit->UnsetStatus(blob_status::needs_beams);
      bit->UnsetStatus(blob_status::internal_flag);
      set=true;
    }
  }
  if (!set) return Return_Value::Nothing;
  if (bloblist->FourMomentumConservation()) return Return_Value::Success;
  msg_Tracking()<<METHOD<<" found four momentum conservation error.\n";
  if (m_vmode) THROW(fatal_error,"Four Momentum not conserved.");
  return Return_Value::New_Event;
}

Return_Value::code Beam_Remnant_Handler::FillBunchBlobs(
        Blob_List* const bloblist)
{
  for (auto* bit : *bloblist) {
    if (bit->Type() == btp::Bunch) return Return_Value::Nothing;
  }
  p_beam->FixPositions();
  if (!m_bunchrescatter) return FillSimpleBunchBlobs(bloblist);
  if (FillRescatterBunchBlobs(bloblist))
    return p_schandler->GenerateBunchRescatter(bloblist);
  return Return_Value::Nothing;
}

bool Beam_Remnant_Handler::FillRescatterBunchBlobs(
        ATOOLS::Blob_List* const bloblist)
{
  bool flag(false);
  m_beam = 0;
  for (auto* bit : *bloblist) {
    if (bit->Has(blob_status::needs_beams) &&
        (bit->Type() == btp::Beam || bit->Type() == btp::Shower)) {
      bit->UnsetStatus(blob_status::needs_beams);
      if (bit->NInP() == 1 && bit->NOutP() == 1 &&
          !bit->InParticle(0)->Flav().IsHadron() &&
          bit->InParticle(0)->Flav() == bit->OutParticle(0)->Flav())
        bit->UnsetStatus(blob_status::needs_softUE);
      Blob* bunch = FillBunchBlob(bit->Beam(), bit->InParticle(0));
      bunch->AddStatus(blob_status::needs_beamRescatter);
      p_schandler->SetPosition(m_beam-1,bunch->Position());
      bloblist->push_front(bunch);
      if (m_beam>2) THROW(fatal_error,"Too many bunch blobs required");
      flag=true;
    }
  }
  return flag;
}

ATOOLS::Return_Value::code Beam_Remnant_Handler::FillSimpleBunchBlobs(
        ATOOLS::Blob_List* const bloblist)
{
  ATOOLS::Return_Value::code flag(Return_Value::Nothing);
  m_beam = 0;
  for (auto* bit : *bloblist) {
    if (bit->Has(blob_status::needs_beams) &&
        (bit->Type() == btp::Beam || bit->Type() == btp::Shower)) {
      bit->UnsetStatus(blob_status::needs_beams);
      bloblist->push_front(FillBunchBlob(bit->Beam(), bit->InParticle(0)));
      if (m_beam>2) THROW(fatal_error,"Too many bunch blobs required");
      flag = Return_Value::Success;
    } else if (bit->Has(blob_status::needs_beams) ||
               bit->Type() == btp::Elastic_Collision ||
               bit->Type() == btp::Soft_Diffractive_Collision ||
               bit->Type() == btp::Quasi_Elastic_Collision) {
      bit->UnsetStatus(blob_status::needs_beams);
      for (auto* part : *bit->InParticles())
        bloblist->push_front(FillBunchBlob(part->Beam(), part));
      flag = Return_Value::Success;
    }
  }
  return flag;
}


Blob * Beam_Remnant_Handler::FillBunchBlob(int beam,Particle * particle)
{
  Blob *blob = new Blob();
  blob->SetType(btp::Bunch);
  blob->SetBeam(beam);
  blob->SetId();
  blob->AddToOutParticles(particle);
  if (particle->Flav()==p_beam->GetBeam(beam)->Beam() &&
      IsEqual(particle->E(),p_beam->GetBeam(beam)->InMomentum()[0],1.e-6)) {
    Particle *p = new Particle(*particle);
    p->SetNumber(0);
    blob->AddToInParticles(p);
    blob->SetPosition(p_remnants->GetRemnant(beam)->Position());
  }
  else {
    Particle* p = new Particle(-1, p_beam->GetBeam(beam)->Beam(),
                               p_beam->GetBeam(beam)->InMomentum());
    p->SetNumber(0);
    p->SetStatus(part_status::decayed);
    p->SetFinalMass();
    blob->AddToInParticles(p);
    p = new Particle(-1, p_beam->GetBeam(beam)->Remnant(),
                     p_beam->GetBeam(beam)->InMomentum() -
                             particle->Momentum());
    p->SetNumber(0);
    p->SetStatus(part_status::decayed);
    p->SetFinalMass();
    blob->AddToOutParticles(p);
    blob->SetPosition(p_beam->GetBeam(beam)->Position());
  }
  m_beam++;
  return blob;
}

void Beam_Remnant_Handler::CleanUp(const size_t & mode)
{
  p_remnants->Reset();
}
