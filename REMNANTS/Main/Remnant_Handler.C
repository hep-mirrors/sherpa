#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Math/Random.H"
#include "BEAM/Main/Beam_Spectra_Handler.H"
#include "PDF/Main/ISR_Handler.H"
#include "REMNANTS/Main/Electron_Remnant.H"
#include "REMNANTS/Main/Hadron_Remnant.H"
#include "REMNANTS/Main/No_Remnant.H"
#include "REMNANTS/Main/Photon_Remnant.H"
#include "REMNANTS/Main/Remnant_Handler.H"


using namespace REMNANTS;
using namespace ATOOLS;
using namespace std;

Remnant_Handler::Remnant_Handler(PDF::ISR_Handler *isr,
                                 BEAM::Beam_Spectra_Handler *beam)
    : p_softblob(nullptr), m_check(true), m_output(true) {
  InitializeRemnants(isr, beam);
  DefineRemnantStrategy();
  InitializeKinematicsAndColours();
}

Remnant_Handler::~Remnant_Handler() {
  for (size_t i(0); i < 2; ++i) {
    if (p_remnants[i] != nullptr)
      delete p_remnants[i];
  }
}

void Remnant_Handler::InitializeRemnants(PDF::ISR_Handler *isr,
                                         BEAM::Beam_Spectra_Handler *beam) {
  for (size_t i = 0; i < 2; ++i) {
    p_remnants[i] = nullptr;
    Flavour flav = isr->Flav(i);
    if (isr->PDF(i) != nullptr) {
      if (flav.IsHadron())
        p_remnants[i] = new Hadron_Remnant(isr->PDF(i), i);
      else if (flav.IsLepton())
        p_remnants[i] = new Electron_Remnant(isr->PDF(i), i);
      else if (flav.IsPhoton()) {
        p_remnants[i] = new Photon_Remnant(isr->PDF(i), i);
      }
    }
    if (p_remnants[i] == nullptr)
      p_remnants[i] = new No_Remnant(i);
  }
  // Finish the initialisation of the Remnant_Bases: make sure they know
  // each other, their beam, the Colour_Generator, and hand them also to
  // the ISR_Handler.
  // TODO: this latter part may become obsolete - I will have to check this.
  for (size_t i = 0; i < 2; ++i) {
    p_remnants[i]->SetPartner(p_remnants[1 - i]);
    p_remnants[i]->SetBeam(beam->GetBeam(i));
    p_remnants[i]->SetColours(&m_colours);
    p_remnants[i]->Reset();
    isr->SetRemnant(p_remnants[i], i);
  }
}

void Remnant_Handler::DefineRemnantStrategy() {
  // Pretty self-explanatory.  This covers the way of how we make sure that
  // potential intrinsic transverse momenta in the beam breakups do not lead
  // to violation of four-momentum conservation.  We have pretty much 4
  // options:
  // - simple, where the remnants do not break up;
  // - ll where the breakup is collinear only and we only need to boost
  //   the system in the end;
  // - DIS, where one breakup is collinear and the other involves transverse
  //   momenta; and
  //   (TODO: this is the one I still need to implement)
  // - hh, where both breakups generate transverse momenta.
  // For the latter two we will realize four-momentum conservation through
  // insertion of a "soft" blob, mainly a garbage collection where we collect
  // particles and shuffle them in a pretty minimal fashion.
  if (p_remnants[0]->Type() == rtp::intact &&
      p_remnants[1]->Type() == rtp::intact)
    m_type = strat::simple;
  else if (p_remnants[0]->Type() == rtp::lepton &&
           p_remnants[1]->Type() == rtp::lepton)
    m_type = strat::ll;
  else if ((p_remnants[0]->Type() == rtp::hadron ||
            p_remnants[0]->Type() == rtp::photon) &&
           (p_remnants[1]->Type() == rtp::lepton ||
            p_remnants[1]->Type() == rtp::intact))
    m_type = strat::DIS1;
  else if ((p_remnants[0]->Type() == rtp::lepton ||
            p_remnants[0]->Type() == rtp::intact) &&
           (p_remnants[1]->Type() == rtp::hadron ||
            p_remnants[1]->Type() == rtp::photon))
    m_type = strat::DIS2;
  else if ((p_remnants[0]->Type() == rtp::hadron ||
            p_remnants[0]->Type() == rtp::photon) &&
           (p_remnants[1]->Type() == rtp::hadron ||
            p_remnants[1]->Type() == rtp::photon))
    m_type = strat::hh;
  else if ((p_remnants[0]->Type() == rtp::lepton &&
            p_remnants[1]->Type() == rtp::intact) ||
           (p_remnants[0]->Type() == rtp::intact &&
            p_remnants[1]->Type() == rtp::lepton))
    m_type = strat::simple;
  else {
    msg_Error() << METHOD << " throws error: no strategy found for remnants "
                << p_remnants[0]->Type() << " & " << p_remnants[1]->Type()
                << "\n";
    exit(1);
  }
}

void Remnant_Handler::InitializeKinematicsAndColours() {
  m_kinematics.Initialize(this);
  m_colours.Initialize(this);
  m_decorrelator.Initialize(this);
}

bool Remnant_Handler::ExtractShowerInitiators(Blob *const showerblob) {
  // This method is called after each successful parton shower off a hard
  // scatter.  It extracts the two initial particles from the shower and
  // extracts them from the corresponding beam remnant.

  // Make sure only shower blobs with exactly two initiators are treated,
  // and only once.
  // (Shower blob with more initiators are typically from hadron decays.)
  if (showerblob->Type() != btp::Shower ||
      m_treatedshowerblobs.find(showerblob) != m_treatedshowerblobs.end())
    return true;
  size_t countIn = 0;
  for (size_t i = 0; i < showerblob->NInP(); ++i) {
    if (!showerblob->InParticle(i)->ProductionBlob()) countIn++;
  }
  if (countIn!=2) return true;
  // Now extract the shower initiators from the remnants - they will get added
  // to the lists of extracted particles for each remnant and their colour will
  // be added to the Colour_Generator in each beam.
  for (size_t i = 0; i < showerblob->NInP(); ++i) {
    Particle *part = showerblob->InParticle(i);
    if (part->ProductionBlob() != nullptr) continue;
    // Make sure extraction works out - mainly subject to energy conservation
    if (!Extract(part, part->Beam()))      return false;
  }
  m_treatedshowerblobs.insert(showerblob);
  return true;
}

void Remnant_Handler::ConnectColours(ATOOLS::Blob *const showerblob) {
  // After each showering step, we try to compensate some of the colours.
  // In the absence of multiple parton interactions this will not involve
  // anything complicated with colours.  In each step, the shower initiators
  // are checked for being a valence quark.  In case they are, remnants
  // (equivalent to recoilers) are generated, usually diquarks for protons,
  // if they are seaquarks, suitable spectators are generated.  The colours
  // of the shower initiators and, possibly, spectators will be added to a
  // stack which will in turn partially replace the new colours.  This is
  // handled in the Colour_Generator.
  m_colours.ConnectColours(showerblob);
}

Return_Value::code
Remnant_Handler::MakeBeamBlobs(Blob_List *const bloblist,
                               Particle_List *const particlelist) {
  // Adding the blobs related to the breakup of incident beams: one for each
  // beam, plus, potentially a third one to balance transverse momenta.
  InitBeamAndSoftBlobs(bloblist);
  // Fill in the transverse momenta through the Kinematics_Generator.
  // Check for colour connected parton-pairs including beam partons and
  // add soft gluons in between them if their invariant mass is too large.
  // This still needs debugging - therefore it is commented out.
  if (!m_kinematics.FillBlobs(bloblist) || !CheckBeamBreakup(bloblist)) {
    // || !m_decorrelator(p_softblob)) {
    Reset();
    msg_Error()
        << "Warning in " << METHOD
        << ": FillBlobs or CheckBeamBreakup failed. Will return new event\n";
    return Return_Value::New_Event;
  }
  Reset();
  return Return_Value::Success;
}

void Remnant_Handler::InitBeamAndSoftBlobs(Blob_List *const bloblist) {
  // Making a new blob (softblob) to locally compensate 4 momentum.
  // Ultimately, it will reflect different strategies of how to compensate
  // intrinsic kperp: hadron colliders vs. DIS
  // For hh collisions, the softblob is inserted after the beam blobs and
  // before the parton shower blobs, to capture the kperp of both beams,
  // pair-by-pair of the shower initiators, to combine them into one
  // transverse momentum, "kicking" the full shower blob around.  In this
  // way we can guarantee that the kperps between the two beams and their
  // recoils are compensated by each other.
  // For DIS we do not want to change the kinematics of the incident lepton,
  // so we cannot use the other beam for the recoils.  Instad we take the
  // full coloured final state after the shower for the kperp compensation.
  // Effectively we distribute the kperp of the shower initiator over all
  // coloured FS particles and then shuffle their momenta together with the
  // beam remnants.  To visualise this better, here the soft blob is
  // inserted after both beam and shower blobs.
  if (!(m_type == strat::simple || m_type == strat::ll)) {
    p_softblob = m_kinematics.MakeSoftBlob();
    if (m_type == strat::DIS1 || m_type == strat::DIS2)
      bloblist->push_back(p_softblob);
    else
      bloblist->push_front(p_softblob);
  }
  // Look for shower blobs that need beams and unset the flag
  for (auto &bit : *bloblist) {
    if (bit->Has(blob_status::needs_beams) && bit->Type() == btp::Shower) {
      bit->UnsetStatus(blob_status::needs_beams);
    }
  }
  // Remnant bases will generate their beam blobs, reset the incoming
  // four-momenta and make the two beam blobs
  m_colours.ResetFlags();
  for (size_t beam = 0; beam < 2; beam++) {
    bloblist->push_front(p_remnants[beam]->MakeBlob());
  }
}

bool Remnant_Handler::CheckBeamBreakup(Blob_List *bloblist) {
  // Final checks on beam breakup: four-momentum and colour conservation
  if (m_type == strat::simple || !m_check)
    return true;
  bool ok = true;
  for (size_t beam = 0; beam < 2; beam++) {
    if (!p_remnants[beam]->GetBlob()->MomentumConserved() ||
        !p_remnants[beam]->GetBlob()->CheckColour()) {
      ok = false;
      if (m_output)
        msg_Error() << "Error in " << METHOD << ": "
                    << "colour or four-momentum not conserved in:\n"
                    << (*p_remnants[beam]->GetBlob()) << "\n";
    }
  }
  if (!p_softblob)
    return ok;
  if (!p_softblob->MomentumConserved() || !p_softblob->CheckColour()) {
    ok = false;
    if (m_output)
      msg_Error() << "Error in " << METHOD << ": "
                  << "colour or four-momentum not conserved in:\n"
                  << (*p_softblob) << "\n";
  }
  return ok;
}

void Remnant_Handler::SetImpactParameter(const double & b) {
  Vec4D  pos = b/2.*Vec4D(0.,1.,0.,0.);
  for (size_t i=0;i<2;i++) p_remnants[i]->SetPosition((i==0?1.:-1.) * pos);      
}

bool Remnant_Handler::Extract(ATOOLS::Particle * part,const unsigned int beam) {
  // Extracting a particle from a remnant only works for positive energies.
  if (part->Momentum()[0] < 0.) {
    msg_Error() << METHOD << " yields shower with negative incoming energies.\n"
                << (*part->DecayBlob()) << "\n";
    return false;
  }
  return p_remnants[beam]->Extract(part);
}

void Remnant_Handler::Reset() {
  bool DIS = m_type == strat::DIS1 || m_type == strat::DIS2;
  for (size_t beam = 0; beam < 2; beam++)
    p_remnants[beam]->Reset(DIS);
  m_treatedshowerblobs.clear();
  m_kinematics.Reset();
  m_colours.Reset();
}
