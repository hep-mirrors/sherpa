#include "REMNANTS/Tools/Frame_Guard.H"

#include "REMNANTS/Main/Remnant_Base.H"
#include "ATOOLS/Phys/Blob_List.H"

using namespace REMNANTS;
using namespace ATOOLS;

Frame_Guard::Frame_Guard(Blob_List* bloblist,
                         const std::array<std::shared_ptr<Remnant_Base>, 2>& remnants,
                         const bool& kinematics)
  : p_bloblist(bloblist), p_remnants(remnants), m_active(false)
{
  if (!kinematics) return;
  if (p_remnants[0]->GetBeam()->Type() == BEAM::beamspectrum::monochromatic &&
      p_remnants[1]->GetBeam()->Type() == BEAM::beamspectrum::monochromatic)
    return;
  const Vec4D ptot =
    p_remnants[0]->IncomingMomentum() + p_remnants[1]->IncomingMomentum();
  // Identity short-circuit, and no transverse pair momenta - the current
  // spectra construct their out-momenta exactly longitudinally.
  if (Vec3D(ptot).Sqr() < 1.e-24*sqr(ptot[0])) return;
  if (ptot.PPerp2() > 1.e-12*sqr(ptot[0])) return;
  m_cms = Poincare(ptot);
  std::set<Particle*> treateds;
  p_bloblist->Boost(m_cms, &treateds);
  for (auto& remnant : p_remnants) remnant->SetWindowFrame(&m_cms);
  m_active = true;
}

Frame_Guard::~Frame_Guard()
{
  if (!m_active) return;
  for (auto& remnant : p_remnants) remnant->SetWindowFrame(nullptr);
  Poincare lab(m_cms);
  lab.Invert();
  std::set<Particle*> treateds;
  p_bloblist->Boost(lab, &treateds);
}
