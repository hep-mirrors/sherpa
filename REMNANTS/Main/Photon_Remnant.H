#ifndef REMNANTS_Main_Photon_Remnant_H
#define REMNANTS_Main_Photon_Remnant_H

#include "PDF/Main/PDF_Base.H"
#include "REMNANTS/Main/Remnant_Base.H"
#include <map>

namespace REMNANTS {
class Photon_Remnant : public Remnant_Base {
private:
  PDF::PDF_Base *p_pdf;
  const ATOOLS::Flavour_Set *p_partons;
  ATOOLS::Flavour m_beamflav;

  ATOOLS::Particle *p_remnant, *p_recoiler, *p_spectator;

  double m_alpha, m_gamma, m_beta, m_invb, m_LambdaQCD;

  inline bool IsValence(ATOOLS::Particle *part) { return false; }
  void MakeSpectator(ATOOLS::Particle *parton);
  void CompensateColours();
  bool MakeRemnants();
  void MakeLongitudinalMomenta(ATOOLS::ParticleMomMap *ktmap,
                               const bool &copy = true);
  [[maybe_unused]] double SelectZ(const ATOOLS::Flavour &flav);
  ATOOLS::Particle *MakeParticle(const ATOOLS::Flavour &flav);
  ATOOLS::Flavour RemnantFlavour(const ATOOLS::Flavour &flav);
  void Output();

public:
  // constructor
  Photon_Remnant(PDF::PDF_Base *pdf, const unsigned int beam);

  // member functions
  bool FillBlob(ATOOLS::Blob *beamblob, ATOOLS::Particle_List *particlelist);
  bool FillBlob(ATOOLS::ParticleMomMap *ktmap, const bool &copy);
  bool AdjustKinematics();

  bool TestExtract(const ATOOLS::Flavour &flav, const ATOOLS::Vec4D &mom);
  void Reset(const bool &DIS = false);

  inline ATOOLS::Particle *GetRecoiler() { return p_recoiler; }
  inline ATOOLS::Particle *GetSpectator() { return p_spectator; }
}; // end of class Remnant_Base

} // end of namespace REMNANTS

#endif
