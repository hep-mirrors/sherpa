#ifndef REMNANTS_Main_Photon_Remnant_H
#define REMNANTS_Main_Photon_Remnant_H

#include "PDF/Main/PDF_Base.H"
#include "REMNANTS/Main/Remnant_Base.H"
#include "REMNANTS/Tools/Form_Factor.H"
#include <map>

namespace REMNANTS {
class Photon_Remnant : public Remnant_Base {
private:
  PDF::PDF_Base *p_pdf;
  const ATOOLS::Flavour_Set *p_partons;
  ATOOLS::Flavour m_beamflav;

  double m_LambdaQCD, m_beta_quark, m_beta_gluon, m_remnant_masses;
  bool m_valence;

  ATOOLS::Particle * p_spectator, * p_recoiler;
  Form_Factor        m_ff;

  void MakeSpectator(ATOOLS::Particle *parton) override;
  void MakeRemnants();
  void MakeLongitudinalMomenta(ATOOLS::ParticleMomMap *ktmap, const bool &copy);
  double SelectZ(const ATOOLS::Flavour &flav);
  ATOOLS::Particle *MakeParticle(const ATOOLS::Flavour &flav);
  double EstimateRequiredEnergy(bool needs_valence_quarks);
  void FindRecoiler();
  void Output();
public:
  // constructor
  Photon_Remnant(PDF::PDF_Base *pdf, const unsigned int beam);

  // member functions
  bool FillBlob(ATOOLS::ParticleMomMap *ktmap, const bool &copy) override;
  bool TestExtract(const ATOOLS::Flavour &flav, const ATOOLS::Vec4D &mom) override;
  void Reset(const bool &DIS) override;
  void CompensateColours();

  inline ATOOLS::Particle *GetRecoiler() override { return p_recoiler; }
  inline ATOOLS::Particle *GetSpectator() override { return p_spectator; }
}; // end of class Remnant_Base

} // end of namespace REMNANTS

#endif
