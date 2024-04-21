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

  double m_LambdaQCD, m_beta_quark, m_beta_gluon;
  bool m_valence;

  ATOOLS::Particle * p_spectator, * p_recoiler;

  void MakeSpectator(ATOOLS::Particle *parton) override;
  void MakeRemnants();
  void MakeLongitudinalMomenta(ATOOLS::ParticleMomMap *ktmap, const bool &copy);
  double SelectZ(const ATOOLS::Flavour &flav, double restmom, double remnant_masses) const;
  ATOOLS::Particle *MakeParticle(const ATOOLS::Flavour &flav);
  double EstimateRequiredEnergy(bool needs_valence_quarks) const;
  void FindRecoiler();
  void Output() const;
public:
  // constructor
  Photon_Remnant(PDF::PDF_Base *pdf, const size_t & beam, const size_t & tag=0);

  // member functions
  bool TestExtract(const ATOOLS::Flavour &flav, const ATOOLS::Vec4D &mom) override;
  bool FillBlob(ATOOLS::ParticleMomMap *ktmap, const bool &copy) override;
  void CompensateColours();

  void Reset(const bool & resc=false, const bool &DIS=false) override;
  inline double MinRequiredEnergy(const ATOOLS::Flavour &flav) {
    return ( EstimateRequiredEnergy(!flav.IsQuark() && !m_valence) +
	     ATOOLS::Max(flav.HadMass(), m_LambdaQCD) );
  }
  
  inline ATOOLS::Particle *GetRecoiler() override { return p_recoiler; }
  inline ATOOLS::Particle *GetSpectator() override { return p_spectator; }
}; // end of class Remnant_Base

} // end of namespace REMNANTS

#endif
