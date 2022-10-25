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

  ATOOLS::Particle *p_remnant = nullptr;
  ATOOLS::Particle *p_recoiler = nullptr;

  void MakeSpectator(ATOOLS::Particle *parton);
  void CompensateColours();
  bool MakeRemnants();
  void MakeLongitudinalMomenta(ATOOLS::ParticleMomMap *ktmap, const bool &copy);
  double SelectZ(const ATOOLS::Flavour &flav);
  ATOOLS::Particle *MakeParticle(const ATOOLS::Flavour &flav);
  void Output();

public:
  // constructor
  Photon_Remnant(PDF::PDF_Base *pdf, const unsigned int beam);

  // member functions
  bool FillBlob(ATOOLS::ParticleMomMap *ktmap, const bool &copy);

  bool TestExtract(const ATOOLS::Flavour &flav, const ATOOLS::Vec4D &mom);
  void Reset(const bool &DIS);

  inline ATOOLS::Particle *GetRecoiler() { return p_recoiler; }
}; // end of class Remnant_Base

} // end of namespace REMNANTS

#endif
