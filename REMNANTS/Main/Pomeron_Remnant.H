#ifndef REMNANTS_Main_Pomeron_Remnant_H
#define REMNANTS_Main_Pomeron_Remnant_H

#include "PDF/Main/PDF_Base.H"
#include "REMNANTS/Main/Remnant_Base.H"
#include <map>

namespace REMNANTS {
  class Pomeron_Remnant final : public Remnant_Base {
  private:
    PDF::PDF_Base*             p_pdf;
    const ATOOLS::Flavour_Set* p_partons;
    ATOOLS::Particle*          p_recoiler;
    const double               m_LambdaQCD;

    void   MakeRemnants(Colour_Generator* colours);
    bool   MakeLongitudinalMomenta(ATOOLS::ParticleMomMap* ktmap, const bool& copy);
    double SelectZ(const ATOOLS::Flavour& flav, double restmom, double remnant_masses) const;
    void   Output() const;

  public:
  // constructor
  Pomeron_Remnant(PDF::PDF_Base *pdf, const size_t& beam, const size_t& tag=0);

  // member functions
  bool FillBlob(Colour_Generator* colours, ATOOLS::ParticleMomMap* ktmap,
                const bool& copy) override;
  bool TestExtract(const ATOOLS::Flavour &flav, const ATOOLS::Vec4D &mom) override;
  void Reset(const bool &resc=false,const bool &DIS=false) override;

  inline ATOOLS::Particle *GetRecoiler() override { return p_recoiler; }
}; // end of class Remnant_Base

} // end of namespace REMNANTS

#endif
