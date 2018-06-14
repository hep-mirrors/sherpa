#ifndef REMNANTS_Main_Proton_Remnant_H
#define REMNANTS_Main_Proton_Remnant_H

#include "REMNANTS/Main/Remnant_Base.H"
#include "PDF/Main/PDF_Base.H"
#include <map>

namespace REMNANTS {
  class Hadron_Remnant: public Remnant_Base {
  public:
  private:
    PDF::PDF_Base             * p_pdf;
    const ATOOLS::Flavour_Set * p_partons;
    ATOOLS::Flavour             m_beamflav;
    FlavourList                 m_constituents;
    
    ATOOLS::Particle * p_valence, * p_remnant, * p_recoiler, * p_spectator;

    bool   m_valence;
    double m_alpha, m_gamma, m_beta, m_invb, m_LambdaQCD;

    void ConstructConstituentFlavours();
    void ConstructAllowedPartons();

    bool IsValence(ATOOLS::Particle * part);
    void MakeSpectator(ATOOLS::Particle * parton);
    bool MakeRemnants();
    void MakeLongitudinalMomenta(ATOOLS::ParticleMomMap *ktmap,const bool & copy=true);
    bool ConstructPartner(ATOOLS::Particle * part);
    double SelectZ(const ATOOLS::Flavour & flav,const bool & isvalence);
    ATOOLS::Particle * MakeParticle(const ATOOLS::Flavour & flav);
    ATOOLS::Flavour    RemnantFlavour(const ATOOLS::Flavour & flav);
    void Output();
  public:
    Hadron_Remnant(PDF::PDF_Base * pdf,const unsigned int beam);    
    ~Hadron_Remnant() {}

    bool FillBlob(ATOOLS::ParticleMomMap *ktmap,const bool & copy=true);
    bool AdjustKinematics();
    bool AdjustColours();
    bool TestExtract(const ATOOLS::Flavour &flav,const ATOOLS::Vec4D &mom);
    void Reset();
    
    inline ATOOLS::Particle * GetRecoiler()  { return p_recoiler; }
    inline ATOOLS::Particle * GetSpectator() { return p_spectator; }
  };
}

#endif
