#ifndef REMNANTS_Main_Proton_Remnant_H
#define REMNANTS_Main_Proton_Remnant_H

#include "REMNANTS/Main/Remnant_Base.H"
#include "REMNANTS/Tools/Form_Factor.H"
#include "PDF/Main/PDF_Base.H"
#include <map>

namespace REMNANTS {
  class Hadron_Remnant: public Remnant_Base {
  private:
    PDF::PDF_Base             * p_pdf;
    const ATOOLS::Flavour_Set * p_partons;
    FlavourList                 m_constituents;
    
    ATOOLS::Particle * p_valence, * p_remnant, * p_recoiler, * p_spectator;

    bool   m_valence;
    double m_alpha, m_gamma, m_beta, m_invb, m_LambdaQCD;

    void ConstructConstituentFlavours();

    bool IsValence(ATOOLS::Particle * part);
    void MakeSpectator(ATOOLS::Particle * parton);
    void CompensateColours();
    void SquashFlavourSinglets();
    void SquashColourSinglets();
    bool ReplaceInBlob(ATOOLS::Blob * shower,const bool & check,
		       const size_t & pos,size_t col[2]);
    bool MakeRemnants();
    void MakeLongitudinalMomenta(ATOOLS::ParticleMomMap *ktmap,
				 const bool & copy=true);
    double SelectZ(const ATOOLS::Flavour &flav, double restmom, double remnant_masses) const;
    ATOOLS::Particle * MakeParticle(const ATOOLS::Flavour & flav);
    ATOOLS::Flavour    RemnantFlavour(const ATOOLS::Flavour & flav);
    double EstimateRequiredEnergy() const;
    void Output();
  public:
    Hadron_Remnant(PDF::PDF_Base * pdf,const unsigned int & beam,
		   const unsigned int & tag=0);    
    ~Hadron_Remnant() = default;

    bool FillBlob(ATOOLS::ParticleMomMap *ktmap,const bool & copy=true) final;
    bool TestExtract(const ATOOLS::Flavour &flav,const ATOOLS::Vec4D &mom) final;
    void Reset(const bool & resc=false,const bool & DIS=false) final;
    
    inline double MinRequiredEnergy(const double energy, const ATOOLS::Flavour &flav) {
      return ( EstimateRequiredEnergy()+energy+ATOOLS::Max(flav.HadMass(), m_LambdaQCD) );
    }
    inline ATOOLS::Particle * GetRecoiler()            { return p_recoiler; }
    inline ATOOLS::Particle * GetSpectator()           { return p_spectator; }
  };
}

#endif
