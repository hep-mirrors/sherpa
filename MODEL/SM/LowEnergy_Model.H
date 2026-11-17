#ifndef MODEL_LowEnergy_Model_H
#define MODEL_LowEnergy_Model_H

#include "MODEL/Main/Model_Base.H"
#include <list>

namespace MODEL {
  class LowEnergy_Model: public Model_Base {
  private:
    std::list<ATOOLS::Flavour> m_Diracs, m_PseudoScalars;
    double m_alpha, m_sinthetaW;
    bool   m_threepion_contact;

    void ParticleInit();

    void InitQEDConstants();

    void InitQEDVertices();
    void InitThreePionVertices(const ATOOLS::Kabbala & cpl);
    //! THREE_PION_FORM_FACTOR: RhoOmegaMixing, the rho-omega mixing chain
    void InitRhoOmegaMixingVertices(const ATOOLS::Kabbala & cpl_grp);
    void InitTwoPhotonVertices();
    //! FOUR_PION_FORM_FACTOR: Active, gamma* -> pi+ pi- 2pi0 and 2pi+ 2pi-
    void InitFourPionVertices(const ATOOLS::Kabbala & cpl);
    void InitEWVertices();
  public :
    //! g_{rho pi pi} from the rho width; THREE_PION_FORM_FACTOR: key overrides it
    static double RhoPiPiCoupling(const ATOOLS::Flavour & rho,const std::string & key);
    //! g_{gamma rho pi} from the WZW anomaly; THREE_PION_FORM_FACTOR: g_gammarhopi overrides it
    static double GammaRhoPiCoupling();
    LowEnergy_Model();
    bool ModelInit();
    void InitVertices();
  };

}

#endif
