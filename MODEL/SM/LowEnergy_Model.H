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
    void InitEWVertices();
  public :
    LowEnergy_Model();
    bool ModelInit();
    void InitVertices();
  };

}

#endif
