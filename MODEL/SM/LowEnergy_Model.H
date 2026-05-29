#ifndef MODEL_LowEnergy_Model_H
#define MODEL_LowEnergy_Model_H

#include "MODEL/Main/Model_Base.H"
#include <list>

namespace MODEL {
  class LowEnergy_Model: public Model_Base {
  private:
    std::list<ATOOLS::Flavour> m_Diracs, m_PseudoScalars;
    double m_alpha, m_sinthetaW;
    
    void ParticleInit();

    void InitQEDConstants();
    
    void InitQEDVertices();
    void InitEWVertices();
  public :
    LowEnergy_Model();
    bool ModelInit();
    void InitVertices();
  };

}

#endif
