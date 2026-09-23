#ifndef MODEL_LowEnergy_Model_H
#define MODEL_LowEnergy_Model_H

#include "MODEL/Main/Model_Base.H"


namespace MODEL {
  class LowEnergy_Model: public Model_Base {
  private:
    double m_alpha, m_sinthetaW;
    
    void ParticleInit();

    void InitQEDVertices();
    void InitEWVertices();
  public :
    LowEnergy_Model();
    bool ModelInit();
    void InitVertices();
  };

}

#endif
