#ifndef MODEL_SM_Model_H
#define MODEL_SM_Model_H

#include "MODEL/Main/Model_Base.H"


namespace MODEL {

  class Standard_Model: public Model_Base {
  private:

    int  m_ckmorder, m_dec_g4;

    void FixEWParameters();
    void FixCKM();

    void ParticleInit();

    void InitQEDVertices();
    void InitQCDVertices();
    void InitEWVertices();

  public :

    Standard_Model();
    bool ModelInit();
    void InitVertices();

  };

}

#endif
