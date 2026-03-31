#ifndef MODEL_SQED_Model_H
#define MODEL_SQED_Model_H

#include "MODEL/Main/Model_Base.H"


namespace MODEL {

  class SQED: public Model_Base {
  private:

    int  m_ckmorder, m_dec_g4;

    void FixEWParameters();
    void FixCKM();

    void ParticleInit();

    void InitQEDVertices();
    void InitQCDVertices();
    void InitEWVertices();

  public :

    SQED();
    bool ModelInit();
    void InitVertices();
    void addPhotonScalarVertex(const ATOOLS::Flavour &a, const ATOOLS::Flavour &b, const ATOOLS::Kabbala &cpl);

  };

}

#endif
