#ifndef MODEL_Main_Standard_Model_H
#define MODEL_Main_Standard_Model_H

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

    Standard_Model(std::string,std::string);
    bool ModelInit(const PDF::ISR_Handler_Map& isr);
    void InitVertices();

  };

}

#endif
