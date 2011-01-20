#ifndef MODEL_Main_FeynRules_Model_h
#define MODEL_Main_FeynRules_Model_h

#include "MODEL/Main/Model_Base.H"

namespace MODEL {
  class FeynRules_Model : public Model_Base {
  private :
    std::string m_partfile;
    void RunSpectrumGenerator(PDF::ISR_Handler *const isr);
  protected :
    void ParticleInit();
    void FillSpectrum(PDF::ISR_Handler *const isr);
  public :
    FeynRules_Model(std::string,std::string,bool);
    ~FeynRules_Model()                   {}
    
    bool CheckFlavours(int nin, int nout, ATOOLS::Flavour* flavs);
    bool ModelInit(PDF::ISR_Handler *const isr);
  };
}

#endif


