#ifndef MODEL_Main_FeynRules_Model_h
#define MODEL_Main_FeynRules_Model_h

#include "MODEL/Main/Model_Base.H"

namespace MODEL {
  class FeynRules_Model : public Model_Base {
  private :
    std::string m_partfile;
    void RunSpectrumGenerator();
  protected :
    void ParticleInit();
    void FillSpectrum();
  public :
    FeynRules_Model(std::string,std::string,bool);
    ~FeynRules_Model()                   {}
    
    bool CheckFlavours(int nin, int nout, ATOOLS::Flavour* flavs);
  };
}

#endif


