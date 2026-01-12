#ifndef UFO_VARIATIONS_TEST_MODEL_SM
#define UFO_VARIATIONS_TEST_MODEL_SM

#include "MODEL/UFO/UFO_Model.H"
#include "ATOOLS/Math/Kabbala.H"


namespace MODEL {

  class Test_Standard_Model: public Model_Base {
  private:

    UFO::Variations* p_variations = nullptr;
    int  m_ckmorder, m_dec_g4;
    std::map<std::string, ATOOLS::Kabbala> k_constants;
    ATOOLS::Kabbala Kabbala_Constant(std::string name) {if (k_constants.count(name)!=1) return ATOOLS::Kabbala(); return k_constants[name];}
    bool CreateKabbalaFromConst(std::string var_name) {return k_constants.insert(std::make_pair(var_name, ATOOLS::Kabbala(var_name, ScalarConstant(var_name)))).second;};

    void FixEWParameters();
    void FixCKM();

    void ParticleInit();

    void InitQEDVertices();
    void InitQCDVertices();
    void InitEWVertices();

  public :

    Test_Standard_Model();
    bool ModelInit();
    void InitVertices();
    
    UFO::Variations* GetParameterVariations() override {return p_variations;}
    bool InitVariations() override {if (!p_variations) p_variations = new UFO::Variations(); return true;}
  };

}

#endif
