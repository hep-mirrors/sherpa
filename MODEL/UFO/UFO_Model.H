#ifndef MODEL__UFO__UFO_Model_H
#define MODEL__UFO__UFO_Model_H

#include "MODEL/Main/Model_Base.H"
#include "MODEL/UFO/UFO_Param_Reader.H"
#include "ATOOLS/Math/Kabbala.H"
#include "MODEL/UFO/Variations.H"

namespace UFO{

  typedef std::map<std::string,std::string> StringMap;

  class UFO_Model : public MODEL::Model_Base
  {

  public:

    UFO_Model(bool elementary);
    ~UFO_Model();
    bool ModelInit();
    void SetSMMass(const kf_code &kf,const double &m);
    void SetSMMasses();
    void SetMassiveFlags();
    void SetStableFlags();
    Variations* GetParameterVariations() override {return p_variations;}
    bool InitVariations() override {if (!p_variations) p_variations = new Variations(); return true;}

  protected:

    StringMap m_lorentz_map;
    std::map<std::string, ATOOLS::Kabbala> m_cpls;
    UFO::UFO_Param_Reader* p_dataread;
    Variations* p_variations = nullptr;

    /*
    static Complex complexconjugate(const Complex& arg);
    static Complex re(const Complex& arg);
    static Complex im(const Complex& arg);
    static Complex complex(double real, double imag);
    static Complex sqrt(const double& arg);
    static Complex sqrt(const Complex& arg);
    static double  ToDouble(const Complex& arg);
    static inline double ToDouble(const ATOOLS::Kabbala& k) {return ATOOLS::ToDouble(k);}
    static inline ATOOLS::Kabbala sqrt(const ATOOLS::Kabbala& k) {return ATOOLS::sqrt(k);}
    static inline ATOOLS::Kabbala complexconjugate(const Complex& c) {return ATOOLS::}*/

    virtual void ParamInit()    {}
    virtual void ParticleInit() {}

    virtual void FillLorentzMap() = 0;
    virtual std::string MappedLorentzName(const std::string& label) const;

  };

}

#endif
