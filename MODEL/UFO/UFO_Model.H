#ifndef MODEL__UFO__UFO_Model_H
#define MODEL__UFO__UFO_Model_H

#include "MODEL/Main/Model_Base.H"
#include "MODEL/UFO/UFO_Param_Reader.H"
#include "ATOOLS/Math/Kabbala.H"

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

  protected:

    StringMap m_lorentz_map;
    std::map<std::string, ATOOLS::Kabbala> m_cpls;
    UFO::UFO_Param_Reader* p_dataread;

    virtual void ParamInit()    {}
    virtual void ParticleInit() {}

    virtual void FillLorentzMap() = 0;
    virtual std::string MappedLorentzName(const std::string& label) const;

  };

}

#endif
