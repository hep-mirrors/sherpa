#ifndef MODEL__UFO__UFO_Model_H
#define MODEL__UFO__UFO_Model_H

#include "MODEL/Main/Model_Base.H"
#include "MODEL/UFO/UFO_Param_Reader.H"

namespace UFO{

  class UFO_Model : public MODEL::Model_Base
  {
  public:
    
    UFO_Model(std::string path, std::string file, bool elementary);
    ~UFO_Model();
    bool ModelInit(const PDF::ISR_Handler_Map& isr);
    void SetSMMass(const kf_code &kf,const double &m);
    void SetSMMasses();

  protected:

    UFO::UFO_Param_Reader* p_dataread;

    const Complex complexconjugate(const Complex& arg);
    const Complex re(const Complex& arg);
    const Complex im(const Complex& arg);
    const Complex complex(double real, double imag);

    virtual void ParamInit()    {}
    virtual void ParticleInit() {}

  };

}

#endif
