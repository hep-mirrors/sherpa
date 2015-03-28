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

  protected:

    UFO::UFO_Param_Reader* p_dataread;

    const Complex complexconjugate(const Complex& arg);
    const Complex re(const Complex& arg);
    const Complex im(const Complex& arg);
    const Complex complex(double real, double imag);

    virtual void ParamInit() {    std::cout << "ParamInit" << std::endl;}
    virtual void ParticleInit() {    std::cout << "ParticleInit" << std::endl;}

  };

}

#endif
