#ifndef SHERPA_PerturbativePhysics_Rescattering_Handler_H
#define SHERPA_PerturbativePhysics_Rescattering_Handler_H

#include "ALPACA/Main/Alpaca.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Org/Return_Value.H"

namespace SHERPA {

  class Rescattering_Handler {
  private:
    std::string m_name;

    ALPACA::Alpaca * p_alpaca;
  public:
    //constructor
    Rescattering_Handler();

    // destructor
    ~Rescattering_Handler();

    // member functions
    ATOOLS::Return_Value::code operator()(ATOOLS::Blob_List * blobs);    
    // inline functions
    inline std::string Name() { return m_name; }
  };
}

#endif
