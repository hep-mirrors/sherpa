#ifndef Model_Handler_2_h
#define Model_Handler_2_h

#include "Model_Base.H"
#include "Flavour.H"
#include "Data_Read.H"

namespace MODEL {
  class Model_Handler {
  public:
    Model_Handler() {};
    Model_Base * GetModel(AORGTOOLS::Data_Read *,std::string);
  };
}

#endif


