// chooses the model

#ifndef Model_Handler_h
#define Model_Handler_h

#include "Interaction_Model_Base.H"

namespace MODEL {
  class Interaction_Model_Handler {
    MODEL::Model_Base * p_model;
  public:
    Interaction_Model_Handler(MODEL::Model_Base *);
    Interaction_Model_Base * GetModel(std::string,std::string,std::string);
  };
}

#endif


