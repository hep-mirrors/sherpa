#ifndef Model_Handler_2_h
#define Model_Handler_2_h

#include "Model_Base.H"
#include "Flavour.H"
#include "Data_Read.H"

namespace MODEL {
  class Model_Handler {
  public:
    Model_Handler() {};
    Model_Base * GetModel(AORGTOOLS::Data_Read *,std::string,std::string);
  };

  /*!
    \class Model_Handler
    \brief to select and return a specific model.
    
    This class consists of only one method which is used to select and return
    a model dependent on the data read in through the AORGTOOLS::Data_Read object
    initialized outside.
  */
  /*!
    \fn Model_Base * Model_Handler::GetModel(AORGTOOLS::Data_Read *,std::string)
    The only method of this class, to select a model and return it.
    Through the AORGTOOLS::Data_Read object two specifications are read in, namely 
    the name of the model (default is "SM") and hte data file associated with it
    (default is "Model.dat"). Then the accorsding model is initialized and the 
    name of the data file is passed into its constructor to read in further
    parameters etc..
  */
}

#endif


