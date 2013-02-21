//%module Initialization_Handler
%{
#include "SHERPA/Initialization/Initialization_Handler.H"
%}

%catches (ATOOLS::Exception) SHERPA::Initialization_Handler::Initialization_Handler(int, char**);
%catches (ATOOLS::Exception) SHERPA::Initialization_Handler::~Initialization_Handler();
%catches (ATOOLS::Exception) SHERPA::Initialization_Handler::GetMatrixElementHandler();


namespace ATOOLS {
  class Data_Reader;
}

namespace SHERPA {

  class Input_Output_Handler;
  class Matrix_Element_Handler;
  class Hard_Decay_Handler;
  class Beam_Remnant_Handler;
  class Fragmentation_Handler;
  class Decay_Handler_Base;
  class MI_Handler;
  class Soft_Photon_Handler;
  class Lund_Interface;
  class Event_Reader_Base;
  class Analysis_Interface;
  class Soft_Collision_Handler;

  class Initialization_Handler: public ATOOLS::Terminator_Object {

  public :

    Initialization_Handler(int argc,char * argv[]);
    ~Initialization_Handler();
    Matrix_Element_Handler       * GetMatrixElementHandler()  const { return p_mehandler; }
  };
}
