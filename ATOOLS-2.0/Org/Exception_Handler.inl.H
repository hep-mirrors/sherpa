#ifndef Exception_Handler_inl_H
#define Exception_Handler_inl_H

#include "Exception_Handler.H"

namespace ATOOLS {
  
  inline void Exception_Handler::AddTesterFunction(bool (*testerfunction)(void)) 
  { s_testerfunctions.push_back(testerfunction); }
  
  inline void Exception_Handler::AddTerminatorFunction(void (*terminatorfunction)(void))      
  { s_terminatorfunctions.push_back(terminatorfunction); }
  
  inline void Exception_Handler::AddTesterObject(Tester_Object *const testerobject)     
  { s_testerobjects.push_back(testerobject); }
  
  inline void Exception_Handler::AddTerminatorObject(Terminator_Object *const terminatorobject)     
  { s_terminatorobjects.push_back(terminatorobject); }

}

#endif
