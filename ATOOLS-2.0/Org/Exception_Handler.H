#ifndef Exception_Handler_H
#define Exception_Handler_H

#ifndef Exception_H
#error The header 'Exception_Handler.H' must not be included directly. \
Please include the header 'Exception.H' instead.
#else

#include <vector>

namespace ATOOLS {

  class Exception_Handler {
  public:

    typedef bool (*Tester_Function)(void);
    typedef void (*Terminator_Function)(void);
    
  private:

    static bool          s_active, s_prepared;
    static unsigned int  s_exitcode;
    static Exception    *s_exception;
    
    static unsigned int s_nbus, s_nsegv;

    static std::vector<Tester_Function>     s_testerfunctions;
    static std::vector<Terminator_Function> s_terminatorfunctions;

    static std::vector<Tester_Object*>     s_testerobjects;
    static std::vector<Terminator_Object*> s_terminatorobjects;

    static void SetExitCode();

    friend class Exception;
    
  public:

    // member functions
    static void SignalHandler(int signal); 
    
    static bool ApproveTerminate();
    static void PrepareTerminate();
    static void Exit(int exitcode);

    static void Terminate();

    static void RemoveTesterObject(Tester_Object *const object);
    static void RemoveTerminatorObject(Terminator_Object *const object);
    
    // inline functions
    static void AddTesterFunction(bool (*testerfunction)(void)); 
    static void AddTerminatorFunction(void (*terminatorfunction)(void));      
    static void AddTesterObject(Tester_Object *const object);     
    static void AddTerminatorObject(Terminator_Object *const object);     

    inline static void SetActive(const bool active) { s_active=active;  }
    inline static bool Active()                     { return s_active;  }

  };// end of class Exception_Handler
  
}// end of namespace ATOOLS

#include "Exception_Handler.inl.H"

#endif
#endif
