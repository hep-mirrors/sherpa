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

    static bool          s_active, s_prepared, s_stacktrace, s_print;
    static unsigned int  s_signal, s_exitcode;
    static Exception    *s_exception;
    
    static unsigned int s_nbus, s_nsegv;

    static std::string s_progname;

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

    static void RemoveTesterFunction(bool (*function)(void));
    static void RemoveTerminatorFunction(void (*function)(void));
    static void RemoveTesterObject(Tester_Object *const object);
    static void RemoveTerminatorObject(Terminator_Object *const object);
    
    static void GenerateStackTrace(std::ostream &ostr,
				   const bool endline=true,
				   const std::string &comment="");

    // inline functions
    inline static void AddTesterFunction(bool (*function)(void))
    { s_testerfunctions.push_back(function); }
    inline static void AddTerminatorFunction(void (*function)(void))
    { s_terminatorfunctions.push_back(function); }
    inline static void AddTesterObject(Tester_Object *const object)
    { s_testerobjects.push_back(object); }
    inline static void AddTerminatorObject(Terminator_Object *const object)
    { s_terminatorobjects.push_back(object); }

    inline static void SetActive(const bool active)    
    { s_active=active; }
    inline static void SetStackTrace(const bool trace) 
    { s_stacktrace=trace;  }
    inline static void SetProgramName(const std::string &name)
    { s_progname=name; }

    inline static bool Active()     { return s_active;      }
    inline static bool StackTrace() { return s_stacktrace;  }

    inline static Exception    *LastException() { return s_exception; }
    inline static unsigned int  LastSignal()    { return s_signal;    }

  };// end of class Exception_Handler
  
}// end of namespace ATOOLS

#endif
#endif
