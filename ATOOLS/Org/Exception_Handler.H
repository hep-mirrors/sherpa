#ifndef ATOOLS_Org_Exception_Handler_H
#define ATOOLS_Org_Exception_Handler_H

#ifndef ATOOLS_Org_Exception_H
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

    bool          m_active, m_prepared, m_stacktrace;
    bool          m_print, m_noremove;
    unsigned int  m_signal, m_exitcode;
    Exception    *m_exception;
    
    unsigned int m_nbus, m_nsegv;

    std::string m_progname;

    std::vector<Tester_Function>     m_testerfunctions;
    std::vector<Terminator_Function> m_terminatorfunctions;

    std::vector<Tester_Object*>     m_testerobjects;
    std::vector<Terminator_Object*> m_terminatorobjects;

    void SetExitCode();

    friend class Exception;
    
    // constructor
    Exception_Handler();

  public:

    ~Exception_Handler();

    // member functions
    void Init();
    void Reset();

    bool ReadInStatus(const std::string &path);
    bool ApproveTerminate();
    void PrepareTerminate();
    void Exit(int exitcode);
    
    void SignalHandler(int signal); 
    void Terminate();

    void GenerateStackTrace(std::ostream &ostr,
			    const bool endline=true,
			    const std::string &comment="");

    void AddTesterFunction(bool (*function)(void));
    void AddTerminatorFunction(void (*function)(void));
    void AddTesterObject(Tester_Object *const object);
    void AddTerminatorObject(Terminator_Object *const object);

    void RemoveTesterFunction(bool (*function)(void));
    void RemoveTerminatorFunction(void (*function)(void));
    void RemoveTesterObject(Tester_Object *const object);
    void RemoveTerminatorObject(Terminator_Object *const object);
    
    void SetActive(const bool active);    
    void SetStackTrace(const bool trace); 
    void SetProgramName(const std::string &name);

    bool Active();
    bool StackTrace();

    Exception    *LastException();
    unsigned int  LastSignal();
    std::string   ProgramName();

  };// end of class Exception_Handler

  extern Exception_Handler *exh;

  void SignalHandler(int signal); 
  void Terminate();
  
}// end of namespace ATOOLS

#endif
#endif
