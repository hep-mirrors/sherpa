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
    
    unsigned int m_nbus, m_nsegv, m_mpi_timeout;

    std::string m_progname;

    std::vector<Tester_Function>     m_testerfunctions;
    std::vector<Terminator_Function> m_terminatorfunctions;

    std::vector<Tester_Object*>     m_testerobjects;
    std::vector<Terminator_Object*> m_terminatorobjects;

    std::vector<int> m_mpi;

    void SetExitCode();

    friend class Exception;
    
  public:

    Exception_Handler();

    ~Exception_Handler();

    // member functions
    void Reset();
    void MPISync();

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
    
    inline void SetActive(const bool active)            { m_active=active;    }
    inline void SetStackTrace(const bool trace)         { m_stacktrace=trace; }
    inline void SetProgramName(const std::string &name) { m_progname=name;    }
    inline void SetMPITimeOut(const int t)              { m_mpi_timeout=t;    }

    inline bool Active()     { return m_active;     }
    inline bool StackTrace() { return m_stacktrace; }

    inline Exception    *LastException() { return m_exception; }
    inline unsigned int  LastSignal()    { return m_signal;    }
    inline std::string   ProgramName()   { return m_progname;  }

    int MPIStat(int rank);

  };// end of class Exception_Handler

  extern Exception_Handler *exh;

  void SignalHandler(int signal); 
  void Terminate();
  
}// end of namespace ATOOLS

#endif
#endif
