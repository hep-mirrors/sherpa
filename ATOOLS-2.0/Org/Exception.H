#ifndef Exception_H
#define Exception_H

#include "Message.H"
#include <vector>
#include <exception>
#include <signal.h>

namespace ATOOLS {

  struct ex {
    
    enum type {
      normal_exit    = 1,
      critical_error = 2,
      fatal_error    = 3,
      unknown        = 0 
    };
    
    friend std::ostream &operator<<(std::ostream &str,const ex::type &type);

  };// end of struct ex

  std::ostream &operator<<(std::ostream &str,const ex::type &type);

  class Exception;

  class Terminator_Base {
  private:

    virtual void Terminate();
    
    friend class Exception;

  };// end of class terminator

  class Exception {
  public:

    typedef bool (*Tester)(const Exception &);
    typedef void (*Terminator)(void);
    
  private:

    static unsigned int s_exitcode;

    static std::vector<Tester>     s_testers;
    static std::vector<Terminator> s_terminators;

    static std::vector<Terminator_Base*> s_objects;
    
    ex::type    m_type;
    std::string m_info, m_class, m_method;

    void SetExitCode();

    friend void Terminate();

  public:

    // constructors
    Exception(const ex::type type,const std::string info);
    Exception(const ex::type type,const std::string info,
	      const std::string cclass,const std::string cmethod);

    // member functions
    static void SignalHandler(int signal); 
    static void PrepareTerminate();
    static void Terminate();
    static void Exit(int exitcode);

    void UpdateLogFile() const;
    bool ApproveTerminate() const;

    static void RemoveObject(Terminator_Base *const object);
    
    // inline functions
    inline static void AddTester(bool (*tester)(const Exception &)) { s_testers.push_back(tester);         }
    inline static void AddTerminator(void (*terminator)(void))      { s_terminators.push_back(terminator); }
    inline static void AddObject(Terminator_Base *const object)     { s_objects.push_back(object);         }

    inline void SetClass(const std::string cclass)     { m_class=cclass;    }
    inline void SetMethod(const std::string cmethod)   { m_method=cmethod;  }

    inline const std::string &Class() const   { return m_class;   }
    inline const std::string &Method() const  { return m_method;  }
    
    inline const std::string &Info() const { return m_info; }

    inline ex::type Type() const { return m_type; }

    friend std::ostream &operator<<(std::ostream &str,const Exception &exception);

  };// end of class Exception
  
  std::ostream &operator<<(std::ostream &str,const Exception &exception);

}// end of namespace ATOOLS

  /*!
    \class Exception
    \brief Provides exception handling

    This class provides exception handling. Each exception thrown within Sherpa is caught
    by a catch (Exception exception) in Main.C.
    The appropriate way to include exception throwing into the code is as follows
    \code
       using namespace ATOOLS;
       ...
       void Example_Class::Example_Function
       {
         ...
	 throw(ex::normal_exit,"normal abort","Example_Class","Example_Function");
         ...
       }
    \endcode
  */

#endif
