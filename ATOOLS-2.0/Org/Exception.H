#ifndef Exception_H
#define Exception_H

#include "Message.H"
#include <exception>
#include <signal.h>

namespace ATOOLS {

  struct ex {
    
    enum type {
      normal_exit         = 1,
      unknown_option      = 2,
      inconsistent_option = 3,
      not_implemented     = 4,
      critical_error      = 5,
      fatal_error         = 6,
      unknown             = 0 
    };
    
    friend std::ostream &operator<<(std::ostream &str,const ex::type &type);

  };// end of struct ex

  std::ostream &operator<<(std::ostream &str,const ex::type &type);

  class Exception_Handler;
  
  class Tester_Object {
  private:

    virtual bool ApproveTerminate();
    friend class Exception_Handler;

  };// end of class Tester_Object

  class Terminator_Object {
  private:

    virtual void PrepareTerminate();
    friend class Exception_Handler;

  };// end of class Terminator_Object

  class Exception {
  private:

    ex::type    m_type;
    std::string m_info, m_class, m_method;

    friend class Exception_Handler;

  public:

    // constructors
    Exception(const ex::type type,const std::string info);
    Exception(const ex::type type,const std::string info,
	      const std::string cclass,const std::string cmethod);

    // destructor
    ~Exception();

    // member functions
    void UpdateLogFile() const;
    
    // inline functions
    inline void SetClass(const std::string cclass)     { m_class=cclass;    }
    inline void SetMethod(const std::string cmethod)   { m_method=cmethod;  }

    inline const std::string &Class() const   { return m_class;   }
    inline const std::string &Method() const  { return m_method;  }
    
    inline const std::string &Info() const { return m_info; }
    inline ex::type           Type() const { return m_type; }

    friend std::ostream &operator<<(std::ostream &str,const Exception &exception);

  };// end of class Exception
  
  std::ostream &operator<<(std::ostream &str,const Exception &exception);

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

}// end of namespace ATOOLS

#include "Exception_Handler.H"

#endif
