#ifndef ATOOLS_Org_Exception_H
#define ATOOLS_Org_Exception_H

#include "ATOOLS/Org/Message.H"
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
      missing_input       = 7,
      missing_module      = 8,
      unknown             = 0 
    };
    
    friend std::ostream &operator<<(std::ostream &str,const ex::type &type);

  };// end of struct ex

  std::ostream &operator<<(std::ostream &str,const ex::type &type);

  class Exception_Handler;
  
  class Tester_Object {
  protected:

    virtual bool ApproveTerminate();
    friend class Exception_Handler;

  public:

    // destructor
    virtual ~Tester_Object();

  };// end of class Tester_Object

  class Terminator_Object {
  protected:

    virtual bool ReadInStatus(const std::string &path);
    virtual void PrepareTerminate();
    friend class Exception_Handler;
    
  public:

    // destructor
    virtual ~Terminator_Object();

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
	      std::string cmethod);
    Exception(const ex::type type,const std::string info,
	      const std::string cclass,const std::string cmethod);

    // destructor
    ~Exception();

    // inline functions
    inline void SetClass(const std::string cclass)     { m_class=cclass;    }
    inline void SetMethod(const std::string cmethod)   { m_method=cmethod;  }

    inline const std::string &Class() const   { return m_class;   }
    inline const std::string &Method() const  { return m_method;  }
    
    inline const std::string &Info() const { return m_info; }
    inline ex::type           Type() const { return m_type; }

    friend std::ostream &operator<<(std::ostream &str,
				    const Exception &exception);

  };// end of class Exception
  
  std::ostream &operator<<(std::ostream &str,const Exception &exception);

  /*!
    \class Exception
    \brief Provides exception handling

    This class provides exception handling. Each exception thrown 
    within Sherpa is caught by a catch (Exception exception) in Main.C.
    The appropriate way to include exception throwing into the code 
    is as follows
    \code
       using namespace ATOOLS;
       ...
       void Example_Class::Example_Function()
       {
         ...
	 // explicit call of constructor
	 throw(ex::normal_exit,"normal abort",
	       "Example_Class","Example_Function");
	 // call of constructor through #define directive 
	 THROW(normal_exit,"normal abort");
         ...
       }
    \endcode
  */

}// end of namespace ATOOLS

#if defined(__sgi) || defined(__GNUC__)
#define THROW(exception,message)			\
  throw(ATOOLS::Exception(ATOOLS::ex::exception,	\
			  message,__PRETTY_FUNCTION__));
#else
#define THROW(exception,message)					\
  throw(ATOOLS::Exception(ATOOLS::ex::exception,message,		\
			  "<unknown class>::<unknown function>"));
#endif

#include "ATOOLS/Org/Exception_Handler.H"

#endif
