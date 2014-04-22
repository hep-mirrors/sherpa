//%module Exception
%{
#include <ATOOLS/Org/Message.H>
#include <ATOOLS/Org/Exception.H>
#include <ATOOLS/Org/MyStrStream.H>
#include <signal.h>
#include <iostream>
  %}

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
    
  };// end of struct ex

  //std::ostream &operator<<(std::ostream &str,const ex::type &type);

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

    %extend {
      PyObject* __str__() {
	MyStrStream conv;
	conv<<*self;
	return PyString_FromString(conv.str().c_str());
      };
    };

  };// end of class Exception
  
}// end of namespace ATOOLS
