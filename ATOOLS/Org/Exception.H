#ifndef ATOOLS_Org_Exception_H
#define ATOOLS_Org_Exception_H

#include "ATOOLS/Org/Message.H"
#include <exception>

namespace ATOOLS {

  class Exception {

  private:

    std::string m_info, m_class, m_method, m_type;

  protected:

    const std::string& TypeName() const { return m_type; }

  public:

    Exception(const std::string& name="",
	      const std::string& info="",
	      const std::string& cmethod="");

    ~Exception() {}

    friend std::ostream &operator<<(std::ostream &str,
				    const Exception &exception);

  };

  class normal_exit : public Exception {
  public:
    normal_exit(const std::string info,
		std::string cmethod)
      : Exception("Normal exit",
		  info, cmethod) {}
  };
  
  class inconsistent_option : public Exception {
  public:
    inconsistent_option(const std::string info,
			std::string cmethod)
    : Exception("Inconsistent option",
		info, cmethod) {}
  };

  class not_implemented : public Exception {
  public:
    not_implemented(const std::string info,
		    std::string cmethod)
      : Exception("Not implemented",
		  info, cmethod) {}
  };

  class critical_error : public Exception {
  public:
    critical_error(const std::string info,
		   std::string cmethod)
      : Exception("Critical error",
		  info, cmethod) {}
  };

  class fatal_error : public Exception {
  public:
    fatal_error(const std::string info,
		   std::string cmethod)
      : Exception("Fatal error",
		  info, cmethod) {}
  };

  class missing_module : public Exception {
  public:
    missing_module(const std::string info,
		   std::string cmethod)
      : Exception("Missing module",
		  info, cmethod) {}
  };

  class missing_input : public Exception {
  public:
    missing_input(const std::string info,
		  std::string cmethod)
      : Exception("Missing input",
		  info, cmethod) {}
  };

  class invalid_input : public Exception {
  public:
    invalid_input(const std::string info,
		  std::string cmethod)
      : Exception("Invalid input",
		  info, cmethod) {}
  };

  class numerical_instability : public Exception {
  public:
    numerical_instability(const std::string info,
			  std::string cmethod)
      : Exception("Numerical instability",
		  info, cmethod) {}
  };

  
  std::ostream &operator<<(std::ostream &str,const Exception &exception);

}

#if defined(__GNUC__)
#define THROW(exception,message)			 \
  throw(ATOOLS::exception(message,__PRETTY_FUNCTION__)); 
#else
#define THROW(exception,message)					\
  throw(ATOOLS::exception(message,"<unknown class>::<unknown function>")); 
#endif

#include "ATOOLS/Org/Stacktrace.H"

#define DO_STACK_TRACE ATOOLS::GenerateStackTrace(std::cout)

#endif
