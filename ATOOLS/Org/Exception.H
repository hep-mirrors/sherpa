#ifndef ATOOLS_Org_Exception_H
#define ATOOLS_Org_Exception_H

#include "ATOOLS/Org/Message.H"

namespace ATOOLS {

  class Exception {

  private:

    std::string m_info, m_class, m_method, m_type, m_file;
    int m_line;

  protected:

    const std::string& TypeName() const { return m_type; }

  public:

    Exception(const std::string& name="",
	      const std::string& info="",
	      const std::string& cmethod="",
        const std::string& file="",
        const int& line = -1);

    ~Exception() = default;

    friend std::ostream &operator<<(std::ostream &str,
				    const Exception &exception);

  };

  class normal_exit : public Exception {
  public:
    normal_exit(const std::string& info,
                const std::string& cmethod, const std::string& file,
                const int& line)
      : Exception("Normal exit",
		  info) {}
  };

  class inconsistent_option : public Exception {
  public:
    inconsistent_option(const std::string& info,
                        const std::string& cmethod, const std::string& file,
                        const int& line)
    : Exception("Inconsistent option",
		info, cmethod) {}
  };

  class not_implemented : public Exception {
  public:
    not_implemented(const std::string& info,
                    const std::string& cmethod, const std::string& file,
                    const int& line)
      : Exception("Not implemented",
		  info, cmethod) {}
  };

  class critical_error : public Exception {
  public:
    critical_error(const std::string& info,
                   const std::string& cmethod, const std::string& file,
                   const int& line)
      : Exception("Critical error",
		  info, cmethod, file, line) {}
  };

  class fatal_error : public Exception {
  public:
    fatal_error(const std::string& info,
                const std::string& cmethod, const std::string& file,
                const int& line)
      : Exception("Fatal error",
		  info, cmethod, file, line) {}
  };

  class missing_module : public Exception {
  public:
    missing_module(const std::string& info,
                   const std::string& cmethod, const std::string& file,
                   const int& line)
      : Exception("Missing module",
		  info, cmethod) {}
  };

  class missing_input : public Exception {
  public:
    missing_input(const std::string& info,
		              const std::string& cmethod, const std::string& file,
                  const int& line)
      : Exception("Missing input",
		  info, cmethod) {}
  };

  class invalid_input : public Exception {
  public:
    invalid_input(const std::string& info,
		              const std::string& cmethod, const std::string& file,
                  const int& line)
      : Exception("Invalid input",
		  info, cmethod) {}
  };

  class numerical_instability : public Exception {
  public:
    numerical_instability(const std::string& info,
			                    const std::string& cmethod, const std::string& file,
                          const int& line)
      : Exception("Numerical instability",
		  info, cmethod, file, line) {}
  };


  std::ostream &operator<<(std::ostream &str,const Exception &exception);

}

#ifdef __GNUC__
#define THROW(exception,message)			 \
  throw(ATOOLS::exception(message, __PRETTY_FUNCTION__, __FILE__ ,__LINE__));
#else
#define THROW(exception,message)					\
  throw(ATOOLS::exception(message,"<unknown class>::<unknown function>", __FILE__, __LINE__));
#endif

#include "ATOOLS/Org/Stacktrace.H"

#define DO_STACK_TRACE ATOOLS::GenerateStackTrace(std::cout)

#endif
