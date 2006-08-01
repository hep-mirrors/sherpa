#ifndef Return_Value_H
#define Return_Value_H

#include<string>
#include<map>

namespace ATOOLS {  

  class Return_Value {
  private:
    static std::map<std::string,unsigned long int> s_warning_counter;
    static std::map<std::string,unsigned long int> s_error_counter;
    static std::map<std::string,unsigned long int> s_retry_method_counter;
    static std::map<std::string,unsigned long int> s_retry_phase_counter;
    static unsigned long int                       s_retry_event_counter;
  public:
    void Statistics();
    void IncWarnings(std::string);
    void IncError(std::string);
    void IncRetryMethod(std::string);
    void IncRetryPhase(std::string);
    void IncRetryEvent(std::string);
    enum code {
      Error        = -10,
      Failure      = -1,
      Undefined    = 0,
      Success      = 1,
      Nothing      = 2,
      Warning      = 10,
      Retry_Method = 20,
      Retry_Phase  = 21,
      Retry_Event  = 22
    };
  };

  extern Return_Value rvalue;
}
#endif
