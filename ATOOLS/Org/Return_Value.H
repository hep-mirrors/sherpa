#ifndef ATOOLS_Org_Return_Value_H
#define ATOOLS_Org_Return_Value_H

#include <string>
#include <map>

namespace ATOOLS {  

  typedef std::map<std::string,unsigned long int> Counter_Map;

  class Return_Value {
  private:
    static Counter_Map s_warning_counter;
    static Counter_Map s_error_counter;
    static Counter_Map s_retry_method_counter;
    static Counter_Map s_retry_phase_counter;
    static Counter_Map s_retry_event_counter;
    static Counter_Map s_retry_showerevent_counter;
    static Counter_Map s_new_event_counter;
    static Counter_Map s_call_counter;
    void PrintSingleStatistics(std::ostream &str,const std::string &type,
			       const Counter_Map &map) const;
  public:
    ~Return_Value();
    void Statistics() const;
    void PrintStatistics(std::ostream &str) const;
    void IncWarning(std::string);
    void IncError(std::string);
    void IncRetryMethod(std::string);
    void IncRetryPhase(std::string);
    void IncRetryShowerEvent(std::string);
    void IncRetryEvent(std::string);
    void IncNewEvent(std::string);
    void IncCall(std::string);
    enum code {
      Error        = -10,
      Failure      = -1,
      Undefined    = 0,
      Success      = 1,
      Nothing      = 2,
      Warning      = 10,
      Retry_Method = 20,
      Retry_Phase  = 21,
      Retry_Event  = 22,
      Retry_ShowerEvent  = 23,
      New_Event    = 100
    };
  };

  std::ostream &operator<<(std::ostream &str,const Return_Value::code &rvc);

  extern Return_Value rvalue;
}
#endif
