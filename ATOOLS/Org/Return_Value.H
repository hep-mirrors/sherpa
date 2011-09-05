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
    static Counter_Map s_new_event_counter;
    static Counter_Map s_call_counter;
    static void PrintSingleStatistics
    (std::ostream &str,const std::string &type,
     const Counter_Map &map);
  public:
    static void PrintStatistics(std::ostream &str);
    static void IncWarning(const std::string &);
    static void IncError(const std::string &);
    static void IncRetryMethod(const std::string &);
    static void IncRetryPhase(const std::string &);
    static void IncRetryEvent(const std::string &);
    static void IncNewEvent(const std::string &);
    static void IncCall(const std::string &);
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
      New_Event    = 100
    };
  };

  std::ostream &operator<<(std::ostream &str,const Return_Value::code &rvc);

}
#endif
