#include "ATOOLS/Org/Return_Value.H"

#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"

using namespace ATOOLS;
using namespace std;


Counter_Map ATOOLS::Return_Value::s_warning_counter;
Counter_Map ATOOLS::Return_Value::s_error_counter;
Counter_Map ATOOLS::Return_Value::s_retry_method_counter;
Counter_Map ATOOLS::Return_Value::s_retry_phase_counter;
Counter_Map ATOOLS::Return_Value::s_new_event_counter;
Counter_Map ATOOLS::Return_Value::s_retry_event_counter;
Counter_Map ATOOLS::Return_Value::s_call_counter;

std::ostream &ATOOLS::operator<<(std::ostream &str,const Return_Value::code &rvc)
{
  switch (rvc) {
  case Return_Value::Error: return str<<"Error";
  case Return_Value::Failure: return str<<"Failure";
  case Return_Value::Undefined: return str<<"Undefined";
  case Return_Value::Success: return str<<"Success";
  case Return_Value::Nothing: return str<<"Nothing";
  case Return_Value::Warning: return str<<"Warning";
  case Return_Value::Retry_Method: return str<<"Retry_Method";
  case Return_Value::Retry_Phase: return str<<"Retry_Phase";
  case Return_Value::Retry_Event: return str<<"Retry_Event";
  case Return_Value::New_Event: return str<<"New_Event";
  }
  return str;
}

void Return_Value::PrintSingleStatistics(std::ostream &str,
					 const std::string &type,
					 const Counter_Map &map)
{
  if (!map.empty()){
    str<<"  "<<type<<":\n";
    for (Counter_Map::const_iterator it=map.begin();it!=map.end();it++) {
      unsigned long int calls(s_call_counter[it->first]);
      str<<"    From \""<<it->first<<"\": "<<it->second<<" ("
	 <<calls<<")";
      if (calls>0) str<<" -> "<<((it->second*1000)/calls)/10.0<<" %";
      str<<'\n';
    }
  }
}

void Return_Value::PrintStatistics(std::ostream &str)
{
  str<<"Event generation statistics:\n";
  str<<"  Generated events: "<<rpa->gen.NumberOfGeneratedEvents()<<endl;
  PrintSingleStatistics(str,"Errors",s_error_counter);
  PrintSingleStatistics(str,"Warnings",s_warning_counter);
  PrintSingleStatistics(str,"New events",s_new_event_counter);
  PrintSingleStatistics(str,"Retried events",s_retry_event_counter);
  PrintSingleStatistics(str,"Retried phases",s_retry_phase_counter);
  PrintSingleStatistics(str,"Retried methods",s_retry_method_counter);
}

void Return_Value::IncWarning(const std::string & name) {
  Counter_Map::iterator cit(s_warning_counter.find(name));
  if (cit!=s_warning_counter.end()) cit->second++;
  else s_warning_counter[name] = 1;
}

void Return_Value::IncError(const std::string & name) {
  Counter_Map::iterator cit(s_error_counter.find(name));
  if (cit!=s_error_counter.end()) cit->second++;
  else s_error_counter[name] = 1;
}

void Return_Value::IncRetryMethod(const std::string & name){
  Counter_Map::iterator cit(s_retry_method_counter.find(name));
  if (cit!=s_retry_method_counter.end()) cit->second++;
  else s_retry_method_counter[name] = 1;
}

void Return_Value::IncRetryPhase(const std::string & name) {
  Counter_Map::iterator cit(s_retry_phase_counter.find(name));
  if (cit!=s_retry_phase_counter.end()) cit->second++;
  else s_retry_phase_counter[name] = 1;
}

void Return_Value::IncNewEvent(const std::string & name) {
  Counter_Map::iterator cit(s_new_event_counter.find(name));
  if (cit!=s_new_event_counter.end()) cit->second++;
  else s_new_event_counter[name] = 1;
}

void Return_Value::IncRetryEvent(const std::string & name) {
  Counter_Map::iterator cit(s_retry_event_counter.find(name));
  if (cit!=s_retry_event_counter.end()) cit->second++;
  else s_retry_event_counter[name] = 1;
}

void Return_Value::IncCall(const std::string & name) {
  if (s_call_counter.find(name)!=s_call_counter.end())
    s_call_counter.find(name)->second++;
  else s_call_counter[name] = 1;
}

