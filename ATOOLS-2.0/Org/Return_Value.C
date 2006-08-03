#include "Return_Value.H"
#include "Message.H"

using namespace ATOOLS;
using namespace std;


std::map<std::string,unsigned long int> ATOOLS::Return_Value::s_warning_counter;
std::map<std::string,unsigned long int> ATOOLS::Return_Value::s_error_counter;
std::map<std::string,unsigned long int> ATOOLS::Return_Value::s_retry_method_counter;
std::map<std::string,unsigned long int> ATOOLS::Return_Value::s_retry_phase_counter;
std::map<std::string,unsigned long int> ATOOLS::Return_Value::s_retry_event_counter;
ATOOLS::Return_Value ATOOLS::rvalue;

void Return_Value::Statistics()
{
  msg.Out()<<METHOD<<": Return value statistics: "<<endl<<endl
	   <<"   Retry events: "<<endl;
  for (map<string,unsigned long int>::iterator it=s_retry_event_counter.begin();
       it!=s_retry_event_counter.end();it++)
    msg.Out()<<"   "<<it->first<<" : "<<it->second<<endl;
  msg.Out()<<endl<<"   Retry phases: "<<endl;
  for (map<string,unsigned long int>::iterator it=s_retry_phase_counter.begin();
       it!=s_retry_phase_counter.end();it++)
    msg.Out()<<"   "<<it->first<<" : "<<it->second<<endl;
  msg.Out()<<endl<<"   Retry methods: "<<endl;
  for (map<string,unsigned long int>::iterator it=s_retry_method_counter.begin();
       it!=s_retry_method_counter.end();it++)
    msg.Out()<<"   "<<it->first<<" : "<<it->second<<endl;
}

void Return_Value::IncWarning(std::string name) {
  if (s_warning_counter.find(name)!=s_warning_counter.end())
    s_warning_counter.find(name)->second++;
  else s_warning_counter[name] = 1;
}

void Return_Value::IncError(std::string name) {
  if (s_error_counter.find(name)!=s_error_counter.end())
    s_error_counter.find(name)->second++;
  else s_error_counter[name] = 1;
}

void Return_Value::IncRetryMethod(std::string name){
  if (s_retry_method_counter.find(name)!=s_retry_method_counter.end())
    s_retry_method_counter.find(name)->second++;
  else s_retry_method_counter[name] = 1;
}

void Return_Value::IncRetryPhase(std::string name) {
  if (s_retry_phase_counter.find(name)!=s_retry_phase_counter.end())
    s_retry_phase_counter.find(name)->second++;
  else s_retry_phase_counter[name] = 1;
}

void Return_Value::IncRetryEvent(std::string name) {
  if (s_retry_event_counter.find(name)!=s_retry_event_counter.end())
    s_retry_event_counter.find(name)->second++;
  else s_retry_event_counter[name] = 1;
}

