#define COMPILE__Debugger
#include "Debugger.H"
#undef COMPILE__Debugger

#include "Data_Reader.H"
#include "Message.H"

using namespace ATOOLS;

Value_Carrier_Base::Value_Carrier_Base():
  m_type("") {}

Value_Carrier_Base::~Value_Carrier_Base() {}

std::vector<Debugger*> Debugger::s_objects=std::vector<Debugger*>();
std::ostream *Debugger::s_output=&std::cout;
unsigned int Debugger::s_print=1;

Debugger ATOOLS::dbg=Debugger(0,"global");

Debugger::Debugger(long unsigned int address,const std::string &name): 
  m_address(address),
  m_name(name) 
{
  s_objects.push_back(this);
}

Debugger::~Debugger()
{
  for (std::vector<Debugger*>::iterator oit=s_objects.begin();
       oit!=s_objects.end();) {
    if (*oit==this) oit=s_objects.erase(oit);
    else ++oit;
  }
  for (std::map<std::string,Value_Carrier_Base*>::const_iterator 
	 vcit=m_values.begin();vcit!=m_values.end();++vcit) {
    delete vcit->second;
  }
}

void Debugger::SetOutputMode(const unsigned int print)
{
  if (print!=std::numeric_limits<unsigned int>::max()) {
    s_print=print;
    return;
  }
  Data_Reader *reader = new Data_Reader();
  if (!reader->ReadFromFile(s_print,"DEBUGGING_OUTPUT")) s_print=1;
  delete reader;
}

void Debugger::PrintStatus(const std::string &type,const std::string &function)
{
  if (s_print==0) return;
  *s_output<<om::bold<<"Debugging output for "<<om::reset<<"\""
	   <<om::green<<type<<om::reset<<"\""<<om::bold<<" at '"
	   <<om::blue<<function<<om::reset<<om::bold<<"': {"
	   <<om::reset<<std::endl;
  for (std::vector<Debugger*>::iterator oit=s_objects.begin();
       oit!=s_objects.end();++oit) {
    if ((*oit)->Type()==type || type.length()==0) *s_output<<*(*oit)<<std::endl;
  }
  *s_output<<om::bold<<"}"<<om::reset<<std::endl;
}

void Debugger::PrintMethodInfo(const std::string &method)
{
  if (s_print==0) return;
  *s_output<<om::blue<<method<<om::reset<<std::endl;
}

void Debugger::PrintLocalInfo(const std::string &method,const std::string &info)
{
  if (s_print==0) return;
  *s_output<<om::blue<<method<<om::reset<<":("	
	   <<om::green<<"\""<<info<<"\""<<om::reset<<")"<<std::endl;
}

std::ostream &ATOOLS::operator<<(std::ostream &ostr,const Debugger &debugger)
{
  ostr<<om::bold<<"   Debugger ["<<om::reset<<"("<<om::red<<debugger.m_address
      <<om::reset<<")"<<om::bold<<","<<om::reset<<"\""<<om::green<<debugger.m_name
      <<om::reset<<"\""<<om::bold<<"]: {"<<om::reset<<std::endl;
  for (std::map<std::string,Value_Carrier_Base*>::const_iterator 
	 vcit=debugger.m_values.begin();vcit!=debugger.m_values.end();++vcit) {
    ostr<<"      [\""<<om::green<<vcit->first<<om::reset<<"\",<"
	<<om::blue<<vcit->second->Type()<<om::reset<<">] -> [("<<om::red;
    vcit->second->PrintAddress(ostr);
    ostr<<om::reset<<"),'"<<om::green;
    vcit->second->PrintValue(ostr);
    ostr<<om::reset<<"']"<<std::endl; 
  }
  return ostr<<om::bold<<"   }"<<om::reset;
}


