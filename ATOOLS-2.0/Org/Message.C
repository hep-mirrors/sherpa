#include "Message.H"
#include "Run_Parameter.H"

namespace ATOOLS {
  Message       msg;
}

using namespace ATOOLS;


Message::Message() {      
  output = &std::cout;
  error  = &std::cerr;
  no     = new std::ofstream("/dev/null",std::ios::app);
  file   = 0;
  level  = 0;
}

void Message::SetNoStream(std::ostream * _no) 
{
  no     = (std::ofstream*)_no;
};

void Message::SetOutStream(std::ostream * _output) {
  output = _output;
};

void Message::SetErrStream(std::ostream * _error) {
  error  = _error;
};



void Message::Init(int _level) { 
  level = _level; 
  if (level>2) (*output)<<"Initialize output module Message. Level "<<_level<<std::endl;
}

void Message::SetFile(char* name="Amegic.log") {
  if (file==0) {
    output = new std::ofstream(name);
    error  = output;
  }
  file = 1;
}

void Message::SetStandard() {
  if (file==1) delete output;
  
  output = &std::cout;
  error  = &std::cerr;
  file   = 0;
}

void Message::SetPrecision(int prec) {
  if (output) output->precision(prec);
}
