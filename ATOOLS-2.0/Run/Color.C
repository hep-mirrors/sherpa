#include "Color.H"
#include "Exception.H"
#include "Data_Reader.H"
#include <termios.h>
#include <unistd.h>

using namespace COLOR;
using namespace ATOOLS;

int main(int argc,char **argv)
{
  std::set_terminate(Exception_Handler::Terminate);
  std::set_unexpected(Exception_Handler::Terminate);
  signal(SIGSEGV,Exception_Handler::SignalHandler);
  signal(SIGINT,Exception_Handler::SignalHandler);
  signal(SIGBUS,Exception_Handler::SignalHandler);
  signal(SIGFPE,Exception_Handler::SignalHandler);
  signal(SIGABRT,Exception_Handler::SignalHandler);
  signal(SIGTERM,Exception_Handler::SignalHandler);
  signal(SIGXCPU,Exception_Handler::SignalHandler);
  try {
    msg.Init(2,"");
    termios testos;
    if (tcgetattr(STDOUT_FILENO,&testos)==0) msg.SetModifiable(true);
    else msg.SetModifiable(false);
    Expression expression("");
    std::string argvs;
    for (int i=1;i<argc;++i) argvs+=std::string(" ")+argv[i];
    Data_Reader reader;
    reader.SetString(argvs);
    int level(2);
    if (reader.ReadFromString(level,"-O")) msg.SetLevel(level);
//     std::string expr;
//     for (int i=1;i<argc;++i) {
//       std::string argvs=argv[i];
//       size_t pos=argvs.find("=");
//       if (pos!=std::string::npos && 
// 	  expr.length()>pos && expr[pos+1]!='=')
// 	interpreter.AddTag(argvs.substr(0,pos),argvs.substr(pos+1));
//       else expr+=argv[i];
//     }
    expression.Print();
    expression.Evaluate();
    std::cout<<"Color: calculating -> "<<" = "
	     <<expression.Result()<<std::endl;
    return 0;
  }
  catch (Exception exception) {
    exception.UpdateLogFile();
    msg.Error()<<exception<<std::endl;
    std::terminate();
  }
  catch (std::exception exception) {
    std::cout<<"Sherpa: throws std::exception "
	     <<exception.what()<<" ..."<<std::endl;
    std::terminate();
  }
}
