#include "Algebra_Interpreter.H"
#include "Exception.H"
#include <termios.h>
#include <unistd.h>

using namespace ATOOLS;

struct TDouble: public Term {
  double m_value;
};// end of struct Double

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
    Algebra_Interpreter interpreter;
    std::string expr;
    bool quiet(false);
    for (int i=1;i<argc;++i) {
      std::string argvs=argv[i];
      if (argvs=="-q") {
	quiet=true;
	msg.SetLevel(0);
	continue;
      }
      size_t pos=argvs.find("=");
      if (pos!=std::string::npos && 
	  expr.length()>pos && expr[pos+1]!='=')
	interpreter.AddTag(argvs.substr(0,pos),argvs.substr(pos+1));
      else expr+=argv[i];
    }
    std::string result(interpreter.Interprete(expr));
    if (!quiet) {
      std::cout<<"Calc: interpreting formula -> "<<expr<<" = "
	       <<result<<std::endl;
    }
    else {
      std::cout<<result<<std::endl;
      return 0;
    }
    if (msg.LevelIsTracking()) interpreter.PrintEquation();
    std::cout.precision(12);
    std::cout<<"Calc: recalculating tree   -> "<<expr<<" = "
	     <<((TDouble*)interpreter.Calculate())->m_value<<std::endl;
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
