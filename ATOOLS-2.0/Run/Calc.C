#include "Algebra_Interpreter.H"
#include "Exception.H"

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
    msg.SetModifiable(true);
    Algebra_Interpreter interpreter;
    std::string expr;
    for (int i=1;i<argc;++i) {
      std::string argvs=argv[i];
      size_t pos=argvs.find("=");
      if (pos==std::string::npos) expr+=argv[i];
      else interpreter.AddTag(argvs.substr(0,pos),
			      argvs.substr(pos+1));
    }
    std::cout<<"Calc: interpreting formula -> "<<expr<<" = "
	     <<interpreter.Interprete(expr)<<std::endl;
    if (msg.LevelIsDebugging()) interpreter.PrintEquation();
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
