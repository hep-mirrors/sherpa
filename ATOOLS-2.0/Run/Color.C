#include "Color.H"
#include "Exception.H"
#include "MyStrStream.H"
#include <termios.h>
#include <unistd.h>
#ifdef MALLOC_TRACE
#error MALLOC_TRACE should not be defined
#include <mcheck.h>
#endif

using namespace ATOOLS;

int main(int argc,char **argv)
{
#ifdef MALLOC_TRACE
  setenv("MALLOC_TRACE","malloc_trace_color.log",1);
  mtrace();
#endif
  std::set_terminate(exh->Terminate);
  std::set_unexpected(exh->Terminate);
  signal(SIGSEGV,exh->SignalHandler);
  signal(SIGINT,exh->SignalHandler);
  signal(SIGBUS,exh->SignalHandler);
  signal(SIGFPE,exh->SignalHandler);
  signal(SIGABRT,exh->SignalHandler);
  signal(SIGTERM,exh->SignalHandler);
  signal(SIGXCPU,exh->SignalHandler);
  try {
    msg->Init(2,"");
    termios testos;
    if (tcgetattr(STDOUT_FILENO,&testos)==0) msg->SetModifiable(true);
    else msg->SetModifiable(false);
    std::string expr;
    for (int i=1;i<argc;++i) {
      std::string argvs=argv[i];
      size_t pos=argvs.find("-O");
      if (pos!=std::string::npos && argvs.length()>pos+2)
	msg->SetLevel(ToType<int>(argvs.substr(pos+2)));
      else expr+=argv[i];
    }
    Expression expression(expr);
    expression.Print();
    expression.Evaluate();
    std::cout<<"Color: calculating -> "<<expr<<" = "
	     <<expression.Result()<<std::endl;
    return 0;
#ifdef MALLOC_TRACE
    muntrace();
#endif
  }
  catch (Exception exception) {
    exception.UpdateLogFile();
    msg_Error()<<exception<<std::endl;
    std::terminate();
  }
  catch (std::exception exception) {
    std::cout<<"Sherpa: throws std::exception "
	     <<exception.what()<<" ..."<<std::endl;
    std::terminate();
  }
}
