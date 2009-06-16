#include "AMEGIC++/Amplitude/Zfunctions/Mom.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Phys/Color.H"
#include <termios.h>
#include <unistd.h>
#include <fstream.h>
#include <iostream.h>
#ifdef MALLOC_TRACE
#error MALLOC_TRACE should not be defined
#include <mcheck.h>
#endif

using namespace AMEGIC;
using namespace ATOOLS;

int main(int argc, char *argv[])
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
   
    int print=0;
    bool single=false;
    std::string file("Amp.dat");
    if (argc>1) {
      std::string option;
      for (int i=1;i<argc;i++) {
	option=argv[i];
	if (option=="-1") print=1;
	if (option=="-2") print=2;
	if (option=="-s") single=true;
      }
      option=argv[argc-1];
      if (option.find("-")!=0 && option.size()!=2) file=option;
    }
    MomentumList list(file.c_str());
    list.Print();
    if (single){
      Fullamplitude_MHV result(&list,list.GetHList(),list.GetPList());
      Complex amp(result.Calculate(print));
      msg_Info()<<endl<<endl<<amp<<endl;
    }
    else {
      Expression expression(list.Size(),&list,list.GetHList(),list.GetPList(),print);
      //expression.Print();
      expression.Evaluate();
      std::cout<<"Amplitude: calculating -> "<<expression.Result()<<std::endl;
    }
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
