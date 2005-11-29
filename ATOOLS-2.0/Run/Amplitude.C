#include "Mom.H"
#include "Exception.H"
#include "MyStrStream.H"
#include "Color.H"
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
   
    int print=0;
    std::string param("");
    std::string file("Amp.dat");
    if (argc>1) {
      std::string option;
      for (int i=1;i<argc;i++) {
	option=argv[i];
	if (option=="-1") print=1;
	if (option=="-2") print=2;
	if (option=="-s") param="s";	
	//if (option=="-o") param="o";
      }
      option=argv[argc-1];
      if (option.find("-")!=0 && option.size()!=2) file=option;
    }
    MomentumList list(file.c_str());
    list.Print();
    if (param=="s"){
      Fullamplitude_MHV_test result(&list,list.GetHList(),list.GetPList());
      Complex amp(result.Calculate(print));
      msg_Info()<<endl<<endl<<amp<<endl; 
    }  
  
    //else if (param=="o"){
    //  Expression expression(list.Size(),&list,list.GetHList(),list.GetPList(),print);
    //  expression.Print();
    //  expression.Evaluate();
    //  std::cout<<"Amplitude: calculating -> "<<expression.Result()<<std::endl;
    //}
    else {
      Fullamplitude_MHV result(list.Size(),list.GetPList(),print);
      Complex amp(result.MSquare(&list,list.GetHList()));
      msg_Info()<<endl<<endl<<amp<<endl;
      //amp = result.MSquare(&list,list.GetHList());
      //msg_Info()<<endl<<endl<<amp<<endl;
    }
    return 0;
#ifdef MALLOC_TRACE
    muntrace();
#endif
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
