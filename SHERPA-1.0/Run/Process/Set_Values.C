#include "String_Handler.H"

#include "Message.H"
#include <dlfcn.h>

using namespace AMEGIC;
using namespace std;

typedef Values* (*Getter_Function)(Basic_Sfuncs*);

Values* String_Handler::Set_Values(std::string& pID,Basic_Sfuncs* BS)
{
  std::string libname=std::string("libProc_")+pID.substr(1)+std::string(".so");
  std::string gettername=std::string("Getter_")+pID;

  char * error;
  void * module;
  Getter_Function GetterFunction;

  // try loading library 
  module = dlopen(libname.c_str(),RTLD_LAZY);
  error  = dlerror();
  if (module==NULL) {
    ATOOLS::msg.Error()<<"String_Handler::Set_Values("<<pID<<","<<BS<<"): "
		       <<"Error while loading library "<<libname<<std::endl<<error<<std::endl;
    return 0;
  }

  GetterFunction = (Getter_Function)dlsym(module,gettername.c_str());
  error  = dlerror();
  if (error!=NULL) {
    ATOOLS::msg.Error()<<"String_Handler::Set_Values("<<pID<<","<<BS<<"): "
		       <<"Error while loading symbol from library "<<libname<<std::endl<<error<<std::endl;
    return 0;
  }

  ATOOLS::msg.Tracking()<<" calling Getter for library"<<endl;
  return GetterFunction(BS);
}
