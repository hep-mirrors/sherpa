#include "String_Handler.H"

#include "Message.H"
#include "Run_Parameter.H"
#include <dlfcn.h>

using namespace AMEGIC;
using namespace std;

typedef Values* (*Getter_Function)(Basic_Sfuncs*);

Values* String_Handler::Set_Values(std::string& pID,Basic_Sfuncs* BS)
{
  std::string libname=ATOOLS::rpa.gen.Variable("SHERPA_LIB_PATH")+
    std::string("/libProc_")+pID.substr(1)+std::string(LIB_SUFFIX);
  std::string gettername=std::string("Getter_")+pID;

  char * error;
  void * module;
  Getter_Function GetterFunction;

  // try loading library 
  module = dlopen(libname.c_str(),RTLD_LAZY);
  error  = dlerror();
  if (module==NULL) return 0;

  GetterFunction = (Getter_Function)dlsym(module,gettername.c_str());
  error  = dlerror();
  if (error!=NULL) return 0;

  return GetterFunction(BS);
}
