#include "String_Handler.H"

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
    cout<<" Error in loading library "<<libname<<endl;
    cout<<error<<endl;
    return 0;
  }

  GetterFunction = (Getter_Function)dlsym(module,gettername.c_str());
  error  = dlerror();
  if (error!=NULL) {
    cout<<" Error in loading symbol from library "<<endl;
    cout<<error<<endl;
    return 0;
  }

  cout<<" calling Getter for library"<<endl;
  return GetterFunction(BS);
}
