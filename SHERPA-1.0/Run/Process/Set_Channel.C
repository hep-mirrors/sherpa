#include "Phase_Space_Generator.H"

#include <dlfcn.h>

using namespace AMEGIC;
using namespace PHASIC;
using namespace std;

typedef Single_Channel * (*Getter_Function)(int nin,int nout,ATOOLS::Flavour* fl,
						   int chn);

Single_Channel * Phase_Space_Generator::SetChannel(int nin,int nout,ATOOLS::Flavour* fl,
						   int chn,string& pID)
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
  return GetterFunction(nin,nout,fl,chn);
}
















































