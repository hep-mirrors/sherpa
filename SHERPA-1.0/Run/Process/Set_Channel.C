#include "Phase_Space_Generator.H"

#include "Run_Parameter.H"
#include <dlfcn.h>

using namespace AMEGIC;
using namespace PHASIC;
using namespace std;

typedef Single_Channel * (*Getter_Function)(int nin,int nout,ATOOLS::Flavour* fl);

Single_Channel * Phase_Space_Generator::SetChannel(int nin,int nout,ATOOLS::Flavour* fl,
						   string& pID)
{
  int pos=pID.find(string("/"));
  string libname=ATOOLS::rpa.gen.Variable("SHERPA_LIB_PATH")+
    string("/libProc_")+pID.substr(0,pos)+string(".so");
  string gettername=string("Getter_")+pID.substr(pos+1); 

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

  return GetterFunction(nin,nout,fl);
}
