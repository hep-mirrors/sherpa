#include "Phase_Space_Generator.H"

#include <dlfcn.h>

using namespace AMEGIC;
using namespace PHASIC;
using namespace std;

typedef Single_Channel * (*Getter_Function)(int nin,int nout,ATOOLS::Flavour* fl);

Single_Channel * Phase_Space_Generator::SetChannel(int nin,int nout,ATOOLS::Flavour* fl,
						   string& pID)
{
  int pos=pID.find(string("/"));
  string libname=string("libProc_")+pID.substr(0,pos)+string(".so");
  string gettername=string("Getter_")+pID.substr(pos+1); 

  char * error;
  void * module;
  Getter_Function GetterFunction;

  // try loading library 
  module = dlopen(libname.c_str(),RTLD_LAZY);
  error  = dlerror();
  if (module==NULL) {
    ATOOLS::msg.Error()<<"Phase_Space_Generator::SetChannel("
		       <<nin<<","<<nout<<","<<fl<<","<<pID<<"): "
		       <<"Error in loading library "<<libname<<std::endl<<error<<std::endl;
    return 0;
  }

  GetterFunction = (Getter_Function)dlsym(module,gettername.c_str());
  error  = dlerror();
  if (error!=NULL) {
    ATOOLS::msg.Error()<<"Phase_Space_Generator::SetChannel("
		       <<nin<<","<<nout<<","<<fl<<","<<pID<<"): "
		       <<"Error while loading symbol from library "<<libname<<std::endl<<error<<std::endl;
    return 0;
  }
  return GetterFunction(nin,nout,fl);
}
















































