#include "ATOOLS/Org/Library_Loader.H"

#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/CXXFLAGS.H"

#include <sys/stat.h>
#include <dlfcn.h>
#include <fstream>
#include <unistd.h>
#include <iomanip>

using namespace ATOOLS;

Library_Loader *ATOOLS::s_loader(NULL);

Library_Loader::Library_Loader(): m_wait(3600)
{
  m_paths=EnvironmentVariable(LD_PATH_NAME);
  m_paths.push_back(rpa.gen.Variable("SHERPA_LIBRARY_PATH"));
}

bool Library_Loader::CreateLockFile(const std::string &lockname)
{
  msg_Debugging()<<"checking lock file '"<<lockname<<"' ... "<<std::flush;
  struct stat buffer;
  if (!stat(lockname.c_str(),&buffer)) {
    msg_Debugging()<<" found"<<std::endl;
    msg_Info()<<METHOD<<"(): Another process created '"<<lockname
	      <<"'. Waiting for unlock ...          "<<std::flush;
    size_t check(1), i(0);
    for (i=0;i<m_wait;i+=check) {
      sleep(check);
      msg_Info()<<mm(9,mm::left)<<std::setw(6)<<i<<" s "<<std::flush;
      if (stat(lockname.c_str(),&buffer)) break;
    }
    msg_Info()<<std::endl;
    if (i==m_wait && !stat(lockname.c_str(),&buffer)) {
      msg->Error()<<METHOD<<"(): '"<<lockname<<"' remains for "
		  <<m_wait<<" s. Timeout."<<std::endl;
      THROW(critical_error,"Library locked");
    }
  }
  msg_Debugging()<<" not found"<<std::endl;
  msg_Debugging()<<"creating lock file '"<<lockname<<"' ... "<<std::flush;
  std::ofstream *lock(new std::ofstream(lockname.c_str()));
  delete lock;
  msg_Debugging()<<" done"<<std::endl;
  return true;
}

bool Library_Loader::RemoveLockFile(const std::string &lockname)
{
  msg_Debugging()<<"deleting lock file '"<<lockname<<"' ... "<<std::flush;
  remove(lockname.c_str());
  msg_Debugging()<<" done"<<std::endl;
  return true;
}

void *Library_Loader::LoadLibrary(const std::string &name)
{
  std::map<std::string,void*>::iterator lit(m_libs.find(name));
  if (lit!=m_libs.end()) return lit->second;
  msg_Debugging()<<"loading library 'lib"<<name<<LIB_SUFFIX<<"' {"<<std::endl;
  for (size_t i(0);i<m_paths.size();++i) {
    std::string libname(m_paths[i]+"/lib"+name+LIB_SUFFIX);
    struct stat buffer;
    msg_Debugging()<<"checking for '"<<libname<<"' ... "<<std::flush;
    if (stat(libname.c_str(),&buffer)) {
      msg_Debugging()<<" not found"<<std::endl;
      continue;
    }
    msg_Debugging()<<" found"<<std::endl;
    std::string lockname(m_paths[i]+"/lib"+name+
			 std::string(LIB_SUFFIX)+".lock");
    if (!CreateLockFile(lockname)) return NULL;
    if (!CreateLockFile(rpa.gen.Variable("HOME")+
 			"/.sherpa/.liblock")) return NULL;
    void *module(dlopen(libname.c_str(),RTLD_LAZY|RTLD_GLOBAL));
    if (!RemoveLockFile(rpa.gen.Variable("HOME")+
 			"/.sherpa/.liblock")) return NULL;
    if (!RemoveLockFile(lockname)) return NULL;
    if (module!=NULL) {
      msg_Debugging()<<"} found in '"<<m_paths[i]<<"'"<<std::endl;
      m_libs[name]=module;
      return module;
    }
    char *err(dlerror());
    if (err!=NULL) {
      msg_Error()<<METHOD<<"(): "<<err<<std::endl;
      break;
    }
  }
  msg_Debugging()<<"} failed"<<std::endl;
  msg_Error()<<METHOD<<"(): Failed to load library 'lib"
	     <<name<<LIB_SUFFIX<<"'."<<std::endl;
  return NULL;
}

void *Library_Loader::GetLibraryFunction(const std::string &libname,
					 const std::string &funcname)
{
  msg_Debugging()<<"executing library function '"<<funcname
		 <<"' from 'lib"<<libname<<LIB_SUFFIX<<"' ... "<<std::flush;
  void *module(LoadLibrary(libname));
  if (module==NULL) return NULL;
  void *func(dlsym(module,funcname.c_str()));
  char *error(dlerror());
  if (error!=NULL) {
    msg_Debugging()<<"failed"<<std::endl;
    msg_Error()<<error<<std::endl;
    msg_Error()<<METHOD<<"(): Failed to load function '"
	       <<funcname<<"'."<<LD_PATH_NAME<<std::endl;
    return NULL;
  }
  msg_Debugging()<<"done"<<std::endl;
  return func;
}

void Library_Loader::AddPath(const std::string &path)
{ 
  for (size_t i(0);i<m_paths.size();++i)
    if (m_paths[i]==path) return;
  m_paths.push_back(path); 
}
