#include "ATOOLS/Org/Library_Loader.H"

#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/My_MPI.H"

#include <dlfcn.h>
#include <fstream>
#include <unistd.h>
#include <iomanip>
#include <assert.h>

using namespace ATOOLS;

Library_Loader *ATOOLS::s_loader(NULL);

Library_Loader::Library_Loader()
{
  m_paths.push_back(rpa->gen.Variable("SHERPA_RUN_PATH"));
  m_paths.push_back(rpa->gen.Variable("SHERPA_LIBRARY_PATH"));
  const std::vector<std::string> &paths(EnvironmentVariable(LD_PATH_NAME));
  m_paths.insert(m_paths.end(),paths.begin(),paths.end());
}

void Library_Loader::UnloadLibrary(const std::string &name,void *module)
{
  std::map<std::string,void*>::iterator lit(m_libs.find(name));
  if (lit!=m_libs.end()) m_libs.erase(lit);
  dlclose(module);
}

bool Library_Loader::LibraryIsLoaded(const std::string &name)
{
  std::map<std::string,void*>::iterator lit(m_libs.find(name));
  if (lit!=m_libs.end()) return true;
  return false;
}

void *Library_Loader::LoadLibrary(const std::string &name)
{
  std::map<std::string,void*>::iterator lit(m_libs.find(name));
  if (lit!=m_libs.end()) return lit->second;
  msg_Debugging()<<METHOD<<"(lib"<<name<<LIB_SUFFIX<<") {"<<std::endl;
  int i(0);
#ifdef USING__MPI
  if (mpi->Rank()==0)
#endif
  for (;i<m_paths.size();++i) {
    std::string libname(m_paths[i]+"/lib"+name+LIB_SUFFIX);
    msg_Debugging()<<"  checking '"<<libname<<"'\n";
    void *module(LoadLibrary(m_paths[i],name));
    if (module==NULL) continue;
    msg_Debugging()<<"} found in '"<<m_paths[i]<<"'"<<std::endl;
#ifdef USING__MPI
    mpi->Bcast(&i,1,MPI_INT,0);
#endif
    m_libs[name]=module;
    return module;
  }
#ifdef USING__MPI
  mpi->Bcast(&i,1,MPI_INT,0);
  if (i<m_paths.size())
    return m_libs[name]=LoadLibrary(m_paths[i],name);
#endif
  msg_Debugging()<<"} failed"<<std::endl;
  msg_Info()<<METHOD<<"(): Failed to load library 'lib"
	    <<name<<LIB_SUFFIX<<"'."<<std::endl;
  return NULL;
}

void *Library_Loader::LoadLibrary(const std::string &path,
				  const std::string &name)
{
  std::string fullpath(path+"/lib"+name+LIB_SUFFIX);
  void *module(dlopen(fullpath.c_str(),RTLD_LAZY|RTLD_GLOBAL));
  return module;
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
	       <<funcname<<"'."<<std::endl;
    return NULL;
  }
  msg_Debugging()<<"done"<<std::endl;
  return func;
}

void *Library_Loader::GetLibraryFunction(const std::string &libname,
					 const std::string &funcname,
					 void *&module)
{
  msg_Debugging()<<"executing library function '"<<funcname
		 <<"' from 'lib"<<libname<<LIB_SUFFIX<<"' ... "<<std::flush;
  if (module==NULL) module=LoadLibrary(libname);
  if (module==NULL) return NULL;
  return GetLibraryFunction(funcname, module);
}

void *Library_Loader::GetLibraryFunction(const std::string &funcname,
					 void * const & module) const
{
  void *func(dlsym(module,funcname.c_str()));
  char *error(dlerror());
  if (error!=NULL) {
    msg_Debugging()<<"failed"<<std::endl;
    msg_Error()<<error<<std::endl;
    msg_Error()<<METHOD<<"(): Failed to load function '"
	       <<funcname<<"'."<<std::endl;
    return NULL;
  }
  msg_Debugging()<<"done"<<std::endl;
  return func;
}

void Library_Loader::AddPath(const std::string &path,const int mode)
{ 
  for (size_t i(0);i<m_paths.size();++i)
    if (m_paths[i]==path) return;
  if (mode) m_paths.push_back(path); 
  else m_paths.insert(m_paths.begin(),path);
}
