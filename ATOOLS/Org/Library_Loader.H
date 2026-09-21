#ifndef ATOOLS_Org_Library_Loader_H
#define ATOOLS_Org_Library_Loader_H

#include <vector>
#include <map>
#include <string>

namespace ATOOLS {

  class Library_Loader {
  private:

    std::vector<std::string>    m_paths;
    std::map<std::string,void*> m_libs;

  public:

    // constructor
    Library_Loader();

    // member functions
    void *LoadLibrary(const std::string &name);
    void *LoadLibrary(const std::string &path,
		      const std::string &name);
    void  UnloadLibrary(const std::string &name,void *module);
    bool  LibraryIsLoaded(const std::string &name);

    void *GetLibraryFunction(const std::string &libname,
			     const std::string &funcname);
    void *GetLibraryFunction(const std::string &libname,
			     const std::string &funcname,
			     void *&module);
    void *GetLibraryFunction(const std::string &funcname,
			     void* const & module) const;

    // inline functions
    void AddPath(const std::string &path,const int mode=0);

  };//end of class Library_Loader

  extern Library_Loader *s_loader;

}// end of namsepace COMPARE

#endif
