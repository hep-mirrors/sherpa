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

    size_t m_wait;
    bool   m_check;

    bool CreateLockFile(const std::string &lockname);
    bool RemoveLockFile(const std::string &lockname);

  public:

    // constructor
    Library_Loader();

    // member functions
    void *LoadLibrary(const std::string &name);

    void *GetLibraryFunction(const std::string &libname,
			     const std::string &funcname);

    // inline functions
    void AddPath(const std::string &path);

    inline void SetWaitTime(const size_t &wait) { m_wait=wait; }

    inline void SetCheck(const bool check) { m_check=check; }

  };//end of class Library_Loader

  extern Library_Loader *s_loader;

}// end of namsepace COMPARE

#endif
