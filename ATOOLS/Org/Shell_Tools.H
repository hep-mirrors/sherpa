#ifndef ATOOLS_Org_Shell_Tools_H
#define ATOOLS_Org_Shell_Tools_H

#include <string>
#include <vector>
#include <sys/types.h>
#include <unistd.h> 

namespace ATOOLS {

  bool MakeDir(std::string path,const bool create_tree=1,
	       const mode_t mode=0777);
  bool ChMod(const std::string &file,const mode_t mode=0777);
  bool Copy(const std::string &oldname,const std::string &newname,
	    const bool recursive=false);
  bool Move(const std::string &oldname,const std::string &newname);
  bool Remove(const std::string &file,const bool recursive=false);

  std::vector<std::string> EnvironmentVariable(const std::string &name,
					       std::string entry="");

  bool FileExists(const std::string &file);
  bool DirectoryExists(const std::string &dir);

  std::string Demangle(const std::string &name);

  std::string GetCWD();

  std::vector<std::string> RegExMatch
  (const std::string &str,const std::string &pat,const size_t nm=1);

}// end of namespace ATOOLS

#endif
