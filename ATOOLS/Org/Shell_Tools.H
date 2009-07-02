#ifndef ATOOLS_Org_Shell_Tools_H
#define ATOOLS_Org_Shell_Tools_H

#include <string>
#include <vector>
#include <sys/types.h>

namespace ATOOLS {

  bool MakeDir(std::string path,const bool create_tree=1,
	       const mode_t mode=0777);
  bool CopyFile(const std::string &oldname,const std::string &newname);
  bool MoveFile(const std::string &oldname,const std::string &newname);

  std::vector<std::string> EnvironmentVariable(const std::string &name,
					       std::string entry="");

  bool FileExists(const std::string &file);

  std::string Demangle(const std::string &name);

}// end of namespace ATOOLS

#endif
