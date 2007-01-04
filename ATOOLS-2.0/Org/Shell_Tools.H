#ifndef Shell_Tools_H
#define Shell_Tools_H

#include <string>
#include <vector>

namespace ATOOLS {

  bool MakeDir(std::string path,const mode_t mode,const bool create_tree=1);
  bool CopyFile(const std::string &oldname,const std::string &newname);

  std::vector<std::string> EnvironmentVariable(const std::string &name,
					       std::string entry="");

}// end of namespace ATOOLS

#endif
