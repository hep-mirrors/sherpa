#include "Shell_Tools.H"

#include <sys/stat.h>
#include <errno.h>
#include <dirent.h>

#ifdef DEBUG__Shell_Tools
#include <iostream>
#endif

using namespace ATOOLS;

bool ATOOLS::MakeDir(std::string path,const mode_t mode)
{
  if (path=="") return false;
#ifdef DEBUG__Shell_Tools
  std::cout<<"ATOOLS::MakeDir(\""<<path<<"\"): {\n";
#endif
  std::string piece;
  size_t pos=std::string::npos;
  if (path[path.length()-1]!='/') path+="/";
  while ((pos=path.find("/"))!=std::string::npos) {
    piece+=path.substr(0,pos)+std::string("/");
    path=path.substr(pos+1);
#ifdef DEBUG__Shell_Tools
    std::cout<<"   Making directory '"<<piece<<"'\n";
#endif
    if (mkdir(piece.c_str(),mode)!=0) {
      if (errno==EEXIST) {
	struct dirent **namelist;
  	if (scandir(piece.c_str(),&namelist,NULL,NULL)>=0) continue;
	else {
#ifdef DEBUG__Shell_Tools
	  std::cout<<"   File exists but is not a directory.\n"
		   <<"   Abort.\n}"<<std::endl;
#endif
	  return false;
	}
      }
#ifdef DEBUG__Shell_Tools
      std::cout<<"   Failed with error'"<<errno<<"'.\n"
	       <<"   Abort.\n}"<<std::endl;
#endif
      return false;
    }
  }
#ifdef DEBUG__Shell_Tools
  std::cout<<"}"<<std::endl;
#endif
  return true;
}
