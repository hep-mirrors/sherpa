#include "Shell_Tools.H"

#include <sys/stat.h>
#include <errno.h>
#include <dirent.h>
#include <fstream>

#ifdef DEBUG__Shell_Tools
#include <iostream>
#endif

using namespace ATOOLS;

bool ATOOLS::MakeDir(std::string path,const mode_t mode,const bool create_tree)
{
  if (path=="") return false;
#ifdef DEBUG__Shell_Tools
  std::cout<<"ATOOLS::MakeDir(\""<<path<<"\"): {\n";
#endif
  std::string piece;
  size_t pos=std::string::npos;
  if (path[path.length()-1]!='/') path+="/";
  if (!create_tree) return !mkdir(path.c_str(),mode);
  while ((pos=path.find("/"))!=std::string::npos) {
    if (pos==0) {
      piece+=std::string("/");
      path=path.substr(pos+1);      
      pos=path.find("/");
    }
    piece+=path.substr(0,pos)+std::string("/");
    path=path.substr(pos+1);
#ifdef DEBUG__Shell_Tools
    std::cout<<"   Making directory '"<<piece<<"'\n";
#endif
    if (mkdir(piece.c_str(),mode)!=0) {
      if (errno==EEXIST) {
	struct dirent **namelist;
	int n(scandir(piece.c_str(),&namelist,NULL,NULL));
	for (int i(0);i<n;++i) free(namelist[i]);
	free(namelist);
  	if (n>=0) continue;
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

bool ATOOLS::CopyFile(const std::string &oldname,const std::string &newname)
{
  std::ifstream oldfile(oldname.c_str());
  if (oldfile.bad()) return false;
  std::ofstream newfile(newname.c_str());
  if (newfile.bad()) return false;
  newfile<<oldfile.rdbuf();
  return true;
}

std::vector<std::string> 
ATOOLS::EnvironmentVariable(const std::string &name,std::string entry)
{
#ifdef DEBUG__Shell_Tools
  std::cout<<"EnvironmentVariable("<<name<<"): {\n";
#endif
  if (entry.length()==0) {
    char *var=NULL;
    entry=(var=getenv(name.c_str()))==NULL?"":var;
  }
  size_t pos=std::string::npos;
  std::vector<std::string> entries;
  if (entry[entry.length()-1]!=':') entry+=":";
  while ((pos=entry.find(":"))!=std::string::npos) {
    if (pos>0) entries.push_back(entry.substr(0,pos));
    entry=entry.substr(pos+1);
#ifdef DEBUG__Shell_Tools
    if (pos>0) std::cout<<"   Extracted entry '"<<entries.back()<<"'\n";
#endif
  }
#ifdef DEBUG__Shell_Tools
  std::cout<<"}"<<std::endl;
#endif
  return entries;
}
