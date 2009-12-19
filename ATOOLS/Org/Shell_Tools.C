#include "ATOOLS/Org/Shell_Tools.H"

#include "ATOOLS/Org/My_File.H"

#include <sys/stat.h>
#include <errno.h>
#include <dirent.h>
#include <cstdlib>
#if __GNUC__
#include <cxxabi.h>
#endif

#ifdef DEBUG__Shell_Tools
#include <iostream>
#endif

using namespace ATOOLS;

bool ATOOLS::MakeDir(std::string path,const bool create_tree,
		     const mode_t mode)
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
	if (n>=0) {
	  free(namelist);
	  continue;
	}
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

bool ATOOLS::ChMod(const std::string &file,const mode_t mode)
{
  if (!FileExists(file)) return false;
  if (chmod(file.c_str(),mode)!=0) {
#ifdef DEBUG__Shell_Tools
    std::cout<<METHOD<<"(): Error "<<errno<<" in setting mode "
	     <<mode<<"on '"<<file<<"'."<<std::endl;
#endif
    return false;
  }
  return true;
}

bool ATOOLS::CopyFile(const std::string &oldname,const std::string &newname)
{
  if (!FileExists(oldname)) return false;
  My_In_File oldfile("",oldname);
  if (!oldfile.Open()) return false;
  My_Out_File newfile("",newname);
  if (!newfile.Open()) return false;
  (*newfile)<<oldfile->rdbuf();
  struct stat fst;
  stat(oldname.c_str(),&fst);
  chmod(newname.c_str(),fst.st_mode);
  return true;
}

bool ATOOLS::MoveFile(const std::string &oldname,const std::string &newname)
{
  if (!CopyFile(oldname,newname)) return false;
  return !remove(oldname.c_str());
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

bool ATOOLS::FileExists(const std::string &file)
{
  struct stat fst;
  if (stat(file.c_str(),&fst)!=-1)
    return (fst.st_mode&S_IFMT)==S_IFREG;
  return false;
}

std::string ATOOLS::Demangle(const std::string &name)
{
#if __GNUC__
  int s;
  size_t len(name.length());
  char *res(abi::__cxa_demangle(name.c_str(),0,&len,&s));
  return s==0?std::string(res):name;
#else
  return name;
#endif
}
