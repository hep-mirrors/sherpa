#include "ATOOLS/Org/Shell_Tools.H"

#include "ATOOLS/Org/My_File.H"

#include <sys/stat.h>
#include <string.h>
#include <errno.h>
#include <dirent.h>
#include <cstdlib>
#include <unistd.h>
#include <regex.h>
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
    if (DirectoryExists(piece)) continue;
    if (mkdir(piece.c_str(),mode)!=0) {
      if (errno==EEXIST) {
#ifdef DEBUG__Shell_Tools
	std::cout<<"   File exists but is not a directory.\n"
		 <<"   Abort.\n}"<<std::endl;
#endif
	return false;
      }
#ifdef DEBUG__Shell_Tools
      std::cout<<"   Error. "<<strerror(errno)<<".\n"
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
    std::cout<<"   Error. "<<strerror(errno)
	     <<" while setting mode "
	     <<mode<<" on '"<<file<<"'."<<std::endl;
#endif
    return false;
  }
  return true;
}

bool ATOOLS::Copy(const std::string &oldname,
		  const std::string &newname,const bool rec)
{
  struct stat fst;
  if (stat(oldname.c_str(),&fst)==-1) return false;
  if ((fst.st_mode&S_IFMT)==S_IFDIR) {
    if (!MakeDir(newname,fst.st_mode)) return false;
    bool stat=true;
    struct dirent **entries;
    int n(scandir(oldname.c_str(),&entries,NULL,NULL));
    for (int i(0);i<n;++i) {
      if (strcmp(".",entries[i]->d_name)!=0 &&
	  strcmp("..",entries[i]->d_name)!=0 && rec)
	stat&=Copy(oldname+"/"+entries[i]->d_name,
		   newname+"/"+entries[i]->d_name,rec);
      free(entries[i]);
    }
    if (n>=0) free(entries);
    return stat;
  }
  if (!FileExists(oldname)) return false;
  std::ifstream oldfile(oldname.c_str());
  if (!oldfile.good()) return false;
  std::ofstream newfile(newname.c_str());
  if (!newfile.good()) return false;
  newfile<<oldfile.rdbuf();
  return chmod(newname.c_str(),fst.st_mode)==0;
}

bool ATOOLS::Move(const std::string &oldname,
		  const std::string &newname)
{
  if (!Copy(oldname,newname,true)) return false;
  return !Remove(oldname);
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

bool ATOOLS::FileExists(const std::string &file,const int mode)
{
  if (My_In_File::FileInDB(file)) return true;
  if (mode) return false;
  struct stat fst;
  if (stat(file.c_str(),&fst)!=-1)
    return (fst.st_mode&S_IFMT)==S_IFREG;
  return false;
}

bool ATOOLS::DirectoryExists(const std::string &dir)
{
  struct stat fst;
  if (stat(dir.c_str(),&fst)!=-1)
    return (fst.st_mode&S_IFMT)==S_IFDIR;
  return false;
}

bool ATOOLS::Remove(const std::string &file,
		    const bool rec)
{
  struct stat fst;
  if (stat(file.c_str(),&fst)==-1) return false;
  if ((fst.st_mode&S_IFMT)==S_IFDIR) {
    bool stat=true;
    struct dirent **entries;
    int n(scandir(file.c_str(),&entries,NULL,NULL));
    for (int i(0);i<n;++i) {
      if (strcmp(".",entries[i]->d_name)!=0 &&
	  strcmp("..",entries[i]->d_name)!=0 && rec)
	stat&=Remove(file+"/"+entries[i]->d_name,rec);
      free(entries[i]);
    }
    if (n>=0) free(entries);
    if (!stat) return false;
    return rmdir(file.c_str())==0;
  }
  return unlink(file.c_str())==0;
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

std::string ATOOLS::GetCWD()
{
  long int size = pathconf(".",_PC_PATH_MAX);
  char *buf = new char[size];
  char *ptr = getcwd(buf, (size_t)size);
  if (ptr==NULL) abort();
  std::string pwd(buf);
  delete [] buf;
  return pwd;
}

std::vector<std::string> ATOOLS::RegExMatch
(const std::string &str,const std::string &pat,const size_t nm)
{
  std::vector<std::string> res;
  regex_t re;
  if (regcomp(&re,pat.c_str(),REG_EXTENDED)!=0) return res;
  std::vector<regmatch_t> pm(nm);
  int stat=regexec(&re,str.c_str(),nm,&pm.front(),0);
  if (stat==0) {
    for (size_t i=0;i<nm;++i)
      res.push_back(str.substr(pm[i].rm_so,pm[i].rm_eo-pm[i].rm_so));
  }
  regfree(&re);
  return res;
}

std::string ATOOLS::ShortenPathName(std::string path)
{
  while (path.length() && path[path.length()-1]=='/')
    path.erase(path.length()-1);
  for (size_t pos=path.find("//");pos!=std::string::npos;
       pos=path.find("//")) path.erase(pos,1);
  for (size_t pos=path.find("./");pos!=std::string::npos;
       pos=path.find("./")) path.erase(pos,2);
  return path;
}
