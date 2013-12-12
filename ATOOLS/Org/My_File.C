#include "ATOOLS/Org/My_File.H"

#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"

#include <typeinfo>
#include <cstdlib>
#ifdef USING__SQLITE
#include <sqlite3.h>
#include <string.h>
#define PTS long unsigned int
#endif

using namespace ATOOLS;

std::ostream &ATOOLS::operator<<(std::ostream &ostr,const fom::code &code)
{
  switch (code) {
  case fom::temporary: return ostr<<"temporary";
  case fom::permanent: return ostr<<"permanent";
  case fom::error:     return ostr<<"error";
  case fom::nosearch:  return ostr<<"nosearch";
  case fom::unknown:   return ostr<<"unknown";
  }
  return ostr;
}
	
namespace ATOOLS {

#ifdef USING__SQLITE
  std::map<std::string,sqlite3*> s_sqldbs;
  typedef std::pair<std::string,sqlite3*> DB_Ref;
#endif
	
  template <> std::ostream &
  operator<<<std::ifstream>(std::ostream &ostr,
			    const My_File<std::ifstream> &file)
  {
    return ostr<<"("<<(&*file)<<") [input] { m_path = "<<file.Path()
	       <<", m_file = "<<file.File()
	       <<", m_mode = "<<file.Mode()<<" }";
  }

  template <> std::ostream &
  operator<<<std::ofstream>(std::ostream &ostr,
			    const My_File<std::ofstream> &file)
  {
    return ostr<<"("<<(&*file)<<") [output] { m_path = "<<file.Path()
	       <<", m_file = "<<file.File()
	       <<", m_mode = "<<file.Mode()<<" }";
  }

}

template <class FileType> int My_File<FileType>::
ListFiles(void *data,int argc,char **argv,char **name)
{
#ifdef USING__SQLITE
  if (argc!=1 || strcmp(name[0],"file")) return 1;
  msg_IODebugging()<<"  '"<<argv[0]<<"' -> '"
		   <<((DB_Ref*)data)->first+argv[0]<<"'\n";
  s_databases[((DB_Ref*)data)->first+argv[0]]=
    std::pair<void*,std::string>
    (((DB_Ref*)data)->second,((DB_Ref*)data)->first);
#endif
  return 0;
}

template <class FileType> int My_File<FileType>::
GetFile(void *data,int argc,char **argv,char **name)
{
#ifdef USING__SQLITE
  if (argc!=1 || strcmp(name[0],"content")) return 1;
  msg_IODebugging()<<argv[0]<<"\n";
  (*(MyStrStream*)((PTS)data))<<argv[0]<<"\n";
#endif
  return 0;
}

template <class FileType>
bool My_File<FileType>::OpenDB(std::string file)
{
#ifdef USING__SQLITE
  DB_Ref dbref(file,NULL);
  while (file.length() && file[file.length()-1]=='/')
    file.erase(file.length()-1,1);
  file+=".db";
  if (s_sqldbs.find(file)!=s_sqldbs.end()) return true;
  sqlite3 *db=NULL;
  for (size_t i(0);i<s_searchpaths.size();++i) {
    if (!FileExists(s_searchpaths[i]+"/"+file)) continue;
    int res=sqlite3_open((s_searchpaths[i]+"/"+file).c_str(),&db);
    if (res==SQLITE_OK) {
      msg_IODebugging()<<METHOD<<"(): '"<<file<<"' found in '"
		       <<s_searchpaths[i]<<"'."<<std::endl;
      break;
    }
    msg_IODebugging()<<METHOD<<"(): '"<<file<<"' returns '"
		   <<sqlite3_errmsg(db)<<"'."<<std::endl;
  }
  if (db==NULL) {
    msg_IODebugging()<<METHOD<<"(): '"<<file
		     <<"' not found."<<std::endl;
    return false;
  }
  dbref.second=db;
  char sql[100], *zErrMsg=0;
  strcpy(sql,"select file from path");
  msg_IODebugging()<<METHOD<<"(\""<<file<<"\"): {\n";
  int rc=sqlite3_exec(db,sql,ListFiles,(void*)&dbref,&zErrMsg);
  if(rc!=SQLITE_OK) {
    msg_Debugging()<<METHOD<<"(): '"<<file
		   <<"' returns '"<<zErrMsg<<"'."<<std::endl;
    sqlite3_free(zErrMsg);
    sqlite3_close(db);
    return false;
  }
  msg_IODebugging()<<"}\n";
  s_sqldbs[file]=db;
  return true;
#else
  return false;
#endif
}

template <class FileType>
bool My_File<FileType>::CloseDB(std::string file)
{
#ifdef USING__SQLITE
  while (file.length() && file[file.length()-1]=='/')
    file.erase(file.length()-1,1);
  file+=".db";
  std::map<std::string,sqlite3*>::iterator
    dbit(s_sqldbs.find(file));
  if (dbit==s_sqldbs.end()) return true;
  msg_Debugging()<<METHOD<<"("<<file
		 <<"): Closing '"<<dbit->second<<"'.";
  int res=sqlite3_close(dbit->second);
  if (res!=SQLITE_OK)
    msg_Error()<<METHOD<<"(): DB '"<<file
	       <<"' returns "<<res<<"."<<std::endl;
  for (DataBase_Map::iterator it(s_databases.begin());
       it!=s_databases.end();)
    if (it->second.first!=dbit->second) ++it;
    else {
      s_databases.erase(it);
      it=s_databases.begin();
    }
  s_sqldbs.erase(dbit);
  return res==SQLITE_OK;
#else
  return false;
#endif
}

template <class FileType>
My_File<FileType>::My_File(const std::string &path,
			   const std::string &file): 
  m_path(path), m_file(file), 
  p_file(NULL), m_mode(fom::permanent) {}

template <class FileType>
My_File<FileType>::~My_File() 
{ 
}

template <class FileType>
FileType *My_File<FileType>::operator()() const 
{ 
  return --p_file; 
}

template <class FileType>
FileType *My_File<FileType>::operator->() const 
{ 
  return --p_file; 
}

template <class FileType>
FileType &My_File<FileType>::operator*() const  
{ 
  return *p_file;  
}

template <class FileType> bool 
My_File<FileType>::FileInDB(const std::string &name)
{
#ifdef USING__SQLITE
  DataBase_Map::const_iterator sit(s_databases.find(name));
  if (sit!=s_databases.end()) return true;
#endif
  return false;
}

template <class FileType>
bool My_File<FileType>::Open() 
{ 
  if (m_path=="" && m_file=="") {
    p_file = new File_Type();
    return false;
  }
  Close();
  if (m_mode&fom::nosearch) {
    p_file = new File_Type();
    p_file->open((m_path+m_file).c_str());
    return p_file->good();
  }
#ifdef USING__SQLITE
  DataBase_Map::const_iterator sit(s_databases.find(m_path+m_file));
  if (sit!=s_databases.end()) {
    sqlite3 *db=(sqlite3*)sit->second.first;
    msg_IODebugging()<<METHOD<<"(): '"<<m_path+m_file
		     <<"' found in '"<<db<<"' {\n";
    p_stream = new MyStrStream();
    std::string fn(m_path+m_file);
    fn.erase(0,sit->second.second.length());
    char *zErrMsg=0, *sql = new char[100+fn.length()];
    sprintf(sql,"select content from path where file = \"%s\"",fn.c_str());
    int rc=sqlite3_exec(db,sql,GetFile,(void*)&*p_stream,&zErrMsg);
    if(rc!=SQLITE_OK) {
      msg_Error()<<METHOD<<"(): '"<<db<<"' returns '"
		     <<zErrMsg<<"'."<<std::endl;
      sqlite3_free(zErrMsg);
    }
    msg_IODebugging()<<"}\n";
    p_file = new File_Type();
    p_file->copyfmt(*p_stream);
    p_file->clear(p_stream->rdstate());
    std::ifstream *is=dynamic_cast<std::ifstream*>(&*p_file);
    if (is) {
      is->basic_ios<char>::rdbuf(p_stream->rdbuf());
      is->seekg(0);
    }
    return true;
  }
#endif
  String_Map::const_iterator fit(s_filelocations.find(m_path+m_file));
  if (fit!=s_filelocations.end()) m_path=fit->second+"/"+m_path;
  if ((m_path+m_file)[0]=='/') {
    msg_IODebugging()<<METHOD<<"(): Relocated '"<<m_file<<"' at '"
		   <<m_path<<"'."<<std::endl;  
    p_file = new File_Type();
    p_file->open((m_path+m_file).c_str());
    return p_file->good();
  }
  else {
    for (size_t i(0);i<s_searchpaths.size();++i) {
      p_file = new File_Type();
      p_file->open((s_searchpaths[i]+"/"+m_path+m_file).c_str());
      if (p_file->good()) {
	if (i>0 || msg_LevelIsDebugging()) {
	  bool output(true);
	  for (size_t j(0);j<s_nocomplains.size();++j)
	    if ((m_path+m_file).find(s_nocomplains[j])!=
		std::string::npos) output=false;
	  if (output) 
	    msg_Out()<<METHOD<<"(): Located '"<<m_file<<"' at '"
		     <<s_searchpaths[i]<<"/"<<m_path<<"'."<<std::endl;
	}
	s_filelocations[m_path+m_file]=s_searchpaths[i];
	m_path=s_searchpaths[i]+"/"+m_path;
	return true;
      }
    }
  }
  return false;
}

template <class FileType>
bool My_File<FileType>::Close()
{
  if (p_file==NULL) return false;
  p_file->close();
  p_stream=NULL;
  p_file=NULL;
  return true;
}

template <class FileType>
void My_File<FileType>::SetPath(const std::string &path) 
{
  m_path=path; 
}

template <class FileType>
void My_File<FileType>::SetFile(const std::string &file) 
{ 
  m_file=file; 
}

template <class FileType>
void My_File<FileType>::SetMode(const fom::code &mode) 
{
  m_mode=mode; 
}

template <class FileType>
const std::string &My_File<FileType>::Path() const 
{ 
  return m_path; 
}

template <class FileType>
const std::string &My_File<FileType>::File() const 
{ 
  return m_file; 
}

template <class FileType>
const fom::code &My_File<FileType>::Mode() const 
{ 
  return m_mode; 
}

template <class FileType>
void My_File<FileType>::SetSearchPaths(const String_Vector &paths)
{ 
  s_searchpaths=paths;
}

template <class FileType>
void My_File<FileType>::SetNoComplains(const String_Vector &names)
{ 
  s_nocomplains=names;
}

template <class FileType> 
typename My_File<FileType>::String_Vector 
My_File<FileType>::s_searchpaths(1,GetCWD());
template <class FileType> 
typename My_File<FileType>::String_Vector 
My_File<FileType>::s_nocomplains;
template <class FileType> 
typename My_File<FileType>::String_Map 
My_File<FileType>::s_filelocations;
template <class FileType> 
typename My_File<FileType>::DataBase_Map 
My_File<FileType>::s_databases;

namespace ATOOLS {

  template class My_In_File;
  template class My_Out_File;

}

