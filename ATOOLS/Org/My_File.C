#include "ATOOLS/Org/My_File.H"

#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include "ATOOLS/Org/My_MPI.H"

#include <typeinfo>
#include <cstdlib>
#include <sqlite3.h>
#include <string.h>
#define PTS long unsigned int

using namespace ATOOLS;

std::ostream &ATOOLS::operator<<(std::ostream &ostr,const fom::code &code)
{
  switch (code) {
  case fom::temporary: return ostr<<"temporary";
  case fom::permanent: return ostr<<"permanent";
  case fom::error:     return ostr<<"error";
  case fom::unknown:   return ostr<<"unknown";
  }
  return ostr;
}
	
namespace ATOOLS {

  typedef std::map<std::string,std::pair
		   <sqlite3*,std::string> > DataBase_Map;

  static DataBase_Map  s_databases;

  typedef std::map<std::string,sqlite3*> SQLDB_Map;
  typedef std::map<sqlite3*,sqlite3_stmt*> SQLS_Map;
  typedef std::pair<std::string,sqlite3*> DB_Ref;

  SQLDB_Map s_sqldbs;
  SQLS_Map s_getfile;
	
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

int ListFiles(void *data,int argc,char **argv,char **name)
{
  if (argc!=1 || strcmp(name[0],"file")) return 1;
  msg_IODebugging()<<"  '"<<argv[0]<<"' -> '"
		   <<((DB_Ref*)data)->first+argv[0]<<"'\n";
  s_databases[((DB_Ref*)data)->first+argv[0]]=
    std::pair<sqlite3*,std::string>
    (((DB_Ref*)data)->second,((DB_Ref*)data)->first);
  return 0;
}

void PrepareStatements(sqlite3 *db)
{
  char sqlget[100];
  sprintf(sqlget,"select content from path where file = ?1");
  sqlite3_stmt *stmt=NULL;
  int rc=sqlite3_prepare_v2(db,sqlget,100,&stmt,NULL);
  if(rc!=SQLITE_OK)
    msg_IODebugging()<<METHOD<<"(): '"<<db<<"' returns '"
		     <<sqlite3_errmsg(db)<<"'."<<std::endl;
  s_getfile[db]=stmt;
}

template <class FileType>
bool My_File<FileType>::OpenDB(std::string file)
{
  DB_Ref dbref(file,NULL);
  while (file.length() && file[file.length()-1]=='/')
    file.erase(file.length()-1,1);
  file+=".db";
  if (s_sqldbs.find(file)!=s_sqldbs.end()) return true;
  sqlite3 *db=NULL;
  int dummy;
#ifdef USING__MPI
  if (MPI::COMM_WORLD.Get_rank()) {
    MPI::COMM_WORLD.Bcast(&dummy,1,MPI::INT,0);
  }
  else {
#endif
  if (FileExists(file)) {
#ifdef USING__MPI
    MPI::COMM_WORLD.Bcast(&dummy,1,MPI::INT,0);
#endif
  }
  else {
    int res=sqlite3_open(file.c_str(),&db);
    if (res!=SQLITE_OK) {
      msg_IODebugging()<<METHOD<<"(): '"<<file<<"' returns '"
		       <<sqlite3_errmsg(db)<<"'."<<std::endl;
    }
    char *zErrMsg=0;
    std::string sql="create table path(file,content);";
    sql+="create index idx_path on path(file); begin";
    msg_IODebugging()<<METHOD<<"(\""<<file<<"\"): Creating table.\n";
    int rc=sqlite3_exec(db,sql.c_str(),NULL,NULL,&zErrMsg);
    if(rc!=SQLITE_OK) {
      msg_IODebugging()<<METHOD<<"(): '"<<file
		     <<"' returns '"<<zErrMsg<<"'."<<std::endl;
      sqlite3_free(zErrMsg);
      sqlite3_close(db);
      return false;
    }
    s_sqldbs[file]=db;
    PrepareStatements(db);
#ifdef USING__MPI
    MPI::COMM_WORLD.Bcast(&dummy,1,MPI::INT,0);
#endif
    return true;
  }
#ifdef USING__MPI
  }
#endif
  int res=sqlite3_open(file.c_str(),&db);
  if (res!=SQLITE_OK) {
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
    msg_IODebugging()<<METHOD<<"(): '"<<file
		   <<"' returns '"<<zErrMsg<<"'."<<std::endl;
    sqlite3_free(zErrMsg);
    sqlite3_close(db);
    return false;
  }
  msg_IODebugging()<<"}\n";
  s_sqldbs[file]=db;
  PrepareStatements(db);
  return true;
}

template <class FileType> bool
My_File<FileType>::ExecDB(std::string file,const std::string &cmd)
{
  while (file.length() && file[file.length()-1]=='/')
    file.erase(file.length()-1,1);
  file+=".db";
  SQLDB_Map::iterator dbit(s_sqldbs.find(file));
  if (dbit==s_sqldbs.end()) return true;
  msg_IODebugging()<<METHOD<<"("<<file<<"): Executing '"<<cmd<<"'.\n";
  char *zErrMsg=0, *sql = new char[cmd.length()+1];
  strcpy(sql,cmd.c_str());
  int rc=sqlite3_exec(dbit->second,sql,NULL,NULL,&zErrMsg);
  delete [] sql;
  if(rc!=SQLITE_OK) {
    msg_IODebugging()<<METHOD<<"(): '"<<file
		   <<"' returns '"<<zErrMsg<<"'."<<std::endl;
    sqlite3_free(zErrMsg);
    return false;
  }
  return true; 
}

void FinalizeStatements(sqlite3 *db)
{
  int rc=sqlite3_finalize(s_getfile[db]);
  if(rc!=SQLITE_OK)
    msg_IODebugging()<<METHOD<<"(): '"<<db<<"' returns '"
		     <<sqlite3_errmsg(db)<<"'."<<std::endl;
  s_getfile.erase(s_getfile.find(db));
}

template <class FileType>
bool My_File<FileType>::CloseDB(std::string file)
{
  while (file.length() && file[file.length()-1]=='/')
    file.erase(file.length()-1,1);
  file+=".db";
  SQLDB_Map::iterator dbit(s_sqldbs.find(file));
  if (dbit==s_sqldbs.end()) return true;
  msg_IODebugging()<<METHOD<<"("<<file
		 <<"): Closing '"<<dbit->second<<"'.";
  FinalizeStatements(dbit->second);
  int res=sqlite3_close(dbit->second);
  if (res!=SQLITE_OK)
    msg_Error()<<METHOD<<"(): DB '"<<file<<"' returns '"
	       <<sqlite3_errmsg(dbit->second)<<"'."<<std::endl;
  for (DataBase_Map::iterator it(s_databases.begin());
       it!=s_databases.end();)
    if (it->second.first!=dbit->second) ++it;
    else {
      s_databases.erase(it);
      it=s_databases.begin();
    }
  s_sqldbs.erase(dbit);
  return res==SQLITE_OK;
}

template <class FileType>
My_File<FileType>::My_File(const std::string &path,
			   const std::string &file): 
  m_path(path), m_file(file), 
  p_file(NULL), m_mode(fom::permanent) {}

template <class FileType>
My_File<FileType>::~My_File() 
{
  Close();
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
  DataBase_Map::const_iterator sit(s_databases.find(name));
  if (sit!=s_databases.end()) return true;
  return false;
}

template <class FileType> bool
My_File<FileType>::CopyInDB(std::string oldfile, std::string newfile)
{
  DataBase_Map::const_iterator sit(s_databases.find(oldfile));
  if (sit!=s_databases.end()) {
    sqlite3 *db=sit->second.first;
    msg_IODebugging()<<METHOD<<"(): '"<<oldfile<<"' found in '"<<db<<"'";
    oldfile.erase(0,sit->second.second.length());
    newfile.erase(0,sit->second.second.length());
    char *zErrMsg=0;
    std::string sql = "insert into path select '"
      +newfile+"',content from path where file='"+oldfile+"'";
    int rc=sqlite3_exec(db,sql.c_str(),NULL,NULL,&zErrMsg);
    if(rc!=SQLITE_OK) {
      msg_Error()<<METHOD<<"(): '"<<db<<"' returns '"
		     <<zErrMsg<<"'."<<std::endl;
      sqlite3_free(zErrMsg);
    }
    return true;
  }
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
  p_file = new File_Type();
  std::ifstream *is=dynamic_cast<std::ifstream*>(&*p_file);
  std::ofstream *os=dynamic_cast<std::ofstream*>(&*p_file);
  DataBase_Map::const_iterator sit(s_databases.find(m_path+m_file));
  if (is && sit!=s_databases.end()) {
    sqlite3 *db=sit->second.first;
    msg_IODebugging()<<METHOD<<"(): '"<<m_path+m_file
		     <<"' found in '"<<db<<"' {\n";
    p_stream = new MyStrStream();
    std::string fn(m_path+m_file);
    fn.erase(0,sit->second.second.length());
    sqlite3_stmt *stmt(s_getfile[db]);
    sqlite3_bind_text(stmt,1,fn.c_str(),-1,SQLITE_TRANSIENT);
    int rc=sqlite3_step(stmt);
    if (rc==SQLITE_ROW) {
      msg_IODebugging()<<sqlite3_column_text(stmt,0)<<"\n";
      (*p_stream)<<sqlite3_column_text(stmt,0)<<"\n";
    }
    else if (rc==SQLITE_DONE) {
      msg_Error()<<METHOD<<"(): No file content for '"
		 <<fn<<"'"<<std::endl;
    }
    else {
      msg_Error()<<METHOD<<"(): '"<<db<<"' returns '"
		 <<sqlite3_errmsg(db)<<"'."<<std::endl;
    }
    sqlite3_clear_bindings(stmt);
    sqlite3_reset(stmt);
    msg_IODebugging()<<"}\n";
    p_file->copyfmt(*p_stream);
    p_file->clear(p_stream->rdstate());
    is->std::ios::rdbuf(p_stream->rdbuf());
    is->seekg(0);
    return true;
  }
  if (os) {
    if (sit!=s_databases.end()) {
      p_stream = new MyStrStream();
      os->std::ios::rdbuf(p_stream->rdbuf());
      os->seekp(0);
      return true;
    }
    std::string fn(m_path+m_file);
    for (SQLDB_Map::const_iterator sit(s_sqldbs.begin());
	 sit!=s_sqldbs.end();++sit) {
      std::string tag(sit->first);
      tag.replace(tag.length()-3,3,"/");
      if (fn.find(tag)==0) {
	sqlite3 *db=sit->second;
	std::string fn(m_path+m_file);
	msg_IODebugging()<<METHOD<<"(): '"<<fn
			 <<"' added to '"<<db<<"'.\n";
	fn.erase(0,tag.length());
	char *zErrMsg=0, *sql = new char[100+fn.length()];
	sprintf(sql,"insert into path values('%s','')",fn.c_str());
	int rc=sqlite3_exec(db,sql,NULL,NULL,&zErrMsg);
	delete [] sql;
	if(rc!=SQLITE_OK) {
	  msg_Error()<<METHOD<<"(): '"<<db<<"' returns '"
		     <<zErrMsg<<"'."<<std::endl;
	  sqlite3_free(zErrMsg);
	}
	s_databases[m_path+m_file]=
	  std::pair<sqlite3*,std::string>(db,tag);
	p_stream = new MyStrStream();
	os->std::ios::rdbuf(p_stream->rdbuf());
	os->seekp(0);
	return true;
      }
    }
  }
  p_file->open((m_path+m_file).c_str());
  return p_file->good();
}

template <class FileType>
bool My_File<FileType>::Close()
{
  if (p_file==NULL) return false;
  std::ofstream *os=dynamic_cast<std::ofstream*>(&*p_file);
  if (os) {
    DataBase_Map::const_iterator sit(s_databases.find(m_path+m_file));
    if (sit!=s_databases.end()) {
      sqlite3 *db=sit->second.first;
      std::string fn(m_path+m_file), fc(p_stream->str());
      msg_IODebugging()<<METHOD<<"(): Write '"<<fn
		       <<"' to '"<<db<<"' {\n"<<fc;
      fn.erase(0,sit->second.second.length());
      if (fc[fc.length()-1]=='\n') fc.erase(fc.length()-1,1);
      for (size_t pos=fc.find("'");pos!=std::string::npos;
	   pos=fc.find("'",pos+2)) fc.replace(pos,1,"''");
      char *zErrMsg=0, *sql = new char[100+fn.length()+fc.length()];
      sprintf(sql,"update path set content = '%s' where file = '%s'",
	      fc.c_str(),fn.c_str());
      int rc=sqlite3_exec(db,sql,NULL,NULL,&zErrMsg);
      delete [] sql;
      if(rc!=SQLITE_OK) {
      	msg_Error()<<METHOD<<"(): '"<<db<<"' returns '"
      		   <<zErrMsg<<"'."<<std::endl;
      	sqlite3_free(zErrMsg);
      }
      msg_IODebugging()<<"}\n";
    }
  }
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

namespace ATOOLS {

  template class My_In_File;
  template class My_Out_File;

}

