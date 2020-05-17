#include "ATOOLS/Org/My_File.H"

#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Org/libzippp.h"

#include <sys/stat.h>
#include <typeinfo>
#include <cstdlib>
#include <string.h>
#include <algorithm>
#define PTS long unsigned int

using namespace libzippp;

namespace ATOOLS {

  typedef std::pair<ZipArchive*,std::vector<std::string> > ZipArchive_Ref;
  typedef std::map<std::string,ZipArchive_Ref> ZipArchive_Map;
  typedef std::pair<std::string,int> Zip_Entry;
  typedef std::map<std::string,Zip_Entry> ZipEntry_Map;

  ZipArchive_Map s_ziparchives;
  ZipEntry_Map s_zipfiles;

}

using namespace ATOOLS;

std::ostream &ATOOLS::operator<<(std::ostream &ostr,const fom::code &code)
{
  switch (code) {
  case fom::temporary: return ostr<<"temporary";
  case fom::permanent: return ostr<<"permanent";
  case fom::unknown:   return ostr<<"unknown";
  }
  return ostr;
}
	
namespace ATOOLS {

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

template <class FileType>
bool My_File<FileType>::OpenDB(std::string file)
{
  std::string path(file);
#ifdef USING__MPI
  if (mpi->Rank()) {
    s_ziparchives[path]=ZipArchive_Ref(NULL,std::vector<std::string>());
    int size;
    mpi->Bcast(&size,1,MPI_INT);
    for (int i=0;i<size;++i) {
      int length;
      mpi->Bcast(&length,1,MPI_INT);
      char *message = new char[length+1];
      mpi->Bcast(message,length+1,MPI_CHAR);
      std::string name, content;
      for (int p=0;p<length;++p)
	if (message[p]=='\n') {
	  name=std::string(message,p);
	  content=std::string(&message[p+1],length-p-1);
	  break;
	}
      s_ziparchives[path].second.push_back(name);
      s_zipfiles[name]=Zip_Entry(content,0);
      delete [] message;
    }
    return true;
  }
#endif
  while (file.length() && file[file.length()-1]=='/')
    file.erase(file.length()-1,1);
  file+=".zip";
  ZipArchive *zf(new ZipArchive(file));
  s_ziparchives[path]=ZipArchive_Ref(zf,std::vector<std::string>());
  int res=zf->open(ZipArchive::WRITE);
  const std::vector<ZipEntry> &entries=zf->getEntries();
  int size=entries.size();
#ifdef USING__MPI
  mpi->Bcast(&size,1,MPI_INT);
#endif
  for(std::vector<ZipEntry>::const_iterator
	it=entries.begin();it!=entries.end();++it) {
    std::string name=path+it->getName();
    std::string content=it->readAsText();
    s_ziparchives[path].second.push_back(name);
    s_zipfiles[name]=Zip_Entry(content,0);
#ifdef USING__MPI
    int length(name.length()+content.length()+1);
    content=name+'\n'+content;
    mpi->Bcast(&length,1,MPI_INT);
    mpi->Bcast(&content[0],length+1,MPI_CHAR);
#endif
  }
  return true;
}

template <class FileType>
bool My_File<FileType>::CloseDB(std::string file,int mode)
{
  if (mpi->Rank()) {
    int success;
    mpi->Bcast(&success,1,MPI_INT);
    return success;
  }
  std::string path(file);
  while (file.length() && file[file.length()-1]=='/')
    file.erase(file.length()-1,1);
  file+=".zip";
  ZipArchive_Map::iterator ait(s_ziparchives.find(path));
  if (ait==s_ziparchives.end()) {
    int success(false);
    mpi->Bcast(&success,1,MPI_INT);
    return success;
  }
  ZipArchive *zf(ait->second.first);
  const std::vector<std::string> &files(ait->second.second);
  for (size_t i(0);i<files.size();++i) {
    ZipEntry_Map::iterator zit(s_zipfiles.find(files[i]));
    if (zf) {
      std::string fn(files[i]);
      fn.erase(0,path.length());
      if (zit->second.second<0) zf->deleteEntry(fn);
      if (zit->second.second>0) {
	char *tmp = new char[zit->second.first.length()+1];
	strcpy(tmp,zit->second.first.c_str());
	zf->addData(fn,tmp,strlen(tmp));
      }
      zit->second.second=0;
    }
    if (mode) s_zipfiles.erase(zit);
  }
  if (mode) s_ziparchives.erase(ait);
  if (zf) {
    zf->close();
    if (mode) delete zf;
    else zf->open(ZipArchive::WRITE);
  }
  int success(true);
  mpi->Bcast(&success,1,MPI_INT);
  return success;
}

template <class FileType>
My_File<FileType>::My_File(const std::string &path,
			   const std::string &file): 
  m_path(path), m_file(file), 
  m_mode(fom::permanent) {}

template <class FileType>
My_File<FileType>::~My_File() 
{
  Close();
}

template <class FileType>
FileType *My_File<FileType>::operator()() const 
{ 
  return p_file.get();
}

template <class FileType>
FileType *My_File<FileType>::operator->() const 
{ 
  return p_file.get();
}

template <class FileType>
FileType &My_File<FileType>::operator*() const  
{ 
  return *p_file;  
}

template <class FileType> bool 
My_File<FileType>::FileInDB(const std::string &name)
{
  return s_zipfiles.find(name)!=s_zipfiles.end();
}

template <class FileType> bool
My_File<FileType>::CopyInDB(std::string oldfile, std::string newfile)
{
  if (!FileExists(oldfile)) return false;
  My_In_File infile(oldfile);
  if (!infile.Open()) return false;
  My_Out_File outfile(newfile);
  if (!outfile.Open()) return false;
  *outfile<<infile.FileContent();
  return true;
}

template <class FileType>
bool My_File<FileType>::Open() 
{ 
  if (m_path=="" && m_file=="") {
    p_file = std::make_shared<File_Type>();
    return false;
  }
  Close();
  p_file = std::make_shared<File_Type>();
  std::ifstream *is=dynamic_cast<std::ifstream*>(p_file.get());
  std::ofstream *os=dynamic_cast<std::ofstream*>(p_file.get());
  if (is) {
    p_stream = std::make_shared<MyStrStream>();
    ZipEntry_Map::const_iterator zit=
      s_zipfiles.find(m_path+m_file);
    if (zit!=s_zipfiles.end()) {
      (*p_stream)<<zit->second.first;
    }
    else {
#ifdef USING__MPI
    if (mpi->Rank()) {
      int fsize;
      mpi->Bcast(&fsize,1,MPI_INT);
      if (fsize<0) return false;
      char *content = new char[fsize+1];
      mpi->Bcast(content,fsize+1,MPI_CHAR);
      (*p_stream)<<content<<"\n";
      delete [] content;
    }
    else {
#endif
      std::ifstream infile((m_path+m_file).c_str());
      int fsize(infile.good()?1:-1);
      if (fsize<0) {
#ifdef USING__MPI
	mpi->Bcast(&fsize,1,MPI_INT);
#endif
	return false;
      }
      msg_IODebugging()<<infile.rdbuf()<<"\n";
      (*p_stream)<<infile.rdbuf();
#ifdef USING__MPI
      std::string content(p_stream->str());
      fsize=content.length();
      mpi->Bcast(&fsize,1,MPI_INT);
      mpi->Bcast(&content[0],fsize+1,MPI_CHAR);
    }
#endif
    }
    msg_IODebugging()<<"}\n";
    p_file->copyfmt(*p_stream);
    p_file->clear(p_stream->rdstate());
    is->std::ios::rdbuf(p_stream->rdbuf());
    is->seekg(0);
    return true;
  }
  if (os) {
    p_stream = std::make_shared<MyStrStream>();
    os->std::ios::rdbuf(p_stream->rdbuf());
    os->seekp(0);
    return true;
  }
  return false;
}

template <class FileType>
bool My_File<FileType>::Close()
{
  if (p_file == nullptr)
    return false;
  auto os = dynamic_cast<std::ofstream*>(p_file.get());
  if (os) {
    bool indb(false);
    for (ZipArchive_Map::iterator zit(s_ziparchives.begin());
	 zit!=s_ziparchives.end();++zit)
      if ((m_path+m_file).find(zit->first)==0) {
	ZipEntry_Map::iterator fit(s_zipfiles.find(m_path+m_file));
	if (fit!=s_zipfiles.end()) fit->second=Zip_Entry(p_stream->str(),2);
	else {
	  s_zipfiles[m_path+m_file]=Zip_Entry(p_stream->str(),2);
	  zit->second.second.push_back(m_path+m_file);
	}
	indb=true;
	break;
      }
#ifdef USING__MPI
    if (mpi->Rank()==0)
#endif
    if (!indb) {
      std::ofstream file(m_path+m_file);
      file<<p_stream->str();
    }
  }
  p_file->close();
  p_stream.reset();
  p_file.reset();
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

