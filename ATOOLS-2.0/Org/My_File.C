#include "My_File.H"

#include "Exception.H"

#include <typeinfo>

using namespace ATOOLS;

namespace ATOOLS {
  INSTANTIATE_SMART_POINTER(std::ifstream)
  INSTANTIATE_SMART_POINTER(std::ofstream)
}

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
	
  template <> std::ostream &
  operator<<<std::ifstream>(std::ostream &ostr,
			    const My_File<std::ifstream> &file)
  {
    return ostr<<"("<<(&*file)<<") [input] { m_path = "<<file.Path()
	       <<", m_file = "<<file.File()
	       <<", m_mode = "<<file.Mode()<<" }";
  }

  template <>	std::ostream &
  operator<<<std::ofstream>(std::ostream &ostr,
			    const My_File<std::ofstream> &file)
  {
    return ostr<<"("<<(&*file)<<") [output] { m_path = "<<file.Path()
	       <<", m_file = "<<file.File()
	       <<", m_mode = "<<file.Mode()<<" }";
  }

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

template <class FileType>
bool My_File<FileType>::Open() 
{ 
  if (m_path=="" && m_file=="") return false;
  Close();
  String_Map::const_iterator fit(s_filelocations.find(m_path+m_file));
  if (fit!=s_filelocations.end()) m_path=fit->second;
  if ((m_path+m_file)[0]=='/') {
    msg_Debugging()<<METHOD<<"(): Relocated '"<<m_file<<"' at '"
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
	if (i>0 || msg_LevelIsDebugging()) 
	  msg_Out()<<METHOD<<"(): Located '"<<m_file<<"' at '"
		   <<s_searchpaths[i]<<"/"<<m_path<<"'."<<std::endl;
	m_path=s_filelocations[m_path+m_file]=s_searchpaths[i]+"/"+m_path;
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
  p_file.Delete();
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
typename My_File<FileType>::String_Vector 
My_File<FileType>::s_searchpaths;
template <class FileType> 
typename My_File<FileType>::String_Map 
My_File<FileType>::s_filelocations;

namespace ATOOLS {

  template class My_In_File;
  template class My_Out_File;

}
