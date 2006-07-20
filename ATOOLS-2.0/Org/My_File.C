#include "My_File.H"

#include "Exception.H"

#include <typeinfo>

INSTANTIATE_SMART_POINTER(std::ifstream)
INSTANTIATE_SMART_POINTER(std::ofstream)

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
My_File<FileType>::My_File(): 
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
  p_file = new File_Type();
  p_file->open((m_path+m_file).c_str());
  return p_file->good();

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

template class My_File<std::ifstream>;
template class My_File<std::ofstream>;

