#include "My_File.H"

#include "Exception.H"

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

template <>
std::ostream &operator<<<std::ifstream>(std::ostream &ostr,
					const My_File<std::ifstream> &file)
{
  return ostr<<"("<<(&*file)<<") [input] { m_path = "<<file.Path()
	     <<", m_file = "<<file.File()
	     <<", m_mode = "<<file.Mode()<<" }";
}

template <>
std::ostream &operator<<<std::ofstream>(std::ostream &ostr,
					const My_File<std::ofstream> &file)
{
  return ostr<<"("<<(&*file)<<") [output] { m_path = "<<file.Path()
	     <<", m_file = "<<file.File()
	     <<", m_mode = "<<file.Mode()<<" }";
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
  delete p_file;
  p_file=NULL;
  return true;
}

template class My_File<std::ifstream>;
template class My_File<std::ofstream>;

