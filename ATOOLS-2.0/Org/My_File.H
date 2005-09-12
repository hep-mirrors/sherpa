#ifndef My_File_H
#define My_File_H

#include "Smart_Pointer.H"

#include <fstream>
#include <string>

namespace ATOOLS {

  const std::string nullstring("");

  struct fom {
    
    enum code {
      temporary =  1,
      permanent =  2,
      error     =  3,
      unknown   = 99
    };

  };// end of struct fom
    
  std::ostream &operator<<(std::ostream &ostr,const fom::code &code);

  template <class FileType>
  class My_File {
  public:

    typedef FileType File_Type;

  protected:

    std::string m_path, m_file;

    SP(File_Type) p_file;
    
    fom::code m_mode;

  public:

    // constructor
    inline My_File(): p_file(NULL), m_mode(fom::permanent) {}
    // destructor
    inline ~My_File() { Close(); }

    // member functions
    bool Open();
    bool Close();

    //inline functions
    inline void SetPath(const std::string &path) { m_path=path; }
    inline void SetFile(const std::string &file) { m_file=file; }

    inline void SetMode(const fom::code &mode) { m_mode=mode; }

    inline const std::string &Path() const { return m_path; }
    inline const std::string &File() const { return m_file; }

    inline const fom::code &Mode() const { return m_mode; }

    inline File_Type *operator()() const { return --p_file; }
    inline File_Type *operator->() const { return --p_file; }
    inline File_Type &operator*() const  { return *p_file;  }

  };// end of class My_File

  template<class FileType>
  std::ostream &operator<<(std::ostream &ostr,
			   const My_File<FileType> &file);

}// end of namespace ATOOLS

#define My_In_File ATOOLS::My_File<std::ifstream>
#define My_Out_File ATOOLS::My_File<std::ofstream>

#endif
