#ifndef My_File_H
#define My_File_H

#include "Smart_Pointer.H"

#include <fstream>
#include <vector>
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

    typedef std::vector<std::string> String_Vector;

  protected:

    static String_Vector s_searchpaths;

    std::string m_path, m_file;

    SP(File_Type) p_file;
    
    fom::code m_mode;

  public:

    // constructor
    My_File(const std::string &path="",const std::string &file="");
    // destructor
    ~My_File();

    // member functions
    bool Open();
    bool Close();

    //inline functions
    void SetPath(const std::string &path);
    void SetFile(const std::string &file);

    void SetMode(const fom::code &mode);

    const std::string &Path() const;
    const std::string &File() const;

    const fom::code &Mode() const;

    File_Type *operator()() const;
    File_Type *operator->() const;
    File_Type &operator*() const;

    static void SetSearchPaths(const String_Vector &paths);

  };// end of class My_File

  template<class FileType>
  std::ostream &operator<<(std::ostream &ostr,
			   const My_File<FileType> &file);

}// end of namespace ATOOLS

#define My_In_File ATOOLS::My_File<std::ifstream>
#define My_Out_File ATOOLS::My_File<std::ofstream>

#endif
