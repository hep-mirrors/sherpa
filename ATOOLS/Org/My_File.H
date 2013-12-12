#ifndef ATOOLS_Org_My_File_H
#define ATOOLS_Org_My_File_H

#include "ATOOLS/Org/Smart_Pointer.H"
#include "ATOOLS/Org/MyStrStream.H"

#include <fstream>
#include <vector>
#include <string>
#include <map>

namespace ATOOLS {

  const std::string nullstring("");

  struct fom {
    
    enum code {
      temporary =  1,
      permanent =  2,
      error     =  4,
      nosearch  =  8,
      unknown   = 99
    };

  };// end of struct fom
    
  inline fom::code operator|(const fom::code c1,const fom::code c2)
  { return (fom::code)((int)c1|(int)c2); }

  std::ostream &operator<<(std::ostream &ostr,const fom::code &code);

  template <class FileType>
  class My_File {
  public:

    typedef FileType File_Type;

    typedef std::vector<std::string> String_Vector;

    typedef std::map<std::string,std::string> String_Map;

    typedef std::map<std::string,std::pair
		     <void*,std::string> > DataBase_Map;

  protected:

    static String_Vector s_searchpaths, s_nocomplains;
    static String_Map    s_filelocations;
    static DataBase_Map  s_databases;

    std::string m_path, m_file;

    SP(File_Type)   p_file;
    SP(MyStrStream) p_stream;
    
    fom::code m_mode;

    static int ListFiles(void *data,int argc,char **argv,char **name);
    static int GetFile(void *data,int argc,char **argv,char **name);

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
    static void SetNoComplains(const String_Vector &names);

    // static member functions
    static bool OpenDB(std::string file);
    static bool CloseDB(std::string file);

    static bool FileInDB(const std::string &name);

  };// end of class My_File

  template<class FileType>
  std::ostream &operator<<(std::ostream &ostr,
			   const My_File<FileType> &file);

}// end of namespace ATOOLS

#define My_In_File ATOOLS::My_File<std::ifstream>
#define My_Out_File ATOOLS::My_File<std::ofstream>

#endif
