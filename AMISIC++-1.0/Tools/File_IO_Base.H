#ifndef File_IO_Base_H
#define File_IO_Base_H

#include "My_String.H"

namespace ATOOLS {

  class File_IO_Base {
  protected:

    std::string m_inputpath, m_inputfile, m_outputpath, m_outputfile;

  public:

    // constructor
    File_IO_Base();

    // destructor
    ~File_IO_Base();

    // setters
    void SetInputPath(const std::string _m_inputpath);
    void SetInputFile(const std::string _m_inputfile);
    void SetOutputPath(const std::string _m_outputpath);
    void SetOutputFile(const std::string _m_outputfile);

    // getters
    const std::string InputPath() const;
    const std::string InputFile() const;
    const std::string OutputPath() const;
    const std::string OutputFile() const;

    // member functions
    bool CheckInputPath();
    bool CheckInputFile();
    bool CheckOutputPath();
    bool CheckOutputFile();

  }; // end of class File_IO_Base

} // end of namespace ATOOLS

#include "File_IO_Base.IC"

#endif
