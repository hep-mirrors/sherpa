#ifndef File_IO_Base_H
#define File_IO_Base_H

#include <fstream>
#include <string>
#include <vector>

namespace ATOOLS {

  const std::string nullstring=std::string("");

  class File_IO_Base {
  public:

    enum OpenModeID { Temporary =  1,
		      Permanent =  2,
		      Unknown   = 99 };
    
  protected:

    std::vector<std::string> m_inputpath, m_inputfile;
    std::vector<std::string> m_outputpath, m_outputfile;

    std::vector<std::ifstream*> m_infile;
    std::vector<std::ofstream*> m_outfile;

    std::vector<OpenModeID> m_infilemode;
    std::vector<OpenModeID> m_outfilemode;

    virtual bool OpenInFile(const unsigned int i);
    virtual bool OpenOutFile(const unsigned int i);

    virtual void CloseInFile(const unsigned int i,const bool force=false);
    virtual void CloseOutFile(const unsigned int i,const bool force=false);

  public:

    // constructor
    File_IO_Base(const unsigned int inputfiles=1,
		 const unsigned int outputfiles=1);

    // destructor
    virtual ~File_IO_Base();

    // setters
    inline void SetInputPath(const std::string _m_inputpath,const unsigned int i=0);
    inline void SetInputFile(const std::string _m_inputfile,const unsigned int i=0);
    inline void SetOutputPath(const std::string _m_outputpath,const unsigned int i=0);
    inline void SetOutputFile(const std::string _m_outputfile,const unsigned int i=0);

    inline void SetInFileMode(const OpenModeID _m_infilemode,const unsigned int i=0);
    inline void SetOutFileMode(const OpenModeID _m_outfilemode,const unsigned int i=0);

    // getters
    inline const std::string InputPath(const unsigned int i=0) const;
    inline const std::string InputFile(const unsigned int i=0) const;
    inline const std::string OutputPath(const unsigned int i=0) const;
    inline const std::string OutputFile(const unsigned int i=0) const;

    inline const OpenModeID InFileMode(const unsigned int i=0) const;
    inline const OpenModeID OutFileMode(const unsigned int i=0) const;

    inline std::ifstream *InFile(const unsigned int i=0);
    inline std::ofstream *OutFile(const unsigned int i=0);

    // member functions
    bool CheckInputPath(const unsigned int i=0);
    bool CheckInputFile(const unsigned int i=0);
    bool CheckOutputPath(const unsigned int i=0);
    bool CheckOutputFile(const unsigned int i=0);

  }; // end of class File_IO_Base

} // end of namespace ATOOLS

#include "File_IO_Base.inl.H"

#endif
