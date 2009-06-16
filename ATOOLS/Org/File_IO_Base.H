#ifndef ATOOLS_Org_File_IO_Base_H
#define ATOOLS_Org_File_IO_Base_H

#include "ATOOLS/Org/My_File.H"
#include <vector>

#ifdef __GNUC__
#if __GNUC__ < 3
#define AT(VECTOR,ELEMENT) VECTOR[ELEMENT]
#else
#define AT(VECTOR,ELEMENT) VECTOR.at(ELEMENT)
#endif
#else
#define AT(VECTOR,ELEMENT) VECTOR.at(ELEMENT)
#endif

namespace ATOOLS {

  class File_IO_Base {
  private:

    std::vector<My_In_File>  m_infiles;
    std::vector<My_Out_File> m_outfiles;

  protected:

    inline My_In_File &InFile(const unsigned int i=0)  
    { return AT(m_infiles,i);  }
    inline My_Out_File &OutFile(const unsigned int i=0)  
    { return AT(m_outfiles,i);  }

  public:

    // constructor
    File_IO_Base(const unsigned int inputfiles=1,
		 const unsigned int outputfiles=1);

    // destructor
    virtual ~File_IO_Base();

    // inline functions
    inline void SetInputPath(const std::string inputpath,
			     const unsigned int i=0)
    { AT(m_infiles,i).SetPath(inputpath); }
    inline void SetInputFile(const std::string inputfile,
			     const unsigned int i=0)
    { AT(m_infiles,i).SetFile(inputfile); }
    inline void SetInFileMode(const fom::code infilemode,
			      const unsigned int i=0)
    { AT(m_infiles,i).SetMode(infilemode); }

    inline void SetOutputPath(const std::string outputpath,
			     const unsigned int i=0)
    { AT(m_outfiles,i).SetPath(outputpath); }
    inline void SetOutputFile(const std::string outputfile,
			     const unsigned int i=0)
    { AT(m_outfiles,i).SetFile(outputfile); }
    inline void SetOutFileMode(const fom::code outfilemode,
			       const unsigned int i=0)
    { AT(m_outfiles,i).SetMode(outfilemode); }

    inline const std::string &InputPath(const unsigned int i=0) const
    { return AT(m_infiles,i).Path(); }
    inline const std::string &InputFile(const unsigned int i=0) const
    { return AT(m_infiles,i).File(); }
    inline fom::code InFileMode(const unsigned int i=0) const  
    { return AT(m_infiles,i).Mode(); }

    inline const std::string &OutputPath(const unsigned int i=0) const
    { return AT(m_outfiles,i).Path(); }
    inline const std::string &OutputFile(const unsigned int i=0) const
    { return AT(m_outfiles,i).File(); }
    inline fom::code OutFileMode(const unsigned int i=0) const
    { return AT(m_outfiles,i).Mode(); }

    inline std::ifstream *InputStream(const unsigned int i=0)  
    { return AT(m_infiles,i)();  }
    inline std::ofstream *OutputStream(const unsigned int i=0) 
    { return AT(m_outfiles,i)(); }

  }; // end of class File_IO_Base

} // end of namespace ATOOLS

#endif
