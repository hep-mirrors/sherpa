#ifndef Read_Write_Base_C
#define Read_Write_Base_C

#include "Read_Write_Base.H"

namespace ATOOLS {

  Read_Write_Base::Read_Write_Base(const unsigned int infiles,
				   const unsigned int outfiles):
    File_IO_Base(infiles,outfiles),
    m_comment(std::vector<std::string>(1,defaultcom)), 
    m_ignore(std::vector<std::string>(1,defaultcut)),
    m_seperator(std::vector<std::string>(1,defaultsep))
  {
    Init();
  }

  Read_Write_Base::Read_Write_Base(const unsigned int infiles,
				   const unsigned int outfiles,
				   const std::string _m_cut,
				   const std::string _m_seperator,
				   const std::string _m_comment):
    File_IO_Base(infiles,outfiles),
    m_comment(std::vector<std::string>(1,_m_comment)), 
    m_ignore(std::vector<std::string>(1,_m_cut)),
    m_seperator(std::vector<std::string>(1,_m_seperator))
  {
    Init();
  }
  
  Read_Write_Base::~Read_Write_Base() {}

  void Read_Write_Base::Init()
  {
    m_blank.push_back(defaultblank);
    m_blank.push_back(defaulttab);
    m_vectortype = VVertical;
    m_matrixtype = MNormal; 
  }
    
  bool Read_Write_Base::OpenInFile(const unsigned int i)
  {  
    if (InputFile(i)==nullstring) return false;
    if (InFileMode(i)==Unknown) SetInFileMode(Temporary);
    if (m_infile[i]==NULL) {
      m_infile[i]=new std::ifstream();	
      m_infile[i]->open((m_inputpath[i]+m_inputfile[i]).c_str()); 
      m_filecontent.clear();
      while (*m_infile[i]) {
	m_filecontent.push_back(std::string(""));
	getline(*m_infile[i],m_filecontent.back());
      }
    }
    return !m_infile[i]->bad();
  }

  bool Read_Write_Base::OpenOutFile(const unsigned int i)
  {  
    if (OutputFile(i)==nullstring) return false;
    if (OutFileMode(i)==Unknown) SetOutFileMode(Permanent);
    if (m_outfile[i]==NULL) {
      m_outfile[i]=new std::ofstream();	
      m_outfile[i]->open((m_outputpath[i]+m_outputfile[i]).c_str());
    }
    return !m_outfile[i]->bad();
  }

  void Read_Write_Base::CloseInFile(const unsigned int i,const bool force)
  { 
    if (m_infile[i]==NULL) return;
    if ((m_infilemode[i]==Permanent)&&(!force)) return;
    m_filecontent.clear();
    m_infile[i]->close(); 
    delete m_infile[i]; 
    m_infile[i]=NULL;
  }

  void Read_Write_Base::CloseOutFile(const unsigned int i,const bool force)
  { 
    if (m_outfile[i]==NULL) return;
    if ((m_outfilemode[i]==Permanent)&&(!force)) return;
    m_outfile[i]->close(); 
    delete m_outfile[i]; 
    m_outfile[i]=NULL;
  }

} // end of namespace ATOOLS

#endif

