#ifndef Read_Write_Base_C
#define Read_Write_Base_C

#include "Read_Write_Base.H"

namespace ATOOLS {

  Read_Write_Base::Read_Write_Base() 
  {
    m_comment.push_back(defaultcom); 
    m_ignore.push_back(defaultcut);
    m_seperator.push_back(defaultsep);
    Init();
  }

  Read_Write_Base::Read_Write_Base(std::string _m_cut, std::string _m_seperator, std::string _m_comment)
  {
    m_comment.push_back(_m_comment); 
    m_ignore.push_back(_m_cut);
    m_seperator.push_back(_m_seperator);
    Init();
  }
  
  Read_Write_Base::Read_Write_Base(const char *_m_cut, const char *_m_seperator, const char *_m_comment)
  {
    m_comment.push_back(std::string(_m_comment)); 
    m_ignore.push_back(std::string(_m_cut));
    m_seperator.push_back(std::string(_m_seperator));
    Init();
  }

  Read_Write_Base::~Read_Write_Base()
  {
    if (p_file!=NULL) CloseFile();
    delete p_file;
  }

  void Read_Write_Base::Init()
  {
    m_blank.push_back(defaultblank);
    m_blank.push_back(defaulttab);
    m_vectortype = VVertical;
    m_matrixtype = MNormal; 
    m_openmode = Temporary; 
    p_file = NULL;
  }
    
  bool Read_Write_Base::OpenFile(std::string filename,std::ios_base::openmode omode,
				 OpenModeID tempomode)
  {  
    if (filename!=nullstring) SetFileName(filename);
    if (tempomode!=Unknown) SetOpenMode(tempomode);
    if (m_filename==nullstring) return false;
    if (m_openmode==Unknown) m_openmode=Temporary;
    if (p_file==NULL) {
      p_file=new std::fstream();	
      p_file->open(m_filename.c_str(),omode); 
      m_filecontent.clear();
      if (omode==std::ios_base::in) {
	while (*p_file) {
	  m_filecontent.push_back(std::string(""));
	  getline(*p_file,m_filecontent.back());
	}
      }
    }
    return !p_file->bad();
  }

  void Read_Write_Base::CloseFile(bool force)
  { 
    if ((m_openmode==Permanent)&&(!force)) return;
    m_filecontent.clear();
    p_file->close(); 
    delete p_file; 
    p_file=NULL;
  }

} // end of namespace ATOOLS

#endif

