#ifndef Read_Write_Base_inl_H
#define Read_Write_Base_inl_H

#include "Read_Write_Base.H"
#include "My_Limits.H"

namespace ATOOLS {

  inline void Read_Write_Base::SetSeparator(const std::string _m_separator) 
  { 
    m_separator.clear(); 
    m_separator.push_back(_m_separator); 
  }
  
  inline void Read_Write_Base::SetSeparator(const std::vector<std::string> &_m_separator) 
  { 
    m_separator.clear(); 
    m_separator=_m_separator; 
  }
  
  inline void Read_Write_Base::AddSeparator(const std::string _m_separator) 
  { 
    m_separator.push_back(_m_separator); 
  }
  
  inline void Read_Write_Base::AddSeparator(const std::vector<std::string> &_m_separator) 
  { 
    for (unsigned int i=0;i<_m_separator.size();m_separator.push_back(_m_separator[i++])); 
  }
  
  inline void Read_Write_Base::SetComment(const std::string _m_comment)       
  { 
    m_comment.clear(); 
    m_comment.push_back(_m_comment); 
  }

  inline void Read_Write_Base::SetComment(const std::vector<std::string> &_m_comment)       
  { 
    m_comment.clear(); 
    m_comment=_m_comment; 
  }
  
  inline void Read_Write_Base::AddComment(const std::string _m_comment)       
  { 
    m_comment.push_back(_m_comment); 
  }
  
  inline void Read_Write_Base::AddComment(const std::vector<std::string> &_m_comment)       
  { 
    for (unsigned int i=0;i<_m_comment.size();m_comment.push_back(_m_comment[i++])); 
  }
  
  inline void Read_Write_Base::SetFileContent(const std::string _m_filecontent)
  { 
    m_filecontent.clear(); 
    m_filecontent.push_back(_m_filecontent); 
  }

  inline void Read_Write_Base::SetFileContent(const std::vector<std::string> &_m_filecontent)
  { 
    m_filecontent.clear(); 
    for (unsigned int i=0;i<_m_filecontent.size();m_filecontent.push_back(_m_filecontent[i++]));
  }

  inline void Read_Write_Base::AddFileContent(const std::string _m_filecontent)
  { 
    m_filecontent.push_back(_m_filecontent); 
  }

  inline void Read_Write_Base::AddFileContent(const std::vector<std::string> &_m_filecontent)
  { 
    for (unsigned int i=0;i<_m_filecontent.size();m_filecontent.push_back(_m_filecontent[i++])); 
  }

  inline void Read_Write_Base::SetCommandLine(const std::string _s_commandline)
  { 
    s_commandline.clear(); 
    s_commandline.push_back(_s_commandline); 
  }

  inline void Read_Write_Base::SetCommandLine(const std::vector<std::string> &_s_commandline)
  { 
    s_commandline.clear(); 
    for (unsigned int i=0;i<_s_commandline.size();s_commandline.push_back(_s_commandline[i++]));
  }

  inline void Read_Write_Base::AddCommandLine(const std::string _s_commandline)
  { 
    s_commandline.push_back(_s_commandline); 
  }

  inline void Read_Write_Base::AddCommandLine(const std::vector<std::string> &_s_commandline)
  { 
    for (unsigned int i=0;i<_s_commandline.size();s_commandline.push_back(_s_commandline[i++])); 
  }

  inline void Read_Write_Base::SetIgnore(const std::string _m_ignore)       
  { 
    m_ignore.clear(); 
    m_ignore.push_back(_m_ignore); 
  }

  inline void Read_Write_Base::SetIgnore(const std::vector<std::string> &_m_ignore)       
  { 
    m_ignore.clear(); 
    m_ignore=_m_ignore; 
  }
  
  inline void Read_Write_Base::AddIgnore(const std::string _m_ignore)       
  { 
    m_ignore.push_back(_m_ignore); 
  }

  inline void Read_Write_Base::AddIgnore(const std::vector<std::string> &_m_ignore)       
  { 
    for (unsigned int i=0;i<_m_ignore.size();m_ignore.push_back(_m_ignore[i++])); 
  }

  inline void Read_Write_Base::SetBlank(const int _m_blank)                 
  { 
    m_blank.clear(); 
    m_blank.push_back(_m_blank); 
  }
  
  inline void Read_Write_Base::SetBlank(const std::vector<int> &_m_blank)    
  { 
    m_blank.clear(); 
    m_blank = _m_blank; 
  }

  inline void Read_Write_Base::AddBlank(const int _m_blank)                 
  { 
    m_blank.push_back(_m_blank); 
  }

  inline void Read_Write_Base::AddBlank(const std::vector<int> &_m_blank)    
  { 
    for (unsigned int i=0; i<_m_blank.size();m_blank.push_back(_m_blank[i++])); 
  }

  inline void Read_Write_Base::SetFileBegin(const std::string _m_filebegin)       
  { 
    m_filebegin.clear();
    m_filebegin.push_back(_m_filebegin); 
  }

  inline void Read_Write_Base::AddFileBegin(const std::string _m_filebegin)       
  { 
    m_filebegin.push_back(_m_filebegin); 
  }

  inline void Read_Write_Base::SetFileBegin(const std::vector<std::string> &_m_filebegin)       
  { 
    m_filebegin = _m_filebegin; 
  }

  inline void Read_Write_Base::AddFileBegin(const std::vector<std::string> &_m_filebegin)       
  { 
    for (unsigned int i=0;i<_m_filebegin.size();m_filebegin.push_back(_m_filebegin[i++])); 
  }

  inline void Read_Write_Base::SetFileEnd(const std::string _m_fileend)       
  { 
    m_fileend.clear();
    m_fileend.push_back(_m_fileend); 
  }

  inline void Read_Write_Base::AddFileEnd(const std::string _m_fileend)       
  { 
    m_fileend.push_back(_m_fileend); 
  }

  inline void Read_Write_Base::SetFileEnd(const std::vector<std::string> &_m_fileend)       
  { 
    m_fileend = _m_fileend; 
  }

  inline void Read_Write_Base::AddFileEnd(const std::vector<std::string> &_m_fileend)       
  { 
    for (unsigned int i=0;i<_m_fileend.size();m_fileend.push_back(_m_fileend[i++])); 
  }

  inline void Read_Write_Base::SetString(const std::string _m_string)       
  { 
    m_string = _m_string; 
  }

  inline void Read_Write_Base::SetOccurrence(const size_t _m_occurrence)
  { 
    m_occurrence = _m_occurrence; 
  }

  inline void Read_Write_Base::SetVectorType(const VectorTypeID _m_vectortype)
  { 
    m_vectortype=_m_vectortype; 
  }

  inline void Read_Write_Base::SetMatrixType(const MatrixTypeID _m_matrixtype)
  { 
    m_matrixtype=_m_matrixtype; 
  }

  inline const std::vector<std::string> &Read_Write_Base::Comment() const           
  { 
    return m_comment; 
  }
  
  inline const std::vector<std::string> &Read_Write_Base::Ignore() const           
  { 
    return m_ignore;     
  }
  
  inline const std::vector<std::string> &Read_Write_Base::Separator() const         
  { 
    return m_separator; 
  }
  
  inline const std::vector<int> &Read_Write_Base::Blank() const            
  { 
    return m_blank; 
  }
  
  inline const std::vector<std::string> &Read_Write_Base::FileBegin() const           
  { 
    return m_filebegin; 
  }
  
  inline const std::vector<std::string> &Read_Write_Base::FileEnd() const           
  { 
    return m_fileend; 
  }
  
  inline const std::string Read_Write_Base::String() const           
  { 
    return m_string; 
  }
  
  inline size_t Read_Write_Base::Occurrence() const           
  { 
    return m_occurrence; 
  }
  
  inline const Read_Write_Base::VectorTypeID Read_Write_Base::VectorType() const       
  { 
    return m_vectortype; 
  }
  
  inline const Read_Write_Base::MatrixTypeID Read_Write_Base::MatrixType() const       
  { 
    return m_matrixtype; 
  }
  
  inline void Read_Write_Base::SetAllowNans(const bool _m_allownans)
  { 
    m_allownans=_m_allownans; 
  }

  inline const bool Read_Write_Base::AllowNans() const
  { 
    return m_allownans; 
  }

  inline void Read_Write_Base::SetAddCommandLine(const bool _m_addcommandline)
  { 
    m_addcommandline=_m_addcommandline; 
  }

  inline const bool Read_Write_Base::AddCommandLine() const
  { 
    return m_addcommandline; 
  }

  template <class Read_Type>
  inline Read_Type Read_Write_Base::NotDefined() 
  { 
    return (Read_Type)0; 
  }
  
  template <class Read_Type>
  inline Read_Type Read_Write_Base::Default() 
  { 
    return (Read_Type)0; 
  }

  inline const std::vector<std::string> &Read_Write_Base::FileContent() const 
  { 
    return m_filecontent; 
  }

  inline const std::vector<std::string> &Read_Write_Base::CommandLine()
  { 
    return s_commandline; 
  }

  template <> inline int Read_Write_Base::NotDefined<int>()          
  { 
    return std::numeric_limits<int>::max();  
  } 
  
  template <> inline unsigned int Read_Write_Base::NotDefined<unsigned int>()
  { 
    return std::numeric_limits<unsigned int>::max();  
  } 

  template <> inline long int Read_Write_Base::NotDefined<long>()        
  { 
    return std::numeric_limits<long int>::max();  
  } 

  template <> inline float Read_Write_Base::NotDefined<float>()       
  { 
    return std::numeric_limits<float>::max();  
  } 

  template <> inline double Read_Write_Base::NotDefined<double>()      
  { 
    return std::numeric_limits<double>::max();  
  } 

  template <> inline std::string Read_Write_Base::NotDefined<std::string>()  
  { 
    return nullstring; 
  } 

  template <> inline int Read_Write_Base::Default<int>()         
  { 
    return std::numeric_limits<int>::max();  
  } 
  
  template <> inline unsigned int Read_Write_Base::Default<unsigned int>()
  { 
    return std::numeric_limits<unsigned int>::max();  
  } 
  
  template <> inline long int Read_Write_Base::Default<long>()     
  { 
    return std::numeric_limits<long int>::max();  
  } 

  template <> inline float Read_Write_Base::Default<float>()    
  { 
    return std::numeric_limits<float>::max();  
  } 

  template <> inline double Read_Write_Base::Default<double>()   
  { 
    return std::numeric_limits<double>::max();  
  } 

  template <> inline std::string Read_Write_Base::Default<std::string>()  
  { 
    return nullstring; 
  } 

  inline bool Read_Write_Base::RereadInFile(const unsigned int i) 
  { 
    CloseInFile(i,true); 
    return OpenInFile(i); 
  }

  inline void Read_Write_Base::SetIgnoreCase(const bool ignorecase) 
  { 
    m_ignorecase=ignorecase; 
  }

  inline bool Read_Write_Base::IgnoreCase() const
  {
    return m_ignorecase;
  }
  
  inline void Read_Write_Base::SetIgnoreBlanks(const bool ignoreblanks) 
  { 
    m_ignoreblanks=ignoreblanks; 
  }
  
  inline bool Read_Write_Base::IgnoreBlanks() const
  {
    return m_ignoreblanks;
  }
  
  inline void Read_Write_Base::SetExactMatch(const bool exact) 
  { 
    m_exactmatch=exact; 
  }

  inline bool Read_Write_Base::ExactMatch() const 
  { 
    return m_exactmatch; 
  }

  inline void Read_Write_Base::SetInterprete(const bool interprete) 
  { 
    m_interprete=interprete; 
  }

  inline bool Read_Write_Base::Interprete() const 
  { 
    return m_interprete; 
  }

  inline void Read_Write_Base::SetCMode(const bool cmode) 
  { 
    m_cmode=cmode; 
  }

  inline bool Read_Write_Base::CMode() const 
  { 
    return m_cmode; 
  }

  inline Algebra_Interpreter *Read_Write_Base::Interpreter() const
  {
    return p_interpreter;
  }

  inline void Read_Write_Base::SetEscape(const char escape)
  { 
    m_escape=escape; 
  }

  inline char Read_Write_Base::Escape() const 
  { 
    return m_escape; 
  }

} // end of namespace ATOOLS

#endif
