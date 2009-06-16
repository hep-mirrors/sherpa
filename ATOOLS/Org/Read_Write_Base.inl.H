#ifndef ATOOLS_Org_Read_Write_Base_inl_H
#define ATOOLS_Org_Read_Write_Base_inl_H

#include "ATOOLS/Org/Read_Write_Base.H"
#include "ATOOLS/Org/My_Limits.H"

namespace ATOOLS {

  inline void Read_Write_Base::AddCommandLine(const std::string _s_commandline)
  { 
    s_commandline.push_back(_s_commandline); 
  }

  inline void Read_Write_Base::AddCommandLine(const std::vector<std::string> &_s_commandline)
  { 
    for (unsigned int i=0;i<_s_commandline.size();s_commandline.push_back(_s_commandline[i++])); 
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

  inline void Read_Write_Base::SetTags(const std::map<std::string,std::string> &tags) 
  { 
    m_tags=tags; 
  }
  
  inline void Read_Write_Base::AddTags(const std::map<std::string,std::string> &tags) 
  { 
    for (std::map<std::string,std::string>::const_iterator tit=tags.begin();
	 tit!=tags.end();++tit) m_tags[tit->first]=tit->second; 
  }
  
  inline void Read_Write_Base::SetTag(const std::string &tag,const std::string &value) 
  { 
    m_tags.clear(); 
    m_tags[tag]=value; 
  }
  
  inline void Read_Write_Base::AddTag(const std::string &tag,const std::string &value) 
  { 
    m_tags[tag]=value; 
  }

  inline const std::string &Read_Write_Base::Tag(const std::string &tag) const 
  { 
    return m_tags.find(tag)->second;
  }
  
  inline const std::map<std::string,std::string> &Read_Write_Base::Tags() const
  { 
    return m_tags;
  }

} // end of namespace ATOOLS

#endif
