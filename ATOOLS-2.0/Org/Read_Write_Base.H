#ifndef Read_Write_Base_H
#define Read_Write_Base_H

#include "File_IO_Base.H"
#include "Algebra_Interpreter.H"

#include <map>

namespace ATOOLS {

  const std::string defaultcut("=");
  const std::string defaultwsep(" ");
  const std::string defaultlsep(";");
  const std::string defaultcom("#");
  const std::string noinputtag("°");

  const int defaultblank = 32;
  const int defaulttab   = 9;

  struct vtc {

    enum code { 
      vertical   =  1,
      horizontal =  2,
      unknown    = 99 
    };

  };// end of struct vtc

  struct mtc {

    enum code { 
      normal     =  1,
      transposed =  2,
      unknown    = 99
    };

  };// end of struct vtc

  class Algebra_Interpreter;

  class Read_Write_Base: public File_IO_Base, public Tag_Replacer {
  private:

    static std::map<std::string,std::vector<std::string>*> s_buffermap;

    static std::vector<std::string> s_commandline;

    std::map<std::string,std::string> m_tags;

    std::vector<std::string> m_comment, m_ignore, m_wordsep, m_linesep;
    std::vector<int>         m_blank;

    std::vector<std::vector<std::string> > m_buffer, m_filecontent;

    std::vector<std::pair<int,size_t> > m_ifresults;

    std::string              m_string;
    std::vector<std::string> m_filebegin, m_fileend;
    
    char m_escape;

    vtc::code m_vectortype;
    mtc::code m_matrixtype;

    bool m_ignorecase, m_ignoreblanks, m_exactmatch, m_interprete, m_cmode;

    size_t m_occurrence;

    Algebra_Interpreter *p_interpreter;    

    void Init();

  protected:

    bool m_allownans, m_addcommandline;

    template <class Type> Type Default();

    size_t Find(std::string input,std::string parameter,
		size_t &length) const;

    std::string ReplaceTags(std::string &expr) const;    
    std::string &Interprete(std::string &lastline);

  public:

    // constructors
    Read_Write_Base(const unsigned int infiles,
		    const unsigned int outfiles);
    Read_Write_Base(const unsigned int infiles,
		    const unsigned int outfiles,
		    const std::string cut,
		    const std::string wordsep,
		    const std::string linesep,
		    const std::string comment);

    // destructor
    ~Read_Write_Base();

    // member functions
    bool OpenInFile(const unsigned int i=0,const bool force=false);
    bool OpenOutFile(const unsigned int i=0);

    void CloseInFile(const unsigned int i=0,const bool force=false);
    void CloseOutFile(const unsigned int i=0,const bool force=false);

    std::string StripEscapes(const std::string &buffer) const;

    inline const std::vector<std::string> 
    &FileContent(const size_t i=0) const
    { return m_filecontent[i]; }

    inline static const std::vector<std::string> &CommandLine()
    { return s_commandline; }

    // setters
    inline void SetWordSeparator(const std::string separator)
    { m_wordsep.clear(); m_wordsep.push_back(separator); }
    inline void SetWordSeparator(const std::vector<std::string> &separator)
    { m_wordsep=separator; }
    inline void AddWordSeparator(const std::string separator)
    { m_wordsep.push_back(separator); }
    inline void AddWordSeparator(const std::vector<std::string> &separator)
    { m_wordsep.insert(m_wordsep.end(),
		       separator.begin(),separator.end()); }

    inline void SetLineSeparator(const std::string separator)
    { m_linesep.clear(); m_linesep.push_back(separator); }
    inline void SetLineSeparator(const std::vector<std::string> &separator)
    { m_linesep=separator; }
    inline void AddLineSeparator(const std::string separator)
    { m_linesep.push_back(separator); }
    inline void AddLineSeparator(const std::vector<std::string> &separator)
    { m_linesep.insert(m_linesep.end(),
		       separator.begin(),separator.end()); }

    inline void SetComment(const std::string comment)
    { m_comment.clear(); m_comment.push_back(comment); }
    inline void SetComment(const std::vector<std::string> &comment)
    { m_comment=comment; }
    inline void AddComment(const std::string comment)
    { m_comment.push_back(comment); }
    inline void AddComment(const std::vector<std::string> &comment)
    { m_comment.insert(m_comment.end(),comment.begin(),comment.end()); }

    inline void SetIgnore(const std::string ignore)
    { m_ignore.clear(); m_ignore.push_back(ignore); }
    inline void SetIgnore(const std::vector<std::string> &ignore)
    { m_ignore=ignore; }
    inline void AddIgnore(const std::string ignore)
    { m_ignore.push_back(ignore); }
    inline void AddIgnore(const std::vector<std::string> &ignore)
    { m_ignore.insert(m_ignore.end(),ignore.begin(),ignore.end()); }

    inline void SetBlank(const int blank)
    { m_blank.clear(); m_blank.push_back(blank); }
    inline void SetBlank(const std::vector<int> &blank)
    { m_blank=blank; }
    inline void AddBlank(const int blank)
    { m_blank.push_back(blank); }
    inline void AddBlank(const std::vector<int> &blank)
    { m_blank.insert(m_blank.end(),blank.begin(),blank.end()); }

    inline void SetFileBegin(const std::string filebegin)
    { m_filebegin.clear(); m_filebegin.push_back(filebegin); }
    inline void SetFileBegin(const std::vector<std::string> &filebegin)
    { m_filebegin=filebegin; }
    inline void AddFileBegin(const std::string filebegin)
    { m_filebegin.push_back(filebegin); }
    inline void AddFileBegin(const std::vector<std::string> &filebegin)
    { m_filebegin.insert(m_filebegin.end(),
			 filebegin.begin(),filebegin.end()); }

    inline void SetFileEnd(const std::string fileend)
    { m_fileend.clear(); m_fileend.push_back(fileend); }
    inline void SetFileEnd(const std::vector<std::string> &fileend)
    { m_fileend=fileend; }
    inline void AddFileEnd(const std::string fileend)
    { m_fileend.push_back(fileend); }
    inline void AddFileEnd(const std::vector<std::string> &fileend)
    { m_fileend.insert(m_fileend.end(),
		       fileend.begin(),fileend.end()); }

    inline void SetString(const std::string string)
    { m_string=string; }

    inline void SetOccurrence(const size_t occurrence)
    { m_occurrence=occurrence; }

    inline void SetVectorType(const vtc::code vectortype)
    { m_vectortype=vectortype; }
    inline void SetMatrixType(const mtc::code matrixtype)
    { m_matrixtype=matrixtype; }

    inline void SetAllowNans(const bool allownans)
    { m_allownans=allownans; }
    inline void SetAddCommandLine(const bool commandline)
    { m_addcommandline=commandline; }

    inline void SetFileContent(const size_t &i,const std::string filecontent)
    { m_filecontent[i].clear(); m_filecontent[i].push_back(filecontent); }
    inline void SetFileContent(const size_t &i,
			       const std::vector<std::string> &filecontent)
    {m_filecontent[i]=filecontent; }
    inline void AddFileContent(const size_t &i,
			       const std::string filecontent)
    { m_filecontent[i].clear(); m_filecontent[i].push_back(filecontent); }
    inline void AddFileContent(const size_t &i,
			       const std::vector<std::string> &filecontent)
    { m_filecontent[i].insert(m_filecontent[i].end(),
			      filecontent.begin(),filecontent.end()); }

    inline static void SetCommandLine(const std::string commandline)
    { s_commandline.clear(); s_commandline.push_back(commandline); }

    inline static void SetCommandLine(const std::vector<std::string> 
				      &commandline)
    { s_commandline=commandline; }
      inline static void AddCommandLine(const std::string commandline)
    { s_commandline.push_back(commandline); }
    inline static void AddCommandLine(const std::vector<std::string> 
				      &commandline)
    { s_commandline.insert(s_commandline.end(),
			   commandline.begin(),commandline.end()); }
  
    inline void SetIgnoreCase(const bool ignorecase)
    { m_ignorecase=ignorecase; }
    inline void SetIgnoreBlanks(const bool ignoreblanks)
    { m_ignoreblanks=ignoreblanks; }

    inline void SetExactMatch(const bool exact)
    { m_exactmatch=exact; }
    inline void SetInterprete(const bool interprete)
    { m_interprete=interprete; }
    inline void SetCMode(const bool cmode)
    { m_cmode=cmode; }

    inline void SetEscape(const char escape)
    { m_escape=escape; }

    inline void SetTags(const std::map<std::string,std::string> &tags)
    { m_tags=tags; }
    inline void AddTags(const std::map<std::string,std::string> &tags)
    { 
      for (std::map<std::string,std::string>::const_iterator 
	     tit=tags.begin();tit!=tags.end();++tit) 
	m_tags[tit->first]=tit->second; 
    }
    inline void SetTag(const std::string &tag,const std::string &value)
    { m_tags.clear(); m_tags[tag]=value; }
    inline void AddTag(const std::string &tag,const std::string &value)
    { m_tags[tag]=value; }

    // getters
    inline const std::vector<std::string> &Comment() const
    { return m_comment; }
    inline const std::vector<std::string> &Ignore() const
    { return m_ignore; }
    inline const std::vector<std::string> &WordSeparator() const
    { return m_wordsep; }
    inline const std::vector<std::string> &LineSeparator() const
    { return m_linesep; }

    inline const std::vector<int> &Blank() const
    { return m_blank; }

    inline const std::vector<std::string> &FileBegin() const
    { return m_filebegin; }
    inline const std::vector<std::string> &FileEnd() const
    { return m_fileend; }

    inline const std::string String() const
    { return m_string; }

    inline size_t Occurrence() const
    { return m_occurrence; }

    inline const vtc::code VectorType() const
    { return m_vectortype; }
    inline const mtc::code MatrixType() const
    { return m_matrixtype; }

    inline const bool AllowNans() const
    { return m_allownans; }
    inline const bool AddCommandLine() const
    { return m_addcommandline; }

    inline bool IgnoreCase() const
    { return m_ignorecase; }
    inline bool IgnoreBlanks() const
    { return m_ignoreblanks; }

    inline bool ExactMatch() const
    { return m_exactmatch; }
    inline bool Interprete() const
    { return m_interprete; }
    inline bool CMode() const
    { return m_cmode; }

    inline char Escape() const
    { return m_escape; }

    inline bool RereadInFile(const unsigned int i=0)
    { CloseInFile(i,true); return OpenInFile(i,true); }
    inline bool RescanInFile(const unsigned int i=0)
    { CloseInFile(i,false); return OpenInFile(i,false); }
    
    inline Algebra_Interpreter *Interpreter() const
    { return p_interpreter; }

    inline const std::string &Tag(const std::string &tag) const
    { 
      std::map<std::string,std::string>::const_iterator 
	tit(m_tags.find(tag));
      return tit!=m_tags.end()?nullstring:tit->second;
    }
    inline const std::map<std::string,std::string> &Tags() const
    { return m_tags; }

  };// end of class Read_Write_Base

}// end of namespace ATOOLS

#endif
