#ifndef ATOOLS_Org_Read_Write_Base_H
#define ATOOLS_Org_Read_Write_Base_H

#include "ATOOLS/Org/File_IO_Base.H"
#include "ATOOLS/Math/Algebra_Interpreter.H"

#include <map>

#include <algorithm> // for find function

namespace ATOOLS {

  const std::string defaultwsep(" ");
  const std::string defaultlsep(";");
  const std::string defaultcom("#");

  const char defaultblank(' ');
  const char defaulttab('\t');

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

  typedef std::vector<char> Char_Vector;

  typedef std::vector<std::string>   String_Vector;
  typedef std::vector<String_Vector> String_Matrix;
  typedef std::vector<String_Matrix> String_Cube;

  typedef std::map<std::string,String_Vector> Buffer_Map;
  typedef std::map<std::string,std::string>   String_Map;

  class Read_Write_Base: public File_IO_Base {
  private:

    static Buffer_Map s_buffermap;

    String_Vector m_comment, m_wordsep, m_linesep;
    String_Vector m_ignore, m_filebegin, m_fileend;
    Char_Vector   m_blank;
    String_Cube   m_filecontent;
    
    char m_escape, m_namesplit;

    vtc::code m_vectortype;
    mtc::code m_matrixtype;

    bool m_ignorecase, m_ignoreblanks, m_exactmatch, m_interprete, m_cmode;
    bool m_allownans;

    size_t m_occurrence;

    Algebra_Interpreter *p_interpreter;    

    void Init();

  protected:

    template <class Type> Type Default() const;

    size_t Find(std::string input,std::string parameter,
		size_t &length) const;

    std::string &KillBlanks(std::string &buffer) const;
    std::string &KillComments(std::string &buffer) const;
    std::string &KillIgnored(std::string &buffer) const;

    inline String_Matrix &FileContent(const size_t i=0)
    { return m_filecontent[i]; }

    char PrevChar(String_Vector &buffer,
		  int &line,int &pos) const;
    char NextChar(String_Vector &buffer,
		  int &line,int &pos) const;
    void InterpreteBuffer(String_Vector &buffer,
			  int &line,int &pos,
			  const int level,const bool keep);
    void InterpreteBuffer(String_Vector &buffer);
    void SplitInFileName(const size_t &i);

  public:

    // constructors
    Read_Write_Base(const unsigned int infiles,
		    const unsigned int outfiles);
    Read_Write_Base(const unsigned int infiles,const unsigned int outfiles,
		    const std::string &wordsep,const std::string &linesep,
		    const std::string &comment,const std::string &ignore="");

    // destructor
    ~Read_Write_Base();

    // member functions
    bool OpenInFile(const unsigned int i=0,const int mode=0);
    bool OpenOutFile(const unsigned int i=0);

    void CloseInFile(const unsigned int i=0,const int mode=0);
    void CloseOutFile(const unsigned int i=0,const int mode=0);

    std::string StripEscapes(const std::string &buffer) const;

    void AddFileContent(std::string line,const unsigned int i=0);

    bool IsBlank(const char &ch) const;

    size_t IsComment(const std::string &ch) const;
    size_t IsWordSeparator(const std::string &ch) const;
    size_t IsLineSeparator(const std::string &ch) const;

    size_t Find(std::string input,std::string parameter) const;

    // setters
    inline void SetBlank(const char blank)
    { m_blank.clear(); m_blank.push_back(blank); }
    inline void SetBlank(const Char_Vector &blank)
    { m_blank=blank; }
    inline void AddBlank(const char blank)
    { m_blank.push_back(blank); }
    inline void AddBlank(const Char_Vector &blank)
    { m_blank.insert(m_blank.end(),blank.begin(),blank.end()); }

    inline void SetComment(const std::string comment)
    { m_comment.clear(); m_comment.push_back(comment); }
    inline void SetComment(const String_Vector &comment)
    { m_comment=comment; }
    inline void AddComment(const std::string comment)
    { m_comment.push_back(comment); }
    inline void AddComment(const String_Vector &comment)
    { m_comment.insert(m_comment.end(),comment.begin(),comment.end()); }

    inline void SetWordSeparator(const std::string separator)
    { m_wordsep.clear(); m_wordsep.push_back(separator); }
    inline void SetWordSeparator(const String_Vector &separator)
    { m_wordsep=separator; }
    inline void AddWordSeparator(const std::string separator)
    { m_wordsep.push_back(separator); }
    inline void AddWordSeparator(const String_Vector &separator)
    { m_wordsep.insert(m_wordsep.end(),
		       separator.begin(),separator.end()); }

    inline void SetLineSeparator(const std::string separator)
    { m_linesep.clear(); m_linesep.push_back(separator); }
    inline void SetLineSeparator(const String_Vector &separator)
    { m_linesep=separator; }
    inline void AddLineSeparator(const std::string separator)
    { m_linesep.push_back(separator); }
    inline void AddLineSeparator(const String_Vector &separator)
    { m_linesep.insert(m_linesep.end(),
		       separator.begin(),separator.end()); }

    inline void SetIgnore(const std::string ignore)
    { m_ignore.clear(); m_ignore.push_back(ignore); }
    inline void SetIgnore(const String_Vector &ignore)
    { m_ignore=ignore; }
    inline void AddIgnore(const std::string ignore)
    { m_ignore.push_back(ignore); }
    inline void AddIgnore(const String_Vector &ignore)
    { m_ignore.insert(m_ignore.end(),ignore.begin(),ignore.end()); }

    inline void SetOccurrence(const size_t occurrence)
    { m_occurrence=occurrence; }

    inline void SetVectorType(const vtc::code vectortype)
    { m_vectortype=vectortype; }
    inline void SetMatrixType(const mtc::code matrixtype)
    { m_matrixtype=matrixtype; }

    inline void SetAllowNans(const bool allownans)
    { m_allownans=allownans; }

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
    inline void SetNameSplit(const char namesplit)
    { m_namesplit=namesplit; }

    // getters
    inline const Char_Vector &Blank() const
    { return m_blank; }

    inline const String_Vector &Comment() const
    { return m_comment; }
    inline const String_Vector &WordSeparator() const
    { return m_wordsep; }
    inline const String_Vector &LineSeparator() const
    { return m_linesep; }

    inline const String_Vector &Ignore() const
    { return m_ignore; }
    inline const String_Vector &FileBegin() const
    { return m_filebegin; }
    inline const String_Vector &FileEnd() const
    { return m_fileend; }

    inline size_t Occurrence() const
    { return m_occurrence; }

    inline vtc::code VectorType() const
    { return m_vectortype; }
    inline mtc::code MatrixType() const
    { return m_matrixtype; }

    inline bool AllowNans() const
    { return m_allownans; }

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
    inline char NameSplit() const
    { return m_namesplit; }

    inline bool RereadInFile(const unsigned int i=0)
    { CloseInFile(i,1); return OpenInFile(i,1); }
    inline bool RescanInFile(const unsigned int i=0)
    { CloseInFile(i,0); return OpenInFile(i,2); }
    
    inline Algebra_Interpreter *Interpreter() const
    { return p_interpreter; }

    inline const String_Vector &BufferContent(const size_t i=0) const
    { return s_buffermap[InputPath(i)+InputFile(i)]; }

  };// end of class Read_Write_Base

}// end of namespace ATOOLS

#endif
