#ifndef Read_Write_Base_H
#define Read_Write_Base_H

#include "File_IO_Base.H"

namespace ATOOLS {

  const std::string defaultcut = std::string("=");
  const std::string defaultsep = std::string(";");
  const std::string defaultcom = std::string("#");
  const std::string noinputtag = std::string("°");

  const int defaultblank = 32;
  const int defaulttab   = 9;

  class Read_Write_Base: public File_IO_Base {
  public:

    enum VectorTypeID { VVertical   =  1,
			VHorizontal =  2,
			VUnknown    = 99 };

    enum MatrixTypeID { MNormal     =  1,
			MTransposed =  2,
			MUnknown    = 99 };

  private:

    static std::vector<std::string> s_commandline;

    std::vector<std::string> m_comment, m_ignore, m_separator, m_filecontent;
    std::vector<int>         m_blank;

    std::string              m_string;
    std::vector<std::string> m_filebegin, m_fileend;

    VectorTypeID m_vectortype;
    MatrixTypeID m_matrixtype;

    void Init();

  protected:

    const std::vector<std::string> &FileContent() const;

    static const std::vector<std::string> &CommandLine();

    bool m_allownans, m_addcommandline;

    bool OpenInFile(const unsigned int i=0);
    bool OpenOutFile(const unsigned int i=0);

    void CloseInFile(const unsigned int i=0,const bool force=false);
    void CloseOutFile(const unsigned int i=0,const bool force=false);

    template <class ReadWrite_Type >
    ReadWrite_Type NotDefined();
    template <class ReadWrite_Type >
    ReadWrite_Type Default();

  public:

    // constructors
    Read_Write_Base(const unsigned int infiles,
		    const unsigned int outfiles);
    Read_Write_Base(const unsigned int infiles,
		    const unsigned int outfiles,
		    const std::string _m_cut,
		    const std::string _m_separator,
		    const std::string _m_comment);

    // destructor
    ~Read_Write_Base();

    // setters
    inline void SetSeparator(const std::string _m_separator);
    inline void SetSeparator(const std::vector<std::string> &_m_separator);
    inline void AddSeparator(const std::string _m_separator);
    inline void AddSeparator(const std::vector<std::string> &_m_separator);

    inline void SetComment(const std::string _m_comment);
    inline void SetComment(const std::vector<std::string> &_m_comment);
    inline void AddComment(const std::string _m_comment);
    inline void AddComment(const std::vector<std::string> &_m_comment);

    inline void SetIgnore(const std::string _m_ignore);
    inline void SetIgnore(const std::vector<std::string> &_m_ignore);
    inline void AddIgnore(const std::string _m_ignore);
    inline void AddIgnore(const std::vector<std::string> &_m_ignore);

    inline void SetBlank(const int _m_blank);
    inline void SetBlank(const std::vector<int> &_m_blank);
    inline void AddBlank(const int _m_blank);
    inline void AddBlank(const std::vector<int> &_m_blank);

    inline void SetFileBegin(const std::string _m_filebegin);
    inline void SetFileBegin(const std::vector<std::string> &_m_filebegin);
    inline void AddFileBegin(const std::string _m_filebegin);
    inline void AddFileBegin(const std::vector<std::string> &_m_filebegin);

    inline void SetFileEnd(const std::string _m_fileend);
    inline void SetFileEnd(const std::vector<std::string> &_m_fileend);
    inline void AddFileEnd(const std::string _m_fileend);
    inline void AddFileEnd(const std::vector<std::string> &_m_fileend);

    inline void SetString(const std::string _m_string);

    inline void SetOpenMode(const OpenModeID _m_openmode);

    inline void SetVectorType(const VectorTypeID _m_vectortype);
    inline void SetMatrixType(const MatrixTypeID _m_matrixtype);

    inline void SetAllowNans(const bool _m_killnans);
    inline void SetAddCommandLine(const bool _m_commandline);

    inline void SetFileContent(const std::string _m_filecontent);
    inline void SetFileContent(const std::vector<std::string> &_m_filecontent);

    inline void AddFileContent(const std::string _m_filecontent);
    inline void AddFileContent(const std::vector<std::string> &_m_filecontent);

    inline static void SetCommandLine(const std::string _s_commandline);
    inline static void SetCommandLine(const std::vector<std::string> &_s_commandline);
  
    inline static void AddCommandLine(const std::string _s_commandline);
    inline static void AddCommandLine(const std::vector<std::string> &_s_commandline);
  
    // getters
    inline const std::vector<std::string> &Comment() const;
    inline const std::vector<std::string> &Ignore() const;
    inline const std::vector<std::string> &Separator() const;

    inline const std::vector<int> &Blank() const;

    inline const std::vector<std::string> &FileBegin() const;
    inline const std::vector<std::string> &FileEnd() const;

    inline const std::string String() const;

    inline const VectorTypeID VectorType() const;
    inline const MatrixTypeID MatrixType() const;

    inline const bool AllowNans() const;
    inline const bool AddCommandLine() const;

  };// end of class Read_Write_Base

  //////////////////
  // Doxygen part //
  //////////////////

  /*!
    \file Read_Write_Base.H
    \brief Declares the class ATOOLS::Read_Write_Base
  */

  /*!
    \class Read_Write_Base
    \brief Reads standard data types from strings or files

    The Read_Write_Base class currently supports the types 
    int, long int, float, double, and std::string. According to the 
    specific purpose the data may either be read from an input string
    or from a file. It is also possible to read column vectors or line 
    vectors, as well as matrices. The output is passed either as a 
    std::vector<Read_Type> or as a std::vector<std::vector<Read_Type>>.
  */

  /*!
    \fn Read_Write_Base::Read_Write_Base()
    \brief Default constructor

    Initializes the Read_Write_Base with the default cut "=", 
    the default column separator ";" and the default comment indicator "#".
  */

  /*!
    \fn Read_Write_Base::Read_Write_Base(std::string,std::string,std::string)
    \brief Special constructor

    Initializes the Read_Write_Base with the given cut, the given separator 
    and the given comment indicator. 
  */

  /*!
    \fn Read_Write_Base::SetSeparator(std::string)
    \brief Sets the separator
  */

  /*!
    \fn Read_Write_Base::SetComment(std::string)
    \brief Sets the comment indicator
  */

  /*!
    \fn Read_Write_Base::SetString(std::string)
    \brief Sets the current input string
  */

  /*!
    \fn Read_Write_Base::Cut()
    \brief Returns the current cut
  */

  /*!
    \fn Read_Write_Base::Separator()
    \brief Returns the current separator
  */

  /*!
    \fn Read_Write_Base::Comment()
    \brief Returns the current comment indicator
  */

  /*!
    \fn Read_Write_Base::String()
    \brief Returns the current filename
  */

  /*!
    \fn Read_Write_Base::String()
    \brief Returns the current input string
  */

  /*!
    \fn Read_Write_Base::NotDefined()
    \brief Identifies unreadable input 

    This template is used to decide whether the input parameter
    is defined or not. In principle it should return the value 
    std::stringstream returns when a null argument has been passed 
  */

  /*!
    \fn Read_Write_Base::Default()
    \brief Returns default values

    This template passes a default value to the user routine
    in the case that reading the input string resulted in an error. 
  */

}// end of namespace ATOOLS

#include "Read_Write_Base.inl.H"

#endif
