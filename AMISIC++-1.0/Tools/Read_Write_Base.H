#ifndef Read_Write_Base_H
#define Read_Write_Base_H

#include "My_Vector.H"
#include "My_String.H"
#include "My_Filestream.H"

namespace ATOOLS {

  const std::string nullstring = std::string("");
  const std::string defaultcut = std::string("=");
  const std::string defaultsep = std::string(";");
  const std::string defaultcom = std::string("#");
  const std::string noinputtag = std::string("°");

  const int defaultblank = 32;
  const int defaulttab   = 9;

  class Read_Write_Base {
  public:

    enum OpenModeID { Temporary =  1,
		      Permanent =  2,
		      Unknown   = 99 };

    enum VectorTypeID { VVertical   =  1,
			VHorizontal =  2,
			VUnknown    = 99 };

    enum MatrixTypeID { MNormal     =  1,
			MTransposed =  2,
			MUnknown    = 99 };

  private:

    std::string              m_filename, m_string;
    std::vector<std::string> m_comment, m_ignore, m_seperator, m_filecontent;
    std::vector<int>         m_blank;

    std::fstream *p_file;

    OpenModeID m_openmode;

    VectorTypeID m_vectortype;
    MatrixTypeID m_matrixtype;

    void Init();

  protected:

    const std::vector<std::string> &FileContent() const;
    std::fstream                   *File();

    bool OpenFile(std::string filename=nullstring,
		  std::ios_base::openmode omode=std::ios_base::in,
		  OpenModeID tempomode=Unknown);
    void CloseFile(bool force=false);

    const OpenModeID OpenMode() const;

    template <class ReadWrite_Type >
    ReadWrite_Type NotDefined();
    template <class ReadWrite_Type >
    ReadWrite_Type Default();

  public:

    // constructors
    Read_Write_Base();
    Read_Write_Base(std::string _m_cut,std::string _m_seperator,std::string _m_comment);
    Read_Write_Base(const char *_m_cut,const char *_m_seperator,const char *_m_comment);

    // destructor
    ~Read_Write_Base();

    // setters
    void SetSeperator(const std::string _m_seperator);
    void SetSeperator(const std::vector<std::string> &_m_seperator);
    void AddSeperator(const std::string _m_seperator);
    void AddSeperator(const std::vector<std::string> &_m_seperator);

    void SetComment(const std::string _m_comment);
    void SetComment(const std::vector<std::string> &_m_comment);
    void AddComment(const std::string _m_comment);
    void AddComment(const std::vector<std::string> &_m_comment);

    void SetIgnore(const std::string _m_ignore);
    void SetIgnore(const std::vector<std::string> &_m_ignore);
    void AddIgnore(const std::string _m_ignore);
    void AddIgnore(const std::vector<std::string> &_m_ignore);

    void SetBlank(const int _m_blank);
    void SetBlank(const std::vector<int> &_m_blank);
    void AddBlank(const int _m_blank);
    void AddBlank(const std::vector<int> &_m_blank);

    void SetFileName(const std::string _m_filename);
    void SetString(const std::string _m_string);

    void SetOpenMode(const OpenModeID _m_openmode);

    void SetVectorType(const VectorTypeID _m_vectortype);
    void SetMatrixType(const MatrixTypeID _m_matrixtype);

    // getters
    const std::vector<std::string> &Comment() const;
    const std::vector<std::string> &Ignore() const;
    const std::vector<std::string> &Seperator() const;

    const std::vector<int> &Blank() const;

    const std::string FileName() const;
    const std::string String() const;

    const VectorTypeID VectorType() const;
    const MatrixTypeID MatrixType() const;

  }; // end of class Read_Write_Base

  //////////////////
  // Doxygen part //
  //////////////////

  /*!
    \file Read_Write_Base.H
    \brief Declares the class ATOOLS::Read_Write_Base
  */

  /*!
    \class Read_Write_Base
    \brief 
  */

  /*!
    \fn Read_Write_Base::Read_Write_Base()
    \brief Default constructor

    Initializes Read_Write_Base with the default cut "=", the default seperator ";" 
    and the default comment indicator "#".
  */

  /*!
    \fn Read_Write_Base::Read_Write_Base(std::string _m_cut,std::string _m_seperator,std::string _m_comment)
    \brief Special constructor

    Initializes Read_Write_Base with the cut _m_cut, the seperator _m_seperator 
    and the comment indicator _m_comment. 
  */

  /*!
    \fn Read_Write_Base::SetSeperator(std::string _m_seperator)
    \brief Sets the seperator
  */

  /*!
    \fn Read_Write_Base::SetComment(std::string _m_comment)
    \brief Sets the comment indicator
  */

  /*!
    \fn Read_Write_Base::SetFileName(std::string _m_filename)
    \brief Sets the current file name
  */

  /*!
    \fn Read_Write_Base::SetString(std::string _m_string)
    \brief Sets the current input string
  */

  /*!
    \fn Read_Write_Base::Seperator()
    \brief Returns the current seperators
  */

  /*!
    \fn Read_Write_Base::Comment()
    \brief Returns the current comment indicators
  */

  /*!
    \fn Read_Write_Base::FileName()
    \brief Returns the current filename
  */

  /*!
    \fn Read_Write_Base::String()
    \brief Returns the current input string
  */

} // end of namespace ATOOLS

#include "Read_Write_Base.IC"

#endif
