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

    enum OpenMode { Temporary =  1,
		    Permanent =  2,
		    Unknown   = 99 };

    enum VectorType { VVertical   =  1,
		      VHorizontal =  2,
		      VUnknown    = 99 };

    enum MatrixType { MNormal     =  1,
		      MTransposed =  2,
		      MUnknown    = 99 };

  private:

    std::vector<std::string> m_comment, m_ignore, m_seperator;

    std::vector<int> m_blank;

    std::string m_filename, m_string;

    std::fstream *p_file;

    OpenMode m_openmode;

    VectorType m_vectortype;
    MatrixType m_matrixtype;

    void Init();

  protected:

    // access functions
    std::vector<std::string> M_Comment();
    std::vector<std::string> M_Ignore();
    std::vector<std::string> M_Seperator();

    std::vector<int> M_Blank();

    std::string M_FileName();
    std::string M_String();

    OpenMode M_OpenMode();

    std::fstream *M_File();

    VectorType M_VectorType();
    MatrixType M_MatrixType();

    // member functions
    bool OpenFile(std::string filename=nullstring,
		  std::ios_base::openmode omode=std::ios_base::in,
		  OpenMode tempomode=Unknown);
    void CloseFile(bool force=false);

    template <class ReadWrite_Type >
    ReadWrite_Type NotDefined();
    template <class ReadWrite_Type >
    ReadWrite_Type Default();

  public:

    // constructors
    Read_Write_Base();
    Read_Write_Base(const char*,const char*,const char*);
    Read_Write_Base(std::string,std::string,std::string);

    // destructor
    ~Read_Write_Base();

    // setters
    void SetCut(std::string);
    void SetSeperator(std::string);
    void SetSeperator(std::vector<std::string>);
    void AddSeperator(std::string);
    void AddSeperator(std::vector<std::string>);
    void SetComment(std::string);
    void SetComment(std::vector<std::string>);
    void AddComment(std::string);
    void AddComment(std::vector<std::string>);
    void SetIgnore(std::string);
    void SetIgnore(std::vector<std::string>);
    void AddIgnore(std::string);
    void AddIgnore(std::vector<std::string>);

    void SetBlank(int);
    void SetBlank(std::vector<int>);
    void AddBlank(int);
    void AddBlank(std::vector<int>);

    void SetFileName(std::string);
    void SetString(std::string);

    void SetOpenMode(OpenMode _m_openmode);

    void SetVectorType(VectorType _m_vectortype);
    void SetMatrixType(MatrixType _m_matrixtype);

    // getters
    std::vector<std::string> GetComment();
    std::vector<std::string> GetIgnore();
    std::vector<std::string> GetSeperator();

    std::vector<int>         GetBlank();

    std::string              GetFileName();
    std::string              GetString();

    OpenMode                 GetOpenMode();

    VectorType               GetVectorType();
    MatrixType               GetMatrixType();

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
    \brief Reads standard data types from strings or files

    The Read_Write_Base class currently supports the types 
    int, long int, float, double, and std::string. According to the 
    specific purpose the data may either be read from an input string
    or from a file. It is also possible to read multiple inputs, 
    which have to be seperated by blanks, 
    or to read arrays where the columns must be seperated 
    by a chooseable seperator. According to that the output is passed 
    either as a std::vector<Read_Type> or as a std::vector<std::vector<Read_Type>>.
    If not indicated otherwise each function using std::string as input parameter type
    will also exist for const char*.
  */

  /*!
    \fn Read_Write_Base::Read_Write_Base()
    \brief Default constructor

    Initializes the Read_Write_Base with the default cut "=", the default seperator ";" 
    and the default comment indicator "#".
  */

  /*!
    \fn Read_Write_Base::Read_Write_Base(std::string,std::string,std::string)
    \brief Special constructor

    Initializes the Read_Write_Base with the given cut, the given seperator 
    and the given comment indicator. 
  */

  /*!
    \fn Read_Write_Base::SetCut(std::string)
    \brief Sets the cut
  */

  /*!
    \fn Read_Write_Base::SetSeperator(std::string)
    \brief Sets the seperator
  */

  /*!
    \fn Read_Write_Base::SetComment(std::string)
    \brief Sets the comment indicator
  */

  /*!
    \fn Read_Write_Base::SetString(std::string)
    \brief Sets the current filename
  */

  /*!
    \fn Read_Write_Base::SetString(std::string)
    \brief Sets the current input string
  */

  /*!
    \fn Read_Write_Base::GetCut()
    \brief Returns the current cut
  */

  /*!
    \fn Read_Write_Base::GetSeperator()
    \brief Returns the current seperator
  */

  /*!
    \fn Read_Write_Base::GetComment()
    \brief Returns the current comment indicator
  */

  /*!
    \fn Read_Write_Base::GetString()
    \brief Returns the current filename
  */

  /*!
    \fn Read_Write_Base::GetString()
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
    On the basis of a comparison to Read_Write_Base::Default() one is thus able to decide
    whether the reading has been succesful.
  */

} // end of namespace ATOOLS

#include "Read_Write_Base.IC"

#endif
