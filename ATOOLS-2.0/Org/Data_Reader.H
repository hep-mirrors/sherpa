#ifndef Data_Reader_H
#define Data_Reader_H

#include "Read_Write_Base.H"
#include <map>

namespace ATOOLS {

  class Data_Reader: public Read_Write_Base {
  private:

    std::map<std::string,std::string> m_tags;

    void KillComments(std::string& buffer);

    std::string KillBlanks(std::string& buffer);
    std::string HighlightSeperator(std::string& buffer);

    std::string ReplaceTags(std::string tag);

  protected:

    template <class Read_Type >
    Read_Type M_ReadFromFile(std::string parameter,std::string filename);
    template <class Read_Type >
    Read_Type M_ReadFromString(std::string parameter,std::string &inputstring);

    template <class Read_Type >
    std::vector<Read_Type> 
    M_VectorFromFile(std::string parameter,std::string filename,VectorTypeID tempvtype);
    template <class Read_Type >
    std::vector<Read_Type>
    M_VectorFromString(std::string parameter,std::string inputstring,VectorTypeID tempvtype);

    template <class Read_Type >
    std::vector< std::vector<Read_Type> > 
    M_MatrixFromFile(std::string parameter,std::string filename,MatrixTypeID tempmtype);
    template <class Read_Type >
    std::vector< std::vector<Read_Type> >
    M_MatrixFromString(std::string parameter,std::string inputstring,MatrixTypeID tempmtype);

  public:

    // constructors
    Data_Reader();
    Data_Reader(const std::string _m_cut,
		const std::string _m_seperator,
		const std::string _m_comment);

    // member functions
    template <class Read_Type >
    bool ReadFromFile(Read_Type& value,std::string parameter=nullstring,std::string filename=noinputtag);
    template <class Read_Type >
    bool ReadFromString(Read_Type& value,std::string parameter=nullstring,std::string inputstring=noinputtag);

    template <class Read_Type >
    bool VectorFromFile(std::vector<Read_Type>& values,std::string parameter=nullstring,
			std::string filename=noinputtag,VectorTypeID tempvtype=VUnknown);
    template <class Read_Type >
    bool VectorFromString(std::vector<Read_Type>& values,std::string parameter=nullstring,
			  std::string inputstring=noinputtag,VectorTypeID tempvtype=VUnknown);

    template <class Read_Type >
    bool MatrixFromFile(std::vector<std::vector<Read_Type> >& values,std::string parameter=nullstring,
			std::string filename=noinputtag,MatrixTypeID tempmtype=MUnknown);
    template <class Read_Type >
    bool MatrixFromString(std::vector<std::vector<Read_Type> >& values,std::string parameter=nullstring,
			  std::string inputstring=noinputtag,MatrixTypeID tempmtype=MUnknown);

    // inline functions
    inline const std::string &Tag(const std::string &tag) const;

    inline const std::map<std::string,std::string> &Tags() const;

    inline void SetTag(const std::string &tag,const std::string &value);
    inline void SetTags(const std::map<std::string,std::string> &tags);

    inline void AddTag(const std::string &tag,const std::string &value);
    inline void AddTags(const std::map<std::string,std::string> &tags);

  }; // end of class Data_Reader

  //////////////////
  // Doxygen part //
  //////////////////

  /*!
    \file Data_Reader.H
    \brief Declares the class ATOOLS::Data_Reader
  */

  /*!
    \class Data_Reader
    \brief Reads standard data types from strings or files

    The Data_Reader class currently supports the types 
    int, long int, float, double, and std::string. According to the 
    specific purpose the data may either be read from an input string
    or from a file. 
    It is also possible to read multiple inputs, such as vectors and matrices.
    According to the specific purpose in a vector the values may either be
    arranged as a line or as a column. 
    A matrix may either be read normally or as a transposed matrix.
  */

  /*!
    \fn Data_Reader::Data_Reader()
    \brief Default constructor

    Initializes the Data_Reader with the default cut "=", 
    the default column seperator ";" and the default comment indicator "#".
  */

  /*!
    \fn Data_Reader::Data_Reader(std::string,std::string,std::string)
    \brief Special constructor

    Initializes the Data_Reader with the given cut, 
    the given array seperator and the given comment indicator. 
  */

  /*!
    \fn Data_Reader::ReadFromFile(std::string = nullstring, std::string = noinputtag);
    \brief Reads a single Read_Type from a file

    ReadFromFile reads a single Read_Type from a file and stores it into the first 
    argument. The second argument specifies the pattern that tags the input value.
    The input file is specified as an absolute or relative path in the third argument.
    A newline terminates the search for the value even if the pattern has been found.
  */

  /*!
    \fn Data_Reader::ReadFromString(std::string = nullstring, std::string = noinputtag);
    \brief Reads a single Read_Type from a string

    ReadFromString reads a single Read_Type from a string.
    The second argument specifies the pattern that is used to tag the input value.
    The input string is specified in the third argument.
  */

  /*!
    \fn Data_Reader::VectorFromFile(std::string = nullstring, std::string = noinputtag, 
                                    std::vector<std::string> = std::vector<std::string>(0));
    \brief Reads a std::vector<Read_Type> from a file

    VectorFromFile reads a std::vector<Read_Type> from a file.
  */

  /*!
    \fn Data_Reader::VectorFromString(std::string = nullstring, std::string = noinputtag, 
                                      std::vector<std::string> = std::vector<std::string>(0));
    \brief Reads a std::vector<Read_Type> from a string

    VectorFromString reads a std::vector<Read_Type> from a string.
  */

  /*!
    \fn Data_Reader::MatrixFromFile(std::string = nullstring, std::string = noinputtag, 
				  std::vector<std::string> = std::vector<std::string>(0));
    \brief Reads a std::vector<std::vector<Read_Type>> from a file

    MatrixFromFile reads a std::vector<std::vector<Read_Type>> from a file.
  */

  /*!
    \fn Data_Reader::MatrixFromString(std::string = nullstring, std::string = noinputtag, 
				    std::vector<std::string> = std::vector<std::string>(0));
    \brief Reads a std::vector<std::vector<Read_Type>> from a string

    MatrixFromString reads a std::vector<std::vector<Read_Type>> from a string
  */

} // end of namespace ATOOLS

#include "Data_Reader.inl.H"

#endif
