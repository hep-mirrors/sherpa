#ifndef Data_Reader_H
#define Data_Reader_H

#include "Read_Write_Base.H"

namespace ATOOLS {

  class Data_Reader: public Read_Write_Base {
  private:

    std::string KillBlanks(std::string&);
    void KillComments(std::string&);

  protected:

    template <class Read_Type >
    Read_Type M_ReadFromFile(std::string parameter,std::string filename);
    template <class Read_Type >
    Read_Type M_ReadFromString(std::string parameter,std::string inputstring);
    template <class Read_Type >
    std::vector<Read_Type> 
    M_VectorFromFile(std::string parameter,std::string filename,VectorType tempvtype);
    template <class Read_Type >
    std::vector<Read_Type>
    M_VectorFromString(std::string parameter,std::string inputstring,VectorType tempvtype);
    template <class Read_Type >
    std::vector< std::vector<Read_Type> > 
    M_ArrayFromFile(std::string parameter,std::string filename,MatrixType tempmtype);
    template <class Read_Type >
    std::vector< std::vector<Read_Type> >
    M_ArrayFromString(std::string parameter,std::string inputstring,MatrixType tempmtype);

  public:

    // constructors
    Data_Reader();
    Data_Reader(const char*,const char*,const char*);
    Data_Reader(std::string,std::string,std::string);

    // member functions
    bool ReadFromFile(int& value,std::string parameter=nullstring,std::string filename=noinputtag);
    bool ReadFromFile(long int& value,std::string parameter=nullstring,std::string filename=noinputtag);
    bool ReadFromFile(float& value,std::string parameter=nullstring,std::string filename=noinputtag);
    bool ReadFromFile(double& value,std::string parameter=nullstring,std::string filename=noinputtag);
    bool ReadFromFile(std::string& value,std::string parameter=nullstring,std::string filename=noinputtag);

    bool ReadFromString(int& value,std::string parameter=nullstring,std::string inputstring=noinputtag);
    bool ReadFromString(long int& value,std::string parameter=nullstring,std::string inputstring=noinputtag);
    bool ReadFromString(float& value,std::string parameter=nullstring,std::string inputstring=noinputtag);
    bool ReadFromString(double& value,std::string parameter=nullstring,std::string inputstring=noinputtag);
    bool ReadFromString(std::string& value,std::string parameter=nullstring,std::string inputstring=noinputtag);

    bool VectorFromFile(std::vector<int>& values,std::string parameter=nullstring,
			std::string filename=noinputtag,VectorType tempvtype=VUnknown);
    bool VectorFromFile(std::vector<long int>& values,std::string parameter=nullstring,
			std::string filename=noinputtag,VectorType tempvtype=VUnknown);
    bool VectorFromFile(std::vector<float>& values,std::string parameter=nullstring,
			std::string filename=noinputtag,VectorType tempvtype=VUnknown);
    bool VectorFromFile(std::vector<double>& values,std::string parameter=nullstring,
			std::string filename=noinputtag,VectorType tempvtype=VUnknown);
    bool VectorFromFile(std::vector<std::string>& values,std::string parameter=nullstring,
			std::string filename=noinputtag,VectorType tempvtype=VUnknown);

    bool VectorFromString(std::vector<int>& values,std::string parameter=nullstring,
			  std::string inputstring=noinputtag,VectorType tempvtype=VUnknown);
    bool VectorFromString(std::vector<long int>& values,std::string parameter=nullstring,
			  std::string inputstring=noinputtag,VectorType tempvtype=VUnknown);
    bool VectorFromString(std::vector<float>& values,std::string parameter=nullstring,
			  std::string inputstring=noinputtag,VectorType tempvtype=VUnknown);
    bool VectorFromString(std::vector<double>& values,std::string parameter=nullstring,
			  std::string inputstring=noinputtag,VectorType tempvtype=VUnknown);
    bool VectorFromString(std::vector<std::string>& values,std::string parameter=nullstring,
			  std::string inputstring=noinputtag,VectorType tempvtype=VUnknown);

    bool ArrayFromFile(std::vector<std::vector<int> >& values,std::string parameter=nullstring,
		       std::string filename=noinputtag,MatrixType tempmtype=MUnknown);
    bool ArrayFromFile(std::vector<std::vector<long int> >& values,std::string parameter=nullstring,
		       std::string filename=noinputtag,MatrixType tempmtype=MUnknown);
    bool ArrayFromFile(std::vector<std::vector<float> >& values,std::string parameter=nullstring,
		       std::string filename=noinputtag,MatrixType tempmtype=MUnknown);
    bool ArrayFromFile(std::vector<std::vector<double> >& values,std::string parameter=nullstring,
		       std::string filename=noinputtag,MatrixType tempmtype=MUnknown);
    bool ArrayFromFile(std::vector<std::vector<std::string> >& values,std::string parameter=nullstring,
		       std::string filename=noinputtag,MatrixType tempmtype=MUnknown);

    bool ArrayFromString(std::vector<std::vector<int> >& values,std::string parameter=nullstring,
			 std::string inputstring=noinputtag,MatrixType tempmtype=MUnknown);
    bool ArrayFromString(std::vector<std::vector<long int> >& values,std::string parameter=nullstring,
			 std::string inputstring=noinputtag,MatrixType tempmtype=MUnknown);
    bool ArrayFromString(std::vector<std::vector<float> >& values,std::string parameter=nullstring,
			 std::string inputstring=noinputtag,MatrixType tempmtype=MUnknown);
    bool ArrayFromString(std::vector<std::vector<double> >& values,std::string parameter=nullstring,
			 std::string inputstring=noinputtag,MatrixType tempmtype=MUnknown);
    bool ArrayFromString(std::vector<std::vector<std::string> >& values,std::string parameter=nullstring,
			 std::string inputstring=noinputtag,MatrixType tempmtype=MUnknown);

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
    the default array seperator ";" 
    and the default comment indicator "#".
  */

  /*!
    \fn Data_Reader::Data_Reader(std::string,std::string,std::string)
    \brief Special constructor

    Initializes the Data_Reader with the given cut, 
    the given array seperator 
    and the given comment indicator. 
  */

  /*!
    \fn Data_Reader::ReadFromFile(std::string = nullstring, std::string = noinputtag);
    \brief Reads a single Read_Type from a file

    ReadFromFile reads a single Read_Type from a file.
    The input file is specified as an absolute or relative path in the first argument.
    The second argument specifies a pattern that is used to tag the input value.
    A newline terminates the search for the value even if the pattern has been found.
  */

  /*!
    \fn Data_Reader::ReadFromString(std::string = nullstring, std::string = noinputtag);
    \brief Reads a single Read_Type from a string

    ReadFromString reads a single Read_Type from a string.
    The input string is specified in the first argument.
    The second argument specifies a pattern that is used to tag the input value.
  */

  /*!
    \fn Data_Reader::VectorFromFile(std::string = nullstring, std::string = noinputtag, 
                                    std::vector<std::string> = std::vector<std::string>(0));
    \brief Reads a std::vector<Read_Type> from a file

    VectorFromFile reads a std::vector<Read_Type> from a file.
    The input file is specified as an absolute or relative path in the first argument.
    The second argument specifies a pattern that is used to tag the input values.
    A newline terminates the search for new values even if the pattern has been found.
    Data_Reader iterates over all lines in the file so that multiple values will
    be succesively attached to the vector.
  */

  /*!
    \fn Data_Reader::VectorFromString(std::string = nullstring, std::string = noinputtag, 
                                      std::vector<std::string> = std::vector<std::string>(0));
    \brief Reads a std::vector<Read_Type> from a string

    VectorFromString reads a std::vector<Read_Type> from a string.
    The input string is specified in the first argument.
    The second argument specifies a pattern that is used to tag the input values.
  */

  /*!
    \fn Data_Reader::ArrayFromFile(std::string = nullstring, std::string = noinputtag, 
				  std::vector<std::string> = std::vector<std::string>(0));
    \brief Reads a std::vector<std::vector<Read_Type>> from a file

    ArrayFromFile reads a std::vector<std::vector<Read_Type>> from a file.
    The input file is specified as an absolute or relative path in the first argument.
    The second argument specifies a pattern that is used to tag the input values.
    Different arrays of values have to be seperated by a seperator ( the default for this is ';' )
    A newline terminates the search for new values even if the pattern has been found.
    Data_Reader iterates over all lines in the file so that multiple values will
    be succesively attached to the vectors.
  */

  /*!
    \fn Data_Reader::ArrayFromString(std::string = nullstring, std::string = noinputtag, 
				    std::vector<std::string> = std::vector<std::string>(0));
    \brief Reads a std::vector<std::vector<Read_Type>> from a string

    ArrayFromString reads a std::vector<std::vector<Read_Type>> from a string
    The input string is specified in the first argument.
    The second argument specifies a pattern that is used to tag the input values.
    Different arrays of values have to be seperated by a seperator ( the default for this is ';' )
  */

} // end of namespace ATOOLS

#endif
