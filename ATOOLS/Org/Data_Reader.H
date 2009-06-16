#ifndef ATOOLS_Org_Data_Reader_H
#define ATOOLS_Org_Data_Reader_H

#include "ATOOLS/Org/Read_Write_Base.H"

namespace ATOOLS {

  class Data_Reader: public Read_Write_Base {
  protected:
    
    std::string m_string;

    template <class Read_Type > Read_Type 
    Convert(std::string cur) const;

    template <class Read_Type > Read_Type 
    ReadValue(const std::string &parameter,const size_t &file);
    template <class Read_Type > std::vector<Read_Type> 
    ReadVector(const std::string &parameter,const size_t &file);
    template <class Read_Type > std::vector< std::vector<Read_Type> > 
    ReadMatrix(const std::string &parameter,const size_t &file);

  public:

    // constructors
    Data_Reader();
    Data_Reader(const std::string &wordsep,const std::string &linesep,
		const std::string &comment,const std::string &ignore="");

    // destructor
    ~Data_Reader();

    // member functions
    template <class Read_Type >
    bool ReadFromFile(Read_Type& value,std::string parameter=nullstring);
    template <class Read_Type >
    bool ReadFromString(Read_Type& value,std::string parameter=nullstring);

    template <class Read_Type >
    bool VectorFromFile(std::vector<Read_Type>& values,
			std::string parameter=nullstring);
    template <class Read_Type >
    bool VectorFromString(std::vector<Read_Type>& values,
			  std::string parameter=nullstring);

    template <class Read_Type >
    bool MatrixFromFile(std::vector<std::vector<Read_Type> >& values,
			std::string parameter=nullstring);
    template <class Read_Type >
    bool MatrixFromString(std::vector<std::vector<Read_Type> >& values,
			  std::string parameter=nullstring);

    void SetString(const std::string string);

    // inline functions
    inline const std::string String() const
    { return m_string; }

    template <class Read_Type> inline Read_Type
    GetValue(const std::string &parameter,const Read_Type &def)
    { 
      Read_Type val; 
      if (!ReadFromFile(val,parameter)) return def; 
      return val; 
    }

  }; // end of class Data_Reader

} // end of namespace ATOOLS

#endif
