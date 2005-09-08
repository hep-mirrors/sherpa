#ifndef Data_Reader_H
#define Data_Reader_H

#include "Read_Write_Base.H"
#include <map>

namespace ATOOLS {

  class Data_Reader: public Read_Write_Base {
  private:

    void KillComments(std::string& buffer);

    std::string KillBlanks(std::string& buffer);
    std::string HighlightSeparator(std::string& buffer);

  protected:

    template <class Read_Type >
    Read_Type M_ReadFromFile(std::string parameter,
			     std::string filename);
    template <class Read_Type >
    Read_Type M_ReadFromString(std::string parameter,
			       std::string &inputstring);

    template <class Read_Type >
    std::vector<Read_Type> 
    M_VectorFromFile(std::string parameter,std::string filename,
		     vtc::code tempvtype);
    template <class Read_Type >
    std::vector<Read_Type>
    M_VectorFromString(std::string parameter,std::string inputstring,
		       vtc::code tempvtype);

    template <class Read_Type >
    std::vector< std::vector<Read_Type> > 
    M_MatrixFromFile(std::string parameter,std::string filename,
		     mtc::code tempmtype);
    template <class Read_Type >
    std::vector< std::vector<Read_Type> >
    M_MatrixFromString(std::string parameter,std::string inputstring,
		       mtc::code tempmtype);

  public:

    // constructors
    Data_Reader();
    Data_Reader(const std::string _m_cut,
		const std::string _m_separator,
		const std::string _m_comment);

    // member functions
    template <class Read_Type >
    bool ReadFromFile(Read_Type& value,std::string parameter=nullstring,
		      std::string filename=noinputtag);
    template <class Read_Type >
    bool ReadFromString(Read_Type& value,std::string parameter=nullstring,
			std::string inputstring=noinputtag);

    template <class Read_Type >
    bool VectorFromFile(std::vector<Read_Type>& values,
			std::string parameter=nullstring,
			std::string filename=noinputtag,
			vtc::code tempvtype=vtc::unknown);
    template <class Read_Type >
    bool VectorFromString(std::vector<Read_Type>& values,
			  std::string parameter=nullstring,
			  std::string inputstring=noinputtag,
			  vtc::code tempvtype=vtc::unknown);

    template <class Read_Type >
    bool MatrixFromFile(std::vector<std::vector<Read_Type> >& values,
			std::string parameter=nullstring,
			std::string filename=noinputtag,
			mtc::code tempmtype=mtc::unknown);
    template <class Read_Type >
    bool MatrixFromString(std::vector<std::vector<Read_Type> >& values,
			  std::string parameter=nullstring,
			  std::string inputstring=noinputtag,
			  mtc::code tempmtype=mtc::unknown);

    std::string StripEscapes(const std::string &buffer) const;

  }; // end of class Data_Reader

} // end of namespace ATOOLS

#endif
