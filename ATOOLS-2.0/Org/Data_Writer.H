#ifndef Data_Writer_H
#define Data_Writer_H

#include "Read_Write_Base.H"

namespace ATOOLS {

  class Data_Writer: public Read_Write_Base {
  public:

    // constructors
    Data_Writer();
    Data_Writer(const std::string _m_cut,
		const std::string _m_seperator,
		const std::string _m_comment);

    // member functions
    bool WriteComment(std::string comment=nullstring,unsigned int tagreference=0,
		      bool endline=true,std::string tempfname=nullstring);
    bool WriteComment(std::vector<std::string> comments=std::vector<std::string>(0),
		      unsigned int tagreference=0,bool endline=true,
		      std::string tempfname=nullstring);

    template <class Write_Type >
    bool WriteToFile(Write_Type value,std::string tag=nullstring,bool endline=true,
		     std::string tempfname=nullstring,int precision=6);

    template <class Write_Type >
    bool VectorToFile(std::vector<Write_Type> &values,std::string tag=nullstring,bool endline=true,
		      std::string tempfname=nullstring,VectorTypeID=VUnknown,int precision=6);

    template <class Write_Type >
    bool MatrixToFile(std::vector<std::vector<Write_Type> > &values,std::string tag=nullstring,
		      bool endline=true,std::string tempfname=nullstring,
		      MatrixTypeID=MUnknown,int precision=6);

  }; // end of class Data_Writer

  //////////////////
  // Doxygen part //
  //////////////////

} // end of namespace ATOOLS

#endif
