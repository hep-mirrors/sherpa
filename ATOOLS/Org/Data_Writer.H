#ifndef ATOOLS_Org_Data_Writer_H
#define ATOOLS_Org_Data_Writer_H

#include "ATOOLS/Org/Read_Write_Base.H"

namespace ATOOLS {

  class Data_Writer: public Read_Write_Base {
  public:

    // constructors
    Data_Writer();
    Data_Writer(const std::string &wordsep,const std::string &linesep,
		const std::string &comment,const std::string &ignore="");

    // destructor
    ~Data_Writer();

    // member functions
    bool WriteComment(std::string comment=nullstring,
		      unsigned int tagreference=0,
		      bool endline=true,std::string tempfname=nullstring);
    bool WriteComment(std::vector<std::string> comments=
		      std::vector<std::string>(0),
		      unsigned int tagreference=0,bool endline=true,
		      std::string tempfname=nullstring);

    template <class Write_Type >
    bool WriteToFile(const Write_Type &value,
		     std::string tag=nullstring,bool endline=true,
		     std::string tempfname=nullstring,int precision=6);

    template <class Write_Type >
    bool VectorToFile(const std::vector<Write_Type> &values,
		      std::string tag=nullstring,bool endline=true,
		      std::string tempfname=nullstring,
		      vtc::code=vtc::unknown,int precision=6);

    template <class Write_Type >
    bool MatrixToFile(const std::vector<std::vector<Write_Type> > &values,
		      std::string tag=nullstring,bool endline=true,
		      std::string tempfname=nullstring,
		      mtc::code=mtc::unknown,int precision=6);

  }; // end of class Data_Writer

} // end of namespace ATOOLS

#endif
