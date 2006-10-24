#ifndef Data_Writer_H
#define Data_Writer_H

#include "Read_Write_Base.H"

namespace ATOOLS {

  class Data_Writer: public Read_Write_Base {
  public:

    // constructors
    Data_Writer();
    Data_Writer(const std::string cut,const std::string wordsep,
		const std::string linesep,const std::string comment);

    // member functions
    bool WriteComment(std::string comment=nullstring,
		      unsigned int tagreference=0,
		      bool endline=true,std::string tempfname=nullstring);
    bool WriteComment(std::vector<std::string> comments=
		      std::vector<std::string>(0),
		      unsigned int tagreference=0,bool endline=true,
		      std::string tempfname=nullstring);

    template <class Write_Type >
    bool WriteToFile(Write_Type value,
		     std::string tag=nullstring,bool endline=true,
		     std::string tempfname=nullstring,int precision=6);

    template <class Write_Type >
    bool VectorToFile(std::vector<Write_Type> &values,
		      std::string tag=nullstring,bool endline=true,
		      std::string tempfname=nullstring,
		      vtc::code=vtc::unknown,int precision=6);

    template <class Write_Type >
    bool MatrixToFile(std::vector<std::vector<Write_Type> > &values,
		      std::string tag=nullstring,bool endline=true,
		      std::string tempfname=nullstring,
		      mtc::code=mtc::unknown,int precision=6);

    template <class Write_Type >
    inline Data_Writer &operator<<(Write_Type value) 
    { WriteToFile(value,nullstring,false); return *this; }

    template <class Write_Type >
    inline Data_Writer &operator<<(std::vector<Write_Type> &value) 
    { VectorToFile(value,nullstring,false); return *this; }

    template <class Write_Type > inline Data_Writer &
    operator<<(std::vector<std::vector<Write_Type> > &value) 
    { MatrixToFile(value,nullstring,false); return *this; }

  }; // end of class Data_Writer

} // end of namespace ATOOLS

#endif
