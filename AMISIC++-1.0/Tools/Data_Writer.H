#ifndef Data_Writer_H
#define Data_Writer_H

#include "Read_Write_Base.H"

namespace ATOOLS {

  class Data_Writer: public Read_Write_Base {
  private:

    bool InitFile(std::string tempfname);

  protected:

    template <class Write_Type >
    bool M_WriteToFile(Write_Type value,std::string tag,
		       bool endline,std::string tempfname,int precision);
    template <class Write_Type >
    bool M_VectorToFile(std::vector<Write_Type> values,std::string tag,
			bool endline,std::string tempfname,
			VectorTypeID tempvtype,int precision);
    template <class Write_Type >
    bool M_MatrixToFile(std::vector<std::vector<Write_Type> > values,std::string tag,
		       bool endline,std::string tempfname,
		       MatrixTypeID tempmtype,int precision);

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

    bool WriteToFile(int value,std::string tag=nullstring,bool endline=true,
		     std::string tempfname=nullstring,int precision=6);
    bool WriteToFile(long int value,std::string tag=nullstring,bool endline=true,
		     std::string tempfname=nullstring,int precision=6);
    bool WriteToFile(float value,std::string tag=nullstring,bool endline=true,
		     std::string tempfname=nullstring,int precision=6);
    bool WriteToFile(double value,std::string tag=nullstring,bool endline=true,
		     std::string tempfname=nullstring,int precision=6);
    bool WriteToFile(std::string value,std::string tag=nullstring,bool endline=true,
		     std::string tempfname=nullstring,int precision=6);

    bool VectorToFile(std::vector<int>& values,std::string tag=nullstring,bool endline=true,
		      std::string tempfname=nullstring,VectorTypeID=VUnknown,int precision=6);
    bool VectorToFile(std::vector<long int>& values,std::string tag=nullstring,bool endline=true,
		      std::string tempfname=nullstring,VectorTypeID=VUnknown,int precision=6);
    bool VectorToFile(std::vector<float>& values,std::string tag=nullstring,bool endline=true,
		      std::string tempfname=nullstring,VectorTypeID=VUnknown,int precision=6);
    bool VectorToFile(std::vector<double>& values,std::string tag=nullstring,bool endline=true,
		      std::string tempfname=nullstring,VectorTypeID=VUnknown,int precision=6);
    bool VectorToFile(std::vector<std::string>& values,std::string tag=nullstring,bool endline=true,
		      std::string tempfname=nullstring,VectorTypeID=VUnknown,int precision=6);

    bool MatrixToFile(std::vector<std::vector<int> >& values,std::string tag=nullstring,
		     bool endline=true,std::string tempfname=nullstring,
		     MatrixTypeID=MUnknown,int precision=6);
    bool MatrixToFile(std::vector<std::vector<long int> >& values,std::string tag=nullstring,
		     bool endline=true,std::string tempfname=nullstring,
		     MatrixTypeID=MUnknown,int precision=6);
    bool MatrixToFile(std::vector<std::vector<float> >& values,std::string tag=nullstring,
		     bool endline=true,std::string tempfname=nullstring,
		     MatrixTypeID=MUnknown,int precision=6);
    bool MatrixToFile(std::vector<std::vector<double> >& values,std::string tag=nullstring,
		     bool endline=true,std::string tempfname=nullstring,
		     MatrixTypeID=MUnknown,int precision=6);
    bool MatrixToFile(std::vector<std::vector<std::string> >& values,std::string tag=nullstring,
		     bool endline=true,std::string tempfname=nullstring,
		     MatrixTypeID=MUnknown,int precision=6);

  }; // end of class Data_Writer

  //////////////////
  // Doxygen part //
  //////////////////

} // end of namespace ATOOLS

#endif
