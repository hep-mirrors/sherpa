#ifndef Data_Writer_C
#define Data_Writer_C

#include "Data_Writer.H"
#include "Message.H"

#include "My_Filestream.H"
#include "My_Stringstream.H"

#ifdef DEBUG__Data_Writer
#include "My_IO_Stream.H"
#endif

namespace ATOOLS {

  Data_Writer::Data_Writer(): 
    Read_Write_Base(0,1) {}

  Data_Writer::Data_Writer(const std::string _m_cut,
			   const std::string _m_seperator,
			   const  std::string _m_comment):
    Read_Write_Base(0,1,_m_cut,_m_seperator,_m_comment) {}
  
  bool Data_Writer::InitFile(std::string tempfname) 
  {
    if (tempfname!=nullstring) SetOutputFile(tempfname);
    if (OutputFile()==nullstring) {
      ATOOLS::msg.Error()<<"Data_Writer: Error in InitFile(..)! No output file specified."<<std::endl
			 <<"             Ignore command."<<std::endl;
      return false;
    }
    if (!OpenOutFile()) {
      ATOOLS::msg.Error()<<"Data_Writer: Error in WriteComment(..)! Could not open output file."<<std::endl
			 <<"             Ignore command."<<std::endl;
      return false;
    }
    return true;
  }

  bool Data_Writer::WriteComment(std::string comment,unsigned int tagreference,
				 bool endline,std::string tempfname)
  {
    std::string tag;
#ifdef DEBUG__Data_Writer
    std::cout<<"Data_Writer::WriteComment("<<comment<<","<<tagreference<<","<<tempfname<<")"<<std::endl;
    std::cout<<"   comment tags are ";
    for (unsigned int i=0;i<Comment().size();std::cout<<"'"<<Comment()[i++]<<"' ");
    std::cout<<std::endl;
#endif
    if (tagreference>=Comment().size()) tag=defaultcom;
    else tag=Comment()[tagreference];
    if (!InitFile(tempfname)) return false;
    if (tag!=nullstring) {
      if (Blank().size()>0) *OutFile()<<tag<<(char)Blank()[0]<<comment;
      else *OutFile()<<tag<<comment;
    }
    else {
      *OutFile()<<comment;
    }
    if (endline) *OutFile()<<std::endl;
    CloseOutFile();
    return true;
  }
  
  bool Data_Writer::WriteComment(std::vector<std::string> comments,unsigned int tagreference,
				 bool endline,std::string tempfname)
  {
    for (unsigned int i=0;i<comments.size();++i) {
      if (!WriteComment(comments[i],tagreference,endline,tempfname)) return false;
    }
    return true;
  }

  template <class Write_Type>  
  bool Data_Writer::M_WriteToFile(Write_Type value,std::string tag,bool endline,
				  std::string tempfname,int precision)
  {
    if (!InitFile(tempfname)) return false;
    const std::ios_base::fmtflags defaultflags=OutFile()->flags();
    OutFile()->precision(precision);
    if (tag!=nullstring) {
      if (Blank().size()>0) *OutFile()<<tag<<(char)Blank()[0]<<value;
      else *OutFile()<<tag<<value;
    }
    else {
      *OutFile()<<value;
    }
    if (endline) *OutFile()<<std::endl;
    OutFile()->flags(defaultflags);
    CloseOutFile();
    return true;    
  }

  template <class Write_Type>  
  bool Data_Writer::M_VectorToFile(std::vector<Write_Type> values,std::string tag,
				   bool endline,std::string tempfname,
				   VectorTypeID tempvtype,int precision)
  {
    if (!InitFile(tempfname)) return false;
    if (tempvtype==VUnknown) tempvtype=VectorType();
    switch (tempvtype) {
    case VHorizontal:
      if (values.size()>0) M_WriteToFile<Write_Type>(values[0],tag,false,tempfname,precision);
      for (unsigned int i=1;i<values.size();++i) {
	M_WriteToFile<Write_Type>(values[i],"",false,tempfname,precision);
	if (Blank().size()>0) *OutFile()<<(char)Blank()[0];
      }
      if (endline) *OutFile()<<std::endl;
      break;
    case VVertical:
    default:
      for (unsigned int i=0;i<values.size();++i) M_WriteToFile<Write_Type>(values[i],tag,true,
									   tempfname,precision);
      break;
    }
    CloseOutFile();
    return true;    
  }

  template <class Write_Type>  
  bool Data_Writer::M_MatrixToFile(std::vector<std::vector<Write_Type> > values,std::string tag,
				  bool endline,std::string tempfname,
				  MatrixTypeID tempmtype,int precision)
  {
    std::vector<std::vector<Write_Type> > tempvalues;
    if (!InitFile(tempfname)) return false;
    if (tempmtype==MUnknown) tempmtype=MatrixType();
    switch (tempmtype) {
    case MTransposed:
      for (unsigned int i=0;i<values.size();++i) {
	M_VectorToFile<Write_Type>(values[i],tag,true,tempfname,VHorizontal,precision);
      }
      break;
    case MNormal:
    default:
      if (values.size()!=0) {
	for (unsigned int i=0;i<values[0].size();++i) {
	  tempvalues.push_back(std::vector<Write_Type>(values.size()));
	}
      }
      for (unsigned int j=0;j<values.size();++j) {
	for (unsigned int k=0;k<values[j].size();++k) {
	  tempvalues[k][j]=values[j][k];
	}
      }
      for (unsigned int i=0;i<tempvalues.size();++i) {
	M_VectorToFile<Write_Type>(tempvalues[i],tag,true,tempfname,VHorizontal,precision);
      }
      break;
    }
    CloseOutFile();
    return true;    
  }

  bool Data_Writer::WriteToFile(int value,std::string tag,bool endline,
				std::string tempfname,int precision)
  { return M_WriteToFile<int>(value,tag,endline,tempfname,precision); }
  bool Data_Writer::WriteToFile(long int value,std::string tag,bool endline,
				std::string tempfname,int precision)
  { return M_WriteToFile<long int>(value,tag,endline,tempfname,precision); }
  bool Data_Writer::WriteToFile(float value,std::string tag,bool endline,
				std::string tempfname,int precision)
  { return M_WriteToFile<float>(value,tag,endline,tempfname,precision); }
  bool Data_Writer::WriteToFile(double value,std::string tag,bool endline,
				std::string tempfname,int precision)
  { return M_WriteToFile<double>(value,tag,endline,tempfname,precision); }
  bool Data_Writer::WriteToFile(std::string value,std::string tag,bool endline,
				std::string tempfname,int precision)
  { return M_WriteToFile<std::string>(value,tag,endline,tempfname,precision); }

  bool Data_Writer::VectorToFile(std::vector<int>& values,std::string tag,
				 bool endline,std::string tempfname,
				 VectorTypeID tempvtype,int precision)
  { return M_VectorToFile<int>(values,tag,endline,tempfname,tempvtype,precision); }
  bool Data_Writer::VectorToFile(std::vector<long int>& values,std::string tag,
				 bool endline,std::string tempfname,
				 VectorTypeID tempvtype,int precision)
  { return M_VectorToFile<long int>(values,tag,endline,tempfname,tempvtype,precision); }
  bool Data_Writer::VectorToFile(std::vector<float>& values,std::string tag,
				 bool endline,std::string tempfname,
				 VectorTypeID tempvtype,int precision)
  { return M_VectorToFile<float>(values,tag,endline,tempfname,tempvtype,precision); }
  bool Data_Writer::VectorToFile(std::vector<double>& values,std::string tag,
				 bool endline,std::string tempfname,
				 VectorTypeID tempvtype,int precision)
  { return M_VectorToFile<double>(values,tag,endline,tempfname,tempvtype,precision); }
  bool Data_Writer::VectorToFile(std::vector<std::string>& values,std::string tag,
				 bool endline,std::string tempfname,
				 VectorTypeID tempvtype,int precision)
  { return M_VectorToFile<std::string>(values,tag,endline,tempfname,tempvtype,precision); }

  bool Data_Writer::MatrixToFile(std::vector<std::vector<int> >& values,std::string tag,
				bool endline,std::string tempfname,
				MatrixTypeID tempmtype,int precision)
  { return M_MatrixToFile<int>(values,tag,endline,tempfname,tempmtype,precision); }
  bool Data_Writer::MatrixToFile(std::vector<std::vector<long int> >& values,std::string tag,
				bool endline,std::string tempfname,
				MatrixTypeID tempmtype,int precision)
  { return M_MatrixToFile<long int>(values,tag,endline,tempfname,tempmtype,precision); }
  bool Data_Writer::MatrixToFile(std::vector<std::vector<float> >& values,std::string tag,
				bool endline,std::string tempfname,
				MatrixTypeID tempmtype,int precision)
  { return M_MatrixToFile<float>(values,tag,endline,tempfname,tempmtype,precision); }
  bool Data_Writer::MatrixToFile(std::vector<std::vector<double> >& values,std::string tag,
				bool endline,std::string tempfname,
				MatrixTypeID tempmtype,int precision)
  { return M_MatrixToFile<double>(values,tag,endline,tempfname,tempmtype,precision); }
  bool Data_Writer::MatrixToFile(std::vector<std::vector<std::string> >& values,std::string tag,
				bool endline,std::string tempfname,
				MatrixTypeID tempmtype,int precision)
  { return M_MatrixToFile<std::string>(values,tag,endline,tempfname,tempmtype,precision); }

} // end of namespace ATOOLS

#endif

