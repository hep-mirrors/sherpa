#ifndef Data_Writer_C
#define Data_Writer_C

#include "Data_Writer.H"
#include "Message.H"
#include "Type.H"

#include "My_Filestream.H"
#include "My_Stringstream.H"

#ifdef DEBUG__Data_Writer
#include "My_IO_Stream.H"
#endif

namespace ATOOLS {

  Data_Writer::Data_Writer(): Read_Write_Base() 
  {
    SetOpenMode(Permanent);
  }

  Data_Writer::Data_Writer(std::string _m_cut, std::string _m_seperator, std::string _m_comment):
    Read_Write_Base(_m_cut,_m_seperator,_m_comment) 
  {
    SetOpenMode(Permanent);
  }
  
  Data_Writer::Data_Writer(const char *_m_cut, const char *_m_seperator, const char *_m_comment):
    Read_Write_Base(_m_cut,_m_seperator,_m_comment) 
  {
    SetOpenMode(Permanent);
  }
  
  bool Data_Writer::InitFile(std::string tempfname) 
  {
    if (tempfname!=nullstring) SetFileName(tempfname);
    if (M_FileName()==nullstring) {
      ATOOLS::msg.Error()<<"Data_Writer: Error in InitFile(..)! No output file specified."<<std::endl
			 <<"             Ignore command."<<std::endl;
      return false;
    }
    if (!OpenFile(M_FileName(),std::ios_base::out)) {
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
    int blank=defaultblank;
#ifdef DEBUG__Data_Writer
    std::cout<<"Data_Writer::WriteComment("<<comment<<","<<tagreference<<","<<tempfname<<")"<<std::endl;
    std::cout<<"   comment tags are ";
    for (unsigned int i=0;i<M_Comment().size();std::cout<<"'"<<M_Comment()[i++]<<"' ");
    std::cout<<std::endl;
#endif
    if (tagreference>=M_Comment().size()) tag=defaultcom;
    else tag=M_Comment()[tagreference];
    if (!InitFile(tempfname)) return false;
    if (M_Blank().size()>0) blank=M_Blank()[0];
    *M_File()<<tag<<(char)blank<<comment;
    if (endline) *M_File()<<std::endl;
    CloseFile();
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
    const std::ios_base::fmtflags defaultflags=M_File()->flags();
    M_File()->precision(precision);
    int blank=defaultblank;
    if (M_Blank().size()>0) blank=M_Blank()[0];
    *M_File()<<tag<<(char)blank<<value;
    if (endline) *M_File()<<std::endl;
    M_File()->flags(defaultflags);
    CloseFile();
    return true;    
  }

  template <class Write_Type>  
  bool Data_Writer::M_VectorToFile(std::vector<Write_Type> values,std::string tag,
				   bool endline,std::string tempfname,
				   VectorType tempvtype,int precision)
  {
    int blank=defaultblank;
    if (!InitFile(tempfname)) return false;
    if (tempvtype==VUnknown) tempvtype=M_VectorType();
    if (M_Blank().size()>0) blank=M_Blank()[0];
    switch (tempvtype) {
    case VHorizontal:
      if (values.size()>0) M_WriteToFile<Write_Type>(values[0],tag,false,tempfname,precision);
      for (unsigned int i=1;i<values.size();++i) {
	M_WriteToFile<Write_Type>(values[i],"",false,tempfname,precision);
	*M_File()<<(char)blank;
      }
      if (endline) *M_File()<<std::endl;
      break;
    case VVertical:
    default:
      for (unsigned int i=0;i<values.size();++i) M_WriteToFile<Write_Type>(values[i],tag,true,
									   tempfname,precision);
      break;
    }
    CloseFile();
    return true;    
  }

  template <class Write_Type>  
  bool Data_Writer::M_ArrayToFile(std::vector<std::vector<Write_Type> > values,std::string tag,
				  bool endline,std::string tempfname,
				  MatrixType tempmtype,int precision)
  {
    std::vector<std::vector<Write_Type> > tempvalues;
    if (!InitFile(tempfname)) return false;
    if (tempmtype==MUnknown) tempmtype=M_MatrixType();
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
    CloseFile();
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
				 VectorType tempvtype,int precision)
  { return M_VectorToFile<int>(values,tag,endline,tempfname,tempvtype,precision); }
  bool Data_Writer::VectorToFile(std::vector<long int>& values,std::string tag,
				 bool endline,std::string tempfname,
				 VectorType tempvtype,int precision)
  { return M_VectorToFile<long int>(values,tag,endline,tempfname,tempvtype,precision); }
  bool Data_Writer::VectorToFile(std::vector<float>& values,std::string tag,
				 bool endline,std::string tempfname,
				 VectorType tempvtype,int precision)
  { return M_VectorToFile<float>(values,tag,endline,tempfname,tempvtype,precision); }
  bool Data_Writer::VectorToFile(std::vector<double>& values,std::string tag,
				 bool endline,std::string tempfname,
				 VectorType tempvtype,int precision)
  { return M_VectorToFile<double>(values,tag,endline,tempfname,tempvtype,precision); }
  bool Data_Writer::VectorToFile(std::vector<std::string>& values,std::string tag,
				 bool endline,std::string tempfname,
				 VectorType tempvtype,int precision)
  { return M_VectorToFile<std::string>(values,tag,endline,tempfname,tempvtype,precision); }

  bool Data_Writer::ArrayToFile(std::vector<std::vector<int> >& values,std::string tag,
				bool endline,std::string tempfname,
				MatrixType tempmtype,int precision)
  { return M_ArrayToFile<int>(values,tag,endline,tempfname,tempmtype,precision); }
  bool Data_Writer::ArrayToFile(std::vector<std::vector<long int> >& values,std::string tag,
				bool endline,std::string tempfname,
				MatrixType tempmtype,int precision)
  { return M_ArrayToFile<long int>(values,tag,endline,tempfname,tempmtype,precision); }
  bool Data_Writer::ArrayToFile(std::vector<std::vector<float> >& values,std::string tag,
				bool endline,std::string tempfname,
				MatrixType tempmtype,int precision)
  { return M_ArrayToFile<float>(values,tag,endline,tempfname,tempmtype,precision); }
  bool Data_Writer::ArrayToFile(std::vector<std::vector<double> >& values,std::string tag,
				bool endline,std::string tempfname,
				MatrixType tempmtype,int precision)
  { return M_ArrayToFile<double>(values,tag,endline,tempfname,tempmtype,precision); }
  bool Data_Writer::ArrayToFile(std::vector<std::vector<std::string> >& values,std::string tag,
				bool endline,std::string tempfname,
				MatrixType tempmtype,int precision)
  { return M_ArrayToFile<std::string>(values,tag,endline,tempfname,tempmtype,precision); }

} // end of namespace ATOOLS

#endif

