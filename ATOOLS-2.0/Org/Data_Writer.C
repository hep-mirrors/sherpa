#include "Data_Writer.H"

#include "Message.H"
#include "MyStrStream.H"
#ifdef DEBUG__Data_Writer
#include <iostream>
#endif

using namespace ATOOLS;

Data_Writer::Data_Writer(): 
  Read_Write_Base(0,1) 
{
  SetOutFileMode(fom::permanent);
}

Data_Writer::Data_Writer(const std::string cut,
			 const std::string seperator,
			 const  std::string comment):
  Read_Write_Base(0,1,cut,seperator,comment) 
{
  SetOutFileMode(fom::permanent);
}

bool Data_Writer::WriteComment(std::string comment,
			       unsigned int tagreference,
			       bool endline,std::string tempfname)
{
  std::string tag;
  if (tagreference>=Comment().size()) tag=defaultcom;
  else tag=Comment()[tagreference];
  if (tempfname!=nullstring) SetOutputFile(tempfname);
  if (!OpenOutFile()) return false;
  if (tag!=nullstring) *OutFile()<<tag;
  if (Blank().size()>0) *OutFile()<<(char)Blank()[0];
  *OutFile()<<comment;
  if (endline) *OutFile()<<std::endl;
  CloseOutFile();
  return true;
}

bool Data_Writer::WriteComment(std::vector<std::string> comments,
			       unsigned int tagreference,
			       bool endline,std::string tempfname)
{
  for (unsigned int i=0;i<comments.size();++i) {
    if (!WriteComment(comments[i],tagreference,endline,tempfname)) 
      return false;
  }
  return true;
}

template <class Write_Type>  
bool Data_Writer::WriteToFile(Write_Type value,std::string tag,bool endline,
			      std::string tempfname,int precision)
{
  if (tempfname!=nullstring) SetOutputFile(tempfname);
  if (!OpenOutFile()) return false;
#ifdef __GNUC__
#if __GNUC__ > 2
  const std::ios_base::fmtflags defaultflags=OutFile()->flags();
#else
  const std::ios::fmtflags defaultflags=OutFile()->flags();
#endif
#else
  const std::ios::fmtflags defaultflags=OutFile()->flags();
#endif
  OutFile()->precision(precision);
  if (tag!=nullstring) *OutFile()<<tag;
  if (Blank().size()>0) *OutFile()<<(char)Blank()[0];
  *OutFile()<<value;
  if (endline) *OutFile()<<std::endl;
  OutFile()->flags(defaultflags);
  CloseOutFile();
  return true;    
}

template <class Write_Type>  
bool Data_Writer::VectorToFile(std::vector<Write_Type> &values,
			       std::string tag,
			       bool endline,std::string tempfname,
			       vtc::code tempvtype,int precision)
{
  if (tempfname!=nullstring) SetOutputFile(tempfname);
  if (!OpenOutFile()) return false;
  if (tempvtype==vtc::unknown) tempvtype=VectorType();
  switch (tempvtype) {
  case vtc::horizontal:
    if (values.size()>0) {
      WriteToFile<Write_Type>(values[0],tag,false,tempfname,precision);
      if (Blank().size()>0) *OutFile()<<(char)Blank()[0];
    }
    for (unsigned int i=1;i<values.size();++i) {
      WriteToFile<Write_Type>(values[i],"",false,tempfname,precision);
      if (Blank().size()>0) *OutFile()<<(char)Blank()[0];
    }
    if (endline) *OutFile()<<std::endl;
    break;
  case vtc::vertical:
  default:
    for (unsigned int i=0;i<values.size();++i) 
      WriteToFile<Write_Type>(values[i],tag,true,tempfname,precision);
    break;
  }
  CloseOutFile();
  return true;    
}

template <class Write_Type>  
bool Data_Writer::MatrixToFile(std::vector<std::vector<Write_Type> > 
			       &values,std::string tag,
			       bool endline,std::string tempfname,
			       mtc::code tempmtype,int precision)
{
  std::vector<std::vector<Write_Type> > tempvalues;
  if (tempfname!=nullstring) SetOutputFile(tempfname);
  if (!OpenOutFile()) return false;
  if (tempmtype==mtc::unknown) tempmtype=MatrixType();
  switch (tempmtype) {
  case mtc::transposed:
    for (unsigned int i=0;i<values.size();++i) {
      VectorToFile<Write_Type>(values[i],tag,true,tempfname,
			       vtc::horizontal,precision);
    }
    break;
  case mtc::normal:
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
      VectorToFile<Write_Type>(tempvalues[i],tag,true,tempfname,
			       vtc::horizontal,precision);
    }
    break;
  }
  CloseOutFile();
  return true;    
}

template bool Data_Writer::WriteToFile<int>
(int,std::string,bool,std::string,int);
template bool Data_Writer::WriteToFile<unsigned int>
(unsigned int,std::string,bool,std::string,int);
template bool Data_Writer::WriteToFile<long int>
(long int,std::string,bool,std::string,int);
template bool Data_Writer::WriteToFile<float>
(float,std::string,bool,std::string,int);
template bool Data_Writer::WriteToFile<double>
(double,std::string,bool,std::string,int);
template bool Data_Writer::WriteToFile<std::string>
(std::string,std::string,bool,std::string,int);
template bool Data_Writer::WriteToFile<const char*>
(const char*,std::string,bool,std::string,int);

template bool Data_Writer::VectorToFile<int>
(std::vector<int>&,std::string,bool,std::string,vtc::code,int);
template bool Data_Writer::VectorToFile<unsigned int>
(std::vector<unsigned int>&,std::string,bool,std::string,vtc::code,int);
template bool Data_Writer::VectorToFile<long int>
(std::vector<long int>&,std::string,bool,std::string,vtc::code,int);
template bool Data_Writer::VectorToFile<float>
(std::vector<float>&,std::string,bool,std::string,vtc::code,int);
template bool Data_Writer::VectorToFile<double>
(std::vector<double>&,std::string,bool,std::string,vtc::code,int);
template bool Data_Writer::VectorToFile<std::string>
(std::vector<std::string>&,std::string,bool,std::string,vtc::code,int);
template bool Data_Writer::VectorToFile<const char*>
(std::vector<const char*>&,std::string,bool,std::string,vtc::code,int);

template bool Data_Writer::MatrixToFile<int>
(std::vector<std::vector<int> >&,std::string,bool,
 std::string,mtc::code,int);
template bool Data_Writer::MatrixToFile<unsigned int>
(std::vector<std::vector<unsigned int> >&,std::string,bool,
 std::string,mtc::code,int);
template bool Data_Writer::MatrixToFile<long int>
(std::vector<std::vector<long int> >&,std::string,bool,
 std::string,mtc::code,int);
template bool Data_Writer::MatrixToFile<float>
(std::vector<std::vector<float> >&,std::string,bool,
 std::string,mtc::code,int);
template bool Data_Writer::MatrixToFile<double>
(std::vector<std::vector<double> >&,std::string,bool,
 std::string,mtc::code,int);
template bool Data_Writer::MatrixToFile<std::string>
(std::vector<std::vector<std::string> >&,std::string,bool,
 std::string,mtc::code,int);
template bool Data_Writer::MatrixToFile<const char*>
(std::vector<std::vector<const char*> >&,std::string,bool,
 std::string,mtc::code,int);

