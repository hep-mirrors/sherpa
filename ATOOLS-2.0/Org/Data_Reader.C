#include "Data_Reader.H"

#include "Message.H"
#include "MyStrStream.H"
#include <typeinfo>
#include <ctype.h>
#ifdef DEBUG__Data_Reader
#include <iostream>
#endif

using namespace ATOOLS;

Data_Reader::Data_Reader(): 
  Read_Write_Base(1,0)
{
  SetInFileMode(Permanent);
}

Data_Reader::Data_Reader(const std::string _m_cut,
			 const std::string _m_separator,
			 const std::string _m_comment):
  Read_Write_Base(1,0,_m_cut,_m_separator,_m_comment)
{
  SetInFileMode(Permanent);
}

void Data_Reader::KillComments(std::string& buffer)
{
#ifdef DEBUG__Data_Reader
  std::cout<<"Data_Reader::KillComments("<<buffer<<")"<<std::endl<<"   comment tags are ";
  for (unsigned int i=0;i<Comment().size();std::cout<<"'"<<Comment()[i++]<<"' ");
  std::cout<<std::endl<<"   ignored tags are ";
  for (unsigned int i=0;i<Ignore().size();std::cout<<"'"<<Ignore()[i++]<<"' ");
  std::cout<<std::endl;
#endif
  size_t pos;
  for (unsigned int i=0;i<Comment().size();++i) {
    size_t next=0;
    while ((pos=buffer.find(Comment()[i],next))!=std::string::npos) {
      if (pos>0 && buffer[pos-1]==Escape()) next=pos+1;
      else buffer=buffer.substr(0,pos);
    }
  }
  KillBlanks(buffer);
  for (unsigned int i=0; i<Ignore().size(); ++i) {
    size_t next=0;
    while ((pos=buffer.find(Ignore()[i],next))!=std::string::npos) {
      if (pos>0 && buffer[pos-1]==Escape()) next=pos+1;
      else {
	buffer=buffer.substr(0,pos)+" "+
	  buffer.substr(pos+Ignore()[i].length());
	KillBlanks(buffer);
      }
    }
  }
#ifdef DEBUG__Data_Reader
  std::cout<<"   returning '"<<buffer<<"'"<<std::endl;
#endif
}

std::string Data_Reader::KillBlanks(std::string& buffer)
{
#ifdef DEBUG__Data_Reader
  std::cout<<"Data_Reader::KillBlanks("<<buffer<<")"<<std::endl;
#endif
  if (buffer==nullstring) return buffer;
  bool hit=true;
  while (hit && buffer.length()>0) { 
    hit=false;
    for (size_t i=0;i<Blank().size();++i) 
      if (int(buffer[0])==Blank()[i]) {
	buffer.erase(0,1); 
	hit=true;
	break;
      }
  }
  hit=true;
  while (hit && buffer.length()>0) { 
    if (buffer.length()>1 && buffer[buffer.length()-1]==Escape()) break;
    hit=false;
    for (size_t i=0;i<Blank().size();++i) 
      if (int(buffer[buffer.length()-1])==Blank()[i]) {
	buffer.erase(buffer.length()-1,1);
	hit=true;
	break;
      }
  }
#ifdef DEBUG__Data_Reader
  std::cout<<"   returning '"<<buffer<<"'"<<std::endl;
#endif
  return buffer;
}

std::string Data_Reader::HighlightSeparator(std::string& buffer)
{
#ifdef DEBUG__Data_Reader
  std::cout<<"Data_Reader::HighlightSeparator("<<buffer<<")"<<std::endl;
#endif
  if (buffer==nullstring) return buffer;
  size_t pos=std::string::npos, next=0, min=pos;
  for (unsigned int j=0; j<Separator().size(); ++j) {
    while (pos!=next &&
	   (next=pos=buffer.find(Separator()[j],next))!=std::string::npos) {
      if (pos>0 && buffer[pos-1]==Escape()) next=pos+1;
      else {
	buffer.insert(pos+1," ",1);
	buffer.insert(pos," ",1);
	if (pos<min) min=pos+1;
      }
    }
  }
#ifdef DEBUG__Data_Reader
  std::cout<<"   returning '"<<buffer<<"'"<<std::endl;
#endif
  return buffer.substr(0,min);
}

std::string Data_Reader::StripEscapes(const std::string &buffer) const
{
#ifdef DEBUG__Data_Reader
  std::cout<<"Data_Reader::StripEscapes("<<buffer<<")"<<std::endl;
#endif
  if (buffer.length()==0) return buffer;
  std::string input=buffer;
  size_t pos, next=0;
  while ((pos=input.find(Escape(),next))!=std::string::npos) {
    input.erase(pos,1);
    if (input.length()>pos && input[pos]==Escape()) next=pos+1;
  }
#ifdef DEBUG__Data_Reader
  std::cout<<"   returning '"<<input<<"'"<<std::endl;
#endif
  return input;
}

template <class Read_Type>
Read_Type Data_Reader::M_ReadFromString(std::string parameter,std::string &inputstring)
{
  size_t pos;
  Read_Type value;
#ifdef DEBUG__Data_Reader
  std::cout<<"Data_Reader::M_ReadFromString("<<parameter<<","<<inputstring<<")"<<std::endl;
#endif
  if (inputstring==noinputtag) {
    if (String()!=nullstring) {
      inputstring=String();
    }
    else {
      msg_Tracking()<<"Data_Reader: No input string specified ! No default available !"<<std::endl
			    <<"   Abort reading."<<std::endl;
      return Default<Read_Type>();
    }
  }
  KillComments(inputstring);
  if (parameter==nullstring) { 
    parameter=std::string("DUMMY_PARAMETER");
    inputstring=parameter+inputstring;
  }
  size_t length=0;
  if(((pos=Find(inputstring,parameter,length))!=std::string::npos)&&
     ((inputstring=inputstring.substr(pos+length)).length()>0)) {
    std::string cur=HighlightSeparator(inputstring);
    std::string tempstring=ReplaceTags(cur);
    if (typeid(value)==typeid(int) || 
	typeid(value)==typeid(unsigned int) ||
	typeid(value)==typeid(float) ||
	typeid(value)==typeid(double)) {
      if (!m_allownans) {
 	if ((pos=tempstring.find("nan"))!=std::string::npos)
 	  tempstring.replace(pos,3,"1");
 	else if ((pos=tempstring.find("inf"))!=std::string::npos) 
 	  tempstring.replace(pos,3,"1");
      }
      if (Interprete()) tempstring=Interpreter()->
	  Interprete(StripEscapes(tempstring));
    }
    value=ATOOLS::ToType<Read_Type>(tempstring);
#ifdef DEBUG__Data_Reader
    std::cout<<"   returning '"<<value<<"'"<<" ( type = "<<typeid(value).name()<<" )"<<std::endl;
#endif
    return value;
  }
#ifdef DEBUG__Data_Reader
  std::cout<<"   returning default value '"<<Default<Read_Type>()<<"'"
	   <<" ( type "<<typeid(value).name()<<" )"<<std::endl;
#endif
  return Default<Read_Type>();
}

template <class Read_Type>
Read_Type Data_Reader::M_ReadFromFile(std::string parameter,std::string filename)
{
  Read_Type value, temp;
  std::string buffer;
  value=Default<Read_Type>();
  if (filename==noinputtag) {
    if (InputFile()!=nullstring) {
      filename=InputFile();
    }
    else {
      msg_Tracking()<<"Data_Reader: No input file specified. "
		    <<"No default available."<<std::endl
		    <<"   Abort reading."<<std::endl;
      return value;
    }
  }
  if(!OpenInFile()) {
    msg.Tracking()<<"Data_Reader: File '"<<filename<<"' not found."<<std::endl;
    return value;
  }
  for (unsigned int i=0;i<FileContent().size();++i) {
    buffer=FileContent()[i];
    if((temp=M_ReadFromString<Read_Type>(parameter,buffer))!=Default<Read_Type>()) value=temp;
  }
  CloseInFile();
  if (value==Default<Read_Type>()) {
    msg_Debugging()<<"Data_Reader: Parameter '"<<parameter
		   <<"' not specified in '"<<filename<<"' !"<<std::endl;
  }
  return value;
}

template <class Read_Type>
std::vector<Read_Type> 
Data_Reader::M_VectorFromString(std::string parameter, std::string inputstring,VectorTypeID tempvtype)
{
  if (tempvtype==VUnknown) tempvtype=VectorType();
  if (tempvtype==VUnknown) tempvtype=VVertical;
  std::string value;
  std::vector<Read_Type> values;
  inputstring = ReplaceTags(inputstring);
#ifdef DEBUG__Data_Reader
  std::cout<<"Data_Reader::M_VectorFromString("<<parameter<<","<<inputstring<<","<<tempvtype<<")"<<std::endl;
#endif
  if (inputstring==noinputtag) {
    if (String()!=nullstring) {
      inputstring=String();
    }
    else {
      msg_Tracking()<<"Data_Reader: No input string specified. "
		    <<"No default available !"<<std::endl
		    <<"   Abort reading."<<std::endl;
      return values;
    }
  }
  value = M_ReadFromString<std::string>(parameter,inputstring);
  while (value != Default<std::string>()) {
#ifdef DEBUG__Data_Reader
    std::cout<<"   newline tags are ";
    for (unsigned int j=0;j<Separator().size();std::cout<<"'"<<Separator()[j++]<<"' ");
    std::cout<<std::endl;
#endif
    bool stop=false;
    for (unsigned int j=0;j<Separator().size();++j) if (value==Separator()[j]) stop=true;
    if (stop) break;
    values.push_back(M_ReadFromString<Read_Type>(nullstring,value));
    if (tempvtype==VVertical) break;
    inputstring=inputstring.substr(inputstring.find(value)+value.length());
    value=M_ReadFromString<std::string>(nullstring,inputstring);
  }
#ifdef DEBUG__Data_Reader
  std::cout<<"   returned values are [ ";
  for (unsigned int j=0;j<values.size();std::cout<<values[j++]<<" ");
  std::cout<<"]"<<std::endl;
#endif
  if (values.size()!=0) return values;
  return values;
}

template <class Read_Type>
std::vector<Read_Type> 
Data_Reader::M_VectorFromFile(std::string parameter, std::string filename,VectorTypeID tempvtype)
{
  if (tempvtype==VUnknown) tempvtype=VectorType();
  if (tempvtype==VUnknown) tempvtype=VVertical;
  std::vector<Read_Type> values, temp;
  if (filename==noinputtag) {
    if (InputFile()!=nullstring) {
      filename=InputFile();
    }
    else {
      msg_Tracking()<<"Data_Reader: No input file specified ! No default available !"<<std::endl
			    <<"   Abort reading."<<std::endl;
      return values;
    }
  }
  if (!OpenInFile()) {
    msg_Tracking()<<"Data_Reader: Error opening "<<filename<<" !"<<std::endl;
    return values;
  }
  for (unsigned int i=0;i<FileContent().size();++i) {
    temp = M_VectorFromString<Read_Type>(parameter,FileContent()[i],tempvtype);
    switch (tempvtype) {
    case VHorizontal:
      if (temp.size()!=0) values=temp;
      break;
    case VUnknown:
    case VVertical:
      for (unsigned int j=0; j<temp.size(); ++j) values.push_back(temp[j]);
      break;
    }
  }
  CloseInFile();
  if (values.size() != 0) return values;
  msg_Tracking()<<"Data_Reader: Parameter "<<parameter<<" not specified in "<<filename<<" !"<<std::endl;
  return values;
}

template <class Read_Type>
std::vector< std::vector<Read_Type> > 
Data_Reader::M_MatrixFromString(std::string parameter,std::string inputstring,MatrixTypeID tempmtype)
{
  if (tempmtype==MUnknown) tempmtype=MatrixType();
  if (tempmtype==MUnknown) tempmtype=MNormal;
  std::vector<Read_Type> value;
  std::vector< std::vector<Read_Type> > transposedvalues;
  std::string name;
  if (inputstring==noinputtag) {
    if (String()!=nullstring) {
      inputstring=String();
    }
    else {
      msg_Tracking()<<"Data_Reader: No input string specified ! No default available !"<<std::endl
			    <<"   Abort reading."<<std::endl;
      return transposedvalues;
    }
  }
  value=M_VectorFromString<Read_Type>(parameter,inputstring,VHorizontal);
  unsigned int mindim=NotDefined<unsigned int>();
  for(unsigned int i=0;value.size()!=0;++i) {
    transposedvalues.push_back(value);
    if (value.size()<mindim) mindim=value.size();
    bool foundseparator=false;
    size_t sep=std::string::npos;
    for (unsigned int i=0;i<Separator().size();++i) {
      if ((sep=inputstring.find(Separator()[i]))!=std::string::npos) foundseparator=true;
    }
    if(foundseparator) inputstring=inputstring.substr(sep+1);
    else inputstring=nullstring;
    value=M_VectorFromString<Read_Type>(nullstring,inputstring,VHorizontal);
  }
  if (transposedvalues.size()!=0) {
    if (tempmtype==MNormal) {
      std::vector< std::vector<Read_Type> > normalvalues;     
      for (unsigned int i=0;i<transposedvalues[0].size();++i) 
	normalvalues.push_back(std::vector<Read_Type>(transposedvalues.size()));
      for (unsigned int j=0;j<transposedvalues.size();++j) {
	for (unsigned int k=0;k<mindim;++k) {
	  normalvalues[k][j]=transposedvalues[j][k];
#ifdef DEBUG__Data_Reader
	  std::cout<<"   transposed values are ["<<j<<"]["<<k<<"] = "<<transposedvalues[j][k]<<std::endl;
#endif
	}
      }
      for (unsigned int i=0;i<transposedvalues.size();++i) transposedvalues[i].clear();
      transposedvalues.clear();
      return normalvalues;
    }
    else return transposedvalues;
  }
  return transposedvalues;
}

template <class Read_Type>
std::vector< std::vector<Read_Type> > 
Data_Reader::M_MatrixFromFile(std::string parameter,std::string filename,MatrixTypeID tempmtype)
{
  if (tempmtype==MUnknown) tempmtype=MatrixType();
  if (tempmtype==MUnknown) tempmtype=MNormal;
  std::vector< std::vector<Read_Type> > transposedvalues, temp;
  if (filename==noinputtag) {
    if (InputFile()!=nullstring) {
      filename=InputFile();
    }
    else {
      msg_Tracking()<<"Data_Reader: No input file specified ! No default available !"<<std::endl
			    <<"    Abort reading."<<std::endl;
      return transposedvalues;
    }
  }
  if (!OpenInFile()) {
    msg_Tracking()<<"Data_Reader: Error opening "<<filename<<" !"<<std::endl;
    return transposedvalues;
  }
  for (unsigned int i=0;i<FileContent().size();++i) {
    temp = M_MatrixFromString<Read_Type>(parameter,FileContent()[i],MTransposed);
    for (unsigned int j=0; j<temp.size(); ++j) {
      transposedvalues.push_back(temp[j]);
    }
  }
  CloseInFile();
  if (transposedvalues.size()!=0) {
    if (tempmtype==MNormal) {
      std::vector< std::vector<Read_Type> > normalvalues;     
      for (unsigned int i=0;i<transposedvalues[0].size();++i) 
	normalvalues.push_back(std::vector<Read_Type>(transposedvalues.size()));
      for (unsigned int j=0;j<transposedvalues.size();++j) {
	for (unsigned int k=0;k<transposedvalues[j].size();++k) {
#ifdef DEBUG__Data_Reader
	  std::cout<<"   transposed values are ["<<j<<"]["<<k<<"] = "<<transposedvalues[j][k]<<std::endl;
#endif
	  normalvalues[k][j]=transposedvalues[j][k]; 
	}
      }
      for (unsigned int i=0;i<transposedvalues.size();++i) transposedvalues[i].clear();
      transposedvalues.clear();
      return normalvalues;
    }
    else return transposedvalues;
  }
  msg_Tracking()<<"Data_Reader: Parameter "<<parameter<<" not specified in "<<filename<<" !"<<std::endl;
  return transposedvalues;
}

template <class Read_Type >
bool Data_Reader::ReadFromFile(Read_Type &result,std::string parameter,std::string filename) 
{ 
  if ((result=M_ReadFromFile<Read_Type>(parameter,filename))!=Default<Read_Type>()) return true; 
  else return false; 
}

template bool Data_Reader::ReadFromFile<int>(int &,std::string,std::string);
template bool Data_Reader::ReadFromFile<unsigned int>(unsigned int &,std::string,std::string);
template bool Data_Reader::ReadFromFile<long int>(long int &,std::string,std::string);
template bool Data_Reader::ReadFromFile<float>(float &,std::string,std::string);
template bool Data_Reader::ReadFromFile<double>(double &,std::string,std::string);
template bool Data_Reader::ReadFromFile<std::string>(std::string &,std::string,std::string);

template <class Read_Type >
bool Data_Reader::ReadFromString(Read_Type &result,std::string parameter,std::string inputstring) 
{ 
  if ((result=M_ReadFromString<Read_Type>(parameter,inputstring))!=Default<Read_Type>()) return true; 
  else return false; 
}

template bool Data_Reader::ReadFromString<int>(int &,std::string,std::string);
template bool Data_Reader::ReadFromString<unsigned int>(unsigned int &,std::string,std::string);
template bool Data_Reader::ReadFromString<long int>(long int &,std::string,std::string);
template bool Data_Reader::ReadFromString<float>(float &,std::string,std::string);
template bool Data_Reader::ReadFromString<double>(double &,std::string,std::string);
template bool Data_Reader::ReadFromString<std::string>(std::string &,std::string,std::string);

template <class Read_Type >
bool Data_Reader::VectorFromFile(std::vector<Read_Type> &result,std::string parameter,
				 std::string filename, VectorTypeID tempvtype) 
{ 
  if ((result=M_VectorFromFile<Read_Type>(parameter,filename,tempvtype)).size()!=0) return true; 
  else return false; 
}

template bool Data_Reader::VectorFromFile<int>(std::vector<int> &,std::string,std::string,VectorTypeID);
template bool Data_Reader::VectorFromFile<unsigned int>(std::vector<unsigned int> &,std::string,std::string,VectorTypeID);
template bool Data_Reader::VectorFromFile<long int>(std::vector<long int> &,std::string,std::string,VectorTypeID);
template bool Data_Reader::VectorFromFile<float>(std::vector<float> &,std::string,std::string,VectorTypeID);
template bool Data_Reader::VectorFromFile<double>(std::vector<double> &,std::string,std::string,VectorTypeID);
template bool Data_Reader::VectorFromFile<std::string>(std::vector<std::string> &,std::string,std::string,VectorTypeID);

template <class Read_Type>
bool Data_Reader::VectorFromString(std::vector<Read_Type> &result,std::string parameter,
				   std::string inputstring, VectorTypeID tempvtype) 
{ 
  if ((result=M_VectorFromString<Read_Type>(parameter,inputstring,tempvtype)).size()!=0) return true; 
  else return false; 
}

template bool Data_Reader::VectorFromString<int>(std::vector<int> &,std::string,std::string,VectorTypeID);
template bool Data_Reader::VectorFromString<unsigned int>(std::vector<unsigned int> &,std::string,std::string,VectorTypeID);
template bool Data_Reader::VectorFromString<long int>(std::vector<long int> &,std::string,std::string,VectorTypeID);
template bool Data_Reader::VectorFromString<float>(std::vector<float> &,std::string,std::string,VectorTypeID);
template bool Data_Reader::VectorFromString<double>(std::vector<double> &,std::string,std::string,VectorTypeID);
template bool Data_Reader::VectorFromString<std::string>(std::vector<std::string> &,std::string,std::string,VectorTypeID);

template <class Read_Type>
bool Data_Reader::MatrixFromFile(std::vector<std::vector<Read_Type> > &result,std::string parameter,
				  std::string filename, MatrixTypeID tempmtype) 
{ 
  if ((result=M_MatrixFromFile<Read_Type>(parameter,filename,tempmtype)).size()!=0) return true; 
  else return false; 
}

template bool Data_Reader::MatrixFromFile<int>(std::vector<std::vector<int> > &,std::string,std::string,MatrixTypeID); 
template bool Data_Reader::MatrixFromFile<unsigned int>(std::vector<std::vector<unsigned int> > &,std::string,std::string,MatrixTypeID);
template bool Data_Reader::MatrixFromFile<long int>(std::vector<std::vector<long int> > &,std::string,std::string,MatrixTypeID);
template bool Data_Reader::MatrixFromFile<float>(std::vector<std::vector<float> > &,std::string,std::string,MatrixTypeID);
template bool Data_Reader::MatrixFromFile<double>(std::vector<std::vector<double> > &,std::string,std::string,MatrixTypeID);
template bool Data_Reader::MatrixFromFile<std::string>(std::vector<std::vector<std::string> > &,std::string,std::string,MatrixTypeID);

template <class Read_Type>
bool Data_Reader::MatrixFromString(std::vector<std::vector<Read_Type> > &result,std::string parameter,
				   std::string inputstring, MatrixTypeID tempmtype) 
{ 
  if ((result=M_MatrixFromString<Read_Type>(parameter,inputstring,tempmtype)).size()!=0) return true; 
  else return false; 
}

template bool Data_Reader::MatrixFromString<int>(std::vector<std::vector<int> > &,std::string,std::string,MatrixTypeID);
template bool Data_Reader::MatrixFromString<unsigned int>(std::vector<std::vector<unsigned int> > &,std::string,std::string,MatrixTypeID);
template bool Data_Reader::MatrixFromString<long int>(std::vector<std::vector<long int> > &,std::string,std::string,MatrixTypeID);
template bool Data_Reader::MatrixFromString<float>(std::vector<std::vector<float> > &,std::string,std::string,MatrixTypeID);
template bool Data_Reader::MatrixFromString<double>(std::vector<std::vector<double> > &,std::string,std::string,MatrixTypeID);
template bool Data_Reader::MatrixFromString<std::string>(std::vector<std::vector<std::string> > &,std::string,std::string,MatrixTypeID);


