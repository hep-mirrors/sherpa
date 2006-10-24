#include "Data_Reader.H"

#include "Message.H"
#include "MyStrStream.H"
#include <typeinfo>
#include <ctype.h>

using namespace ATOOLS;

Data_Reader::Data_Reader(): 
  Read_Write_Base(1,0)
{
  SetInFileMode(fom::permanent);
}

Data_Reader::Data_Reader(const std::string cut,
			 const std::string wordsep,
			 const std::string linesep,
			 const std::string comment):
  Read_Write_Base(1,0,cut,wordsep,linesep,comment)
{
  SetInFileMode(fom::permanent);
}

void Data_Reader::KillComments(std::string& buffer)
{
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
}

std::string Data_Reader::KillBlanks(std::string& buffer)
{
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
  return buffer;
}

template <class Read_Type>
Read_Type Data_Reader::M_ReadFromString(std::string parameter,
					std::string &inputstring,
					int &septype)
{
  size_t pos;
  Read_Type value;
  if (inputstring==noinputtag) {
    if (String()!=nullstring) {
      inputstring=String();
    }
    else {
      msg_Debugging()<<"Data_Reader: No input string specified !\n";
      return Default<Read_Type>();
    }
  }
  KillComments(inputstring);
  if (parameter==nullstring) { 
    parameter=std::string("DUMMY_PARAMETER");
    inputstring=parameter+" "+inputstring;
  }
  size_t length=0;
  if(((pos=Find(inputstring,parameter,length))!=std::string::npos)&&
     ((inputstring=inputstring.substr(pos+length)).length()>0)) {
    std::string cur=inputstring;
    cur=KillBlanks(cur);
    size_t wi(0), wpos(0);
    for (;wi<WordSeparator().size();++wi) {
      while (true) {
	wpos=cur.find(WordSeparator()[wi],wpos==0?0:wpos+1);
	if (wpos==std::string::npos ||
	    wpos<1 || cur[wpos-1]!=Escape()) break;
      }
      if (wpos!=std::string::npos) break;
    }
    size_t li(0), lpos(0);
    for (;li<LineSeparator().size();++li) {
      while (true) {
	lpos=cur.find(LineSeparator()[li],lpos==0?0:lpos+1);
	if (lpos==std::string::npos || 
	    lpos<1 || cur[lpos-1]!=Escape()) break;
      }
      if (lpos!=std::string::npos) break;
    }
    septype=0;
    if (lpos==std::string::npos && wpos==std::string::npos) {
      inputstring=nullstring;
    }
    else if (wpos<lpos) {
      inputstring=cur.substr(wpos+WordSeparator()[wi].length());
      cur=cur.substr(0,wpos);
      septype=1;
    }
    else {
      inputstring=cur.substr(lpos+LineSeparator()[li].length());
      cur=cur.substr(0,lpos);
      septype=2;
    }
    cur=ReplaceTags(cur);
    if (typeid(value)==typeid(int) || typeid(value)==typeid(unsigned int) ||
	typeid(value)==typeid(float) ||	typeid(value)==typeid(double)) {
      if (!m_allownans) {
 	if ((pos=cur.find("nan"))!=std::string::npos) 
	  cur.replace(pos,3,"1");
 	else if ((pos=cur.find("inf"))!=std::string::npos) 
 	  cur.replace(pos,3,"1");
      }
      if (Interprete()) cur=Interpreter()->Interprete(StripEscapes(cur));
    }
    value=ATOOLS::ToType<Read_Type>(cur);
    return value;
  }
  return Default<Read_Type>();
}

template <class Read_Type>
Read_Type Data_Reader::M_ReadFromFile(std::string parameter,
				      std::string filename)
{
  Read_Type value, temp;
  std::string buffer;
  value=Default<Read_Type>();
  if (filename==noinputtag) {
    if (InputFile()!=nullstring) {
      filename=InputFile();
    }
    else {
      msg_Debugging()<<"Data_Reader: No input file specified. "
		     <<"No default available."<<std::endl
		     <<"   Abort reading."<<std::endl;
      return value;
    }
  }
  if(!OpenInFile()) {
    msg.Tracking()<<"Data_Reader: File '"<<filename
		  <<"' not found."<<std::endl;
    return value;
  }
  for (unsigned int i=0;i<FileContent().size();++i) {
    buffer=FileContent()[i];
    int septype(0);
    if((temp=M_ReadFromString<Read_Type>(parameter,buffer,septype))!=
       Default<Read_Type>()) value=temp;
  }
  CloseInFile();
  if (value==Default<Read_Type>()) {
    msg_Debugging()<<"Data_Reader: Parameter '"<<parameter
		   <<"' not specified in '"<<filename<<"' !"<<std::endl;
  }
  return value;
}

template <class Read_Type> std::vector<Read_Type> 
Data_Reader::M_VectorFromString(std::string parameter, 
				std::string inputstring,vtc::code tempvtype)
{
  if (tempvtype==vtc::unknown) tempvtype=VectorType();
  if (tempvtype==vtc::unknown) tempvtype=vtc::vertical;
  std::string value;
  std::vector<Read_Type> values;
  inputstring = ReplaceTags(inputstring);
  if (inputstring==noinputtag) {
    if (String()!=nullstring) {
      inputstring=String();
    }
    else {
      msg_Debugging()<<"Data_Reader: No input string specified. "
		     <<"No default available !"<<std::endl
		     <<"   Abort reading."<<std::endl;
      return values;
    }
  }
  int septype(0);
  value = M_ReadFromString<std::string>(parameter,inputstring,septype);
  while (value != Default<std::string>()) {
    int dseptype(0);
    values.push_back(M_ReadFromString<Read_Type>(nullstring,value,dseptype));
    if (septype==2 || tempvtype==vtc::vertical) break;
    value=M_ReadFromString<std::string>(nullstring,inputstring,septype);
  }
  if (values.size()!=0) return values;
  return values;
}

template <class Read_Type> std::vector<Read_Type> 
Data_Reader::M_VectorFromFile(std::string parameter, 
			      std::string filename,vtc::code tempvtype)
{
  if (tempvtype==vtc::unknown) tempvtype=VectorType();
  if (tempvtype==vtc::unknown) tempvtype=vtc::vertical;
  std::vector<Read_Type> values, temp;
  if (filename==noinputtag) {
    if (InputFile()!=nullstring) {
      filename=InputFile();
    }
    else {
      msg_Debugging()<<"Data_Reader: No input file specified !\n";
      return values;
    }
  }
  if (!OpenInFile()) {
    msg_Tracking()<<"Data_Reader: Error opening "<<filename
		  <<" !"<<std::endl;
    return values;
  }
  for (unsigned int i=0;i<FileContent().size();++i) {
    temp = M_VectorFromString<Read_Type>
      (parameter,FileContent()[i],tempvtype);
    switch (tempvtype) {
    case vtc::horizontal:
      if (temp.size()!=0) values=temp;
      break;
    case vtc::unknown:
    case vtc::vertical:
      for (unsigned int j=0;j<temp.size();++j) values.push_back(temp[j]);
      break;
    }
  }
  CloseInFile();
  if (values.size() != 0) return values;
  msg_Debugging()<<"Data_Reader: Parameter "
		 <<parameter<<" not specified in "
		 <<filename<<" !"<<std::endl;
  return values;
}

template <class Read_Type> std::vector< std::vector<Read_Type> > 
Data_Reader::M_MatrixFromString(std::string parameter,
				std::string inputstring,mtc::code tempmtype)
{
  if (tempmtype==mtc::unknown) tempmtype=MatrixType();
  if (tempmtype==mtc::unknown) tempmtype=mtc::normal;
  std::vector<Read_Type> value;
  std::vector< std::vector<Read_Type> > transposedvalues;
  std::string name;
  if (inputstring==noinputtag) {
    if (String()!=nullstring) {
      inputstring=String();
    }
    else {
      msg_Debugging()<<"Data_Reader: No input string specified !\n";
      return transposedvalues;
    }
  }
  value=M_VectorFromString<Read_Type>
    (parameter,inputstring,vtc::horizontal);
  unsigned int mindim(0);
  for(unsigned int i=0;value.size()!=0;++i) {
    transposedvalues.push_back(value);
    if (value.size()<mindim) mindim=value.size();
    size_t i(0), sep(0);
    for (;i<LineSeparator().size();++i) {
      while (true) {
	sep=inputstring.find(LineSeparator()[i],sep==0?0:sep+1);
	if (sep==std::string::npos || 
	    sep<1 || inputstring[sep-1]!=Escape()) break;
      }
      if (sep!=std::string::npos) {
	inputstring=inputstring.substr(sep+LineSeparator()[i].length());
	break;
      }
    }
    if(i==LineSeparator().size()) inputstring=nullstring;
    value=M_VectorFromString<Read_Type>
      (nullstring,inputstring,vtc::horizontal);
  }
  if (transposedvalues.size()!=0) {
    if (tempmtype==mtc::normal) {
      std::vector< std::vector<Read_Type> > normalvalues;     
      for (unsigned int i=0;i<transposedvalues[0].size();++i) 
	normalvalues.push_back(std::vector<Read_Type>
			       (transposedvalues.size()));
      for (unsigned int j=0;j<transposedvalues.size();++j) {
	for (unsigned int k=0;k<mindim;++k) {
	  normalvalues[k][j]=transposedvalues[j][k];
	}
      }
      for (unsigned int i=0;i<transposedvalues.size();++i) 
	transposedvalues[i].clear();
      transposedvalues.clear();
      return normalvalues;
    }
    else return transposedvalues;
  }
  return transposedvalues;
}

template <class Read_Type> std::vector< std::vector<Read_Type> > 
Data_Reader::M_MatrixFromFile(std::string parameter,
			      std::string filename,mtc::code tempmtype)
{
  if (tempmtype==mtc::unknown) tempmtype=MatrixType();
  if (tempmtype==mtc::unknown) tempmtype=mtc::normal;
  std::vector< std::vector<Read_Type> > transposedvalues, temp;
  if (filename==noinputtag) {
    if (InputFile()!=nullstring) {
      filename=InputFile();
    }
    else {
      msg_Debugging()<<"Data_Reader: No input file specified !\n";
      return transposedvalues;
    }
  }
  if (!OpenInFile()) {
    msg_Tracking()<<"Data_Reader: Error opening "<<filename
		  <<" !"<<std::endl;
    return transposedvalues;
  }
  for (unsigned int i=0;i<FileContent().size();++i) {
    temp = M_MatrixFromString<Read_Type>
      (parameter,FileContent()[i],mtc::transposed);
    for (unsigned int j=0; j<temp.size(); ++j) {
      transposedvalues.push_back(temp[j]);
    }
  }
  CloseInFile();
  if (transposedvalues.size()!=0) {
    if (tempmtype==mtc::normal) {
      std::vector< std::vector<Read_Type> > normalvalues;     
      for (unsigned int i=0;i<transposedvalues[0].size();++i) 
	normalvalues.push_back(std::vector<Read_Type>
			       (transposedvalues.size()));
      for (unsigned int j=0;j<transposedvalues.size();++j) {
	for (unsigned int k=0;k<transposedvalues[j].size();++k) {
	  normalvalues[k][j]=transposedvalues[j][k]; 
	}
      }
      for (unsigned int i=0;i<transposedvalues.size();++i) 
	transposedvalues[i].clear();
      transposedvalues.clear();
      return normalvalues;
    }
    else return transposedvalues;
  }
  msg_Debugging()<<"Data_Reader: Parameter "
		 <<parameter<<" not specified in "<<filename
		 <<" !"<<std::endl;
  return transposedvalues;
}

template <class Read_Type > bool Data_Reader::
ReadFromFile(Read_Type &result,
	     std::string parameter,std::string filename) 
{ 
  if ((result=M_ReadFromFile<Read_Type>(parameter,filename))!=
      Default<Read_Type>()) return true; 
  else return false; 
}

template bool Data_Reader::ReadFromFile<int>
(int &,std::string,std::string);
template bool Data_Reader::ReadFromFile<unsigned int>
(unsigned int &,std::string,std::string);
template bool Data_Reader::ReadFromFile<long int>
(long int &,std::string,std::string);
template bool Data_Reader::ReadFromFile<float>
(float &,std::string,std::string);
template bool Data_Reader::ReadFromFile<double>
(double &,std::string,std::string);
template bool Data_Reader::ReadFromFile<std::string>
(std::string &,std::string,std::string);

template <class Read_Type > bool Data_Reader::
ReadFromString(Read_Type &result,
	       std::string parameter,std::string inputstring) 
{ 
  int septype(0);
  if ((result=M_ReadFromString<Read_Type>(parameter,inputstring,septype))!=
      Default<Read_Type>()) return true; 
  else return false; 
}

template bool Data_Reader::ReadFromString<int>
(int &,std::string,std::string);
template bool Data_Reader::ReadFromString<unsigned int>
(unsigned int &,std::string,std::string);
template bool Data_Reader::ReadFromString<long int>
(long int &,std::string,std::string);
template bool Data_Reader::ReadFromString<float>
(float &,std::string,std::string);
template bool Data_Reader::ReadFromString<double>
(double &,std::string,std::string);
template bool Data_Reader::ReadFromString<std::string>
(std::string &,std::string,std::string);

template <class Read_Type > bool Data_Reader::
VectorFromFile(std::vector<Read_Type> &result,std::string parameter,
	       std::string filename, vtc::code tempvtype) 
{ 
  if ((result=M_VectorFromFile<Read_Type>
       (parameter,filename,tempvtype)).size()!=0) return true; 
  else return false; 
}

template bool Data_Reader::VectorFromFile<int>
(std::vector<int> &,std::string,std::string,vtc::code);
template bool Data_Reader::VectorFromFile<unsigned int>
(std::vector<unsigned int> &,std::string,std::string,vtc::code);
template bool Data_Reader::VectorFromFile<long int>
(std::vector<long int> &,std::string,std::string,vtc::code);
template bool Data_Reader::VectorFromFile<float>
(std::vector<float> &,std::string,std::string,vtc::code);
template bool Data_Reader::VectorFromFile<double>
(std::vector<double> &,std::string,std::string,vtc::code);
template bool Data_Reader::VectorFromFile<std::string>
(std::vector<std::string> &,std::string,std::string,vtc::code);

template <class Read_Type> bool Data_Reader::
VectorFromString(std::vector<Read_Type> &result,std::string parameter,
		 std::string inputstring, vtc::code tempvtype) 
{ 
  if ((result=M_VectorFromString<Read_Type>
       (parameter,inputstring,tempvtype)).size()!=0) return true; 
  else return false; 
}

template bool Data_Reader::VectorFromString<int>
(std::vector<int> &,std::string,std::string,vtc::code);
template bool Data_Reader::VectorFromString<unsigned int>
(std::vector<unsigned int> &,std::string,std::string,vtc::code);
template bool Data_Reader::VectorFromString<long int>
(std::vector<long int> &,std::string,std::string,vtc::code);
template bool Data_Reader::VectorFromString<float>
(std::vector<float> &,std::string,std::string,vtc::code);
template bool Data_Reader::VectorFromString<double>
(std::vector<double> &,std::string,std::string,vtc::code);
template bool Data_Reader::VectorFromString<std::string>
(std::vector<std::string> &,std::string,std::string,vtc::code);

template <class Read_Type> bool Data_Reader::
MatrixFromFile(std::vector<std::vector<Read_Type> > &result,
	       std::string parameter,
	       std::string filename, mtc::code tempmtype) 
{ 
  if ((result=M_MatrixFromFile<Read_Type>
       (parameter,filename,tempmtype)).size()!=0) return true; 
  else return false; 
}

template bool Data_Reader::MatrixFromFile<int>
(std::vector<std::vector<int> > &,std::string,std::string,mtc::code); 
template bool Data_Reader::MatrixFromFile<unsigned int>
(std::vector<std::vector<unsigned int> > &,std::string,
 std::string,mtc::code);
template bool Data_Reader::MatrixFromFile<long int>
(std::vector<std::vector<long int> > &,std::string,std::string,mtc::code);
template bool Data_Reader::MatrixFromFile<float>
(std::vector<std::vector<float> > &,std::string,std::string,mtc::code);
template bool Data_Reader::MatrixFromFile<double>
(std::vector<std::vector<double> > &,std::string,std::string,mtc::code);
template bool Data_Reader::MatrixFromFile<std::string>
(std::vector<std::vector<std::string> > &,std::string,
 std::string,mtc::code);

template <class Read_Type> bool Data_Reader::
MatrixFromString(std::vector<std::vector<Read_Type> > &result,
		 std::string parameter,
		 std::string inputstring, mtc::code tempmtype) 
{ 
  if ((result=M_MatrixFromString<Read_Type>
       (parameter,inputstring,tempmtype)).size()!=0) return true; 
  else return false; 
}

template bool Data_Reader::MatrixFromString<int>
(std::vector<std::vector<int> > &,std::string,std::string,mtc::code);
template bool Data_Reader::MatrixFromString<unsigned int>
(std::vector<std::vector<unsigned int> > &,std::string,
 std::string,mtc::code);
template bool Data_Reader::MatrixFromString<long int>
(std::vector<std::vector<long int> > &,std::string,std::string,mtc::code);
template bool Data_Reader::MatrixFromString<float>
(std::vector<std::vector<float> > &,std::string,std::string,mtc::code);
template bool Data_Reader::MatrixFromString<double>
(std::vector<std::vector<double> > &,std::string,std::string,mtc::code);
template bool Data_Reader::MatrixFromString<std::string>
(std::vector<std::vector<std::string> > &,std::string,
 std::string,mtc::code);

