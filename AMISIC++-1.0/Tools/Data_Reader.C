#ifndef Data_Reader_C
#define Data_Reader_C

#include "Data_Reader.H"
#include "Message.H"
#include "Type.H"

#include "My_Filestream.H"
#include "My_Stringstream.H"

#ifdef DEBUG__Data_Reader
#include "My_IO_Stream.H"
#endif

namespace ATOOLS {

  Data_Reader::Data_Reader(): 
    Read_Write_Base() {}

  Data_Reader::Data_Reader(std::string _m_cut, std::string _m_seperator, std::string _m_comment):
    Read_Write_Base(_m_cut,_m_seperator,_m_comment) {}
  
  Data_Reader::Data_Reader(const char *_m_cut, const char *_m_seperator, const char *_m_comment):
    Read_Write_Base(_m_cut,_m_seperator,_m_comment) {}
  
  void Data_Reader::KillComments(std::string& buffer)
  {
#ifdef DEBUG__Data_Reader
    std::cout<<"Data_Reader::KillComments("<<buffer<<")"<<std::endl;
    std::cout<<"   comment tags are ";
    for (unsigned int i=0;i<M_Comment().size();std::cout<<"'"<<M_Comment()[i++]<<"' ");
    std::cout<<std::endl;
    std::cout<<"   ignored tags are ";
    for (unsigned int i=0;i<M_Ignore().size();std::cout<<"'"<<M_Ignore()[i++]<<"' ");
    std::cout<<std::endl;
#endif
    int pos;
    for (unsigned int i=0; i<M_Comment().size(); ++i) {
      if ((pos = buffer.find(M_Comment()[i]))!=-1) {
	buffer=buffer.substr(0,pos);
      }
    }
    KillBlanks(buffer);
    for (unsigned int i=0; i<M_Ignore().size(); ++i) {
      while ((pos = buffer.find(M_Ignore()[i]))!=-1) {
	buffer=buffer.substr(0,pos)+buffer.substr(pos+M_Ignore()[i].length());
	KillBlanks(buffer);
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
    if (buffer==ATOOLS::nullstring) return buffer;
    bool hit;
    do { 
      hit = false;
      for (unsigned int i=0; i<M_Blank().size(); ++i) hit = hit || (int(buffer[0])==M_Blank()[i]);
      if (hit) {
	buffer=buffer.substr(1);
      }
    } while (hit);
    do { 
      hit = false;
      for (unsigned int i=0; i<M_Blank().size(); ++i) hit = hit || (int(buffer[buffer.length()-1])==M_Blank()[i]);
      if (hit) buffer=buffer.substr(0,buffer.length()-1);
    } while (hit);
#ifdef DEBUG__Data_Reader
    std::cout<<"   returning '"<<buffer<<"'"<<std::endl;
#endif
    return buffer;
  }
  
  template <class Read_Type>
  Read_Type Data_Reader::M_ReadFromString(std::string parameter,std::string inputstring)
  {
    int pos;
    Read_Type value;
#ifdef DEBUG__Data_Reader
    std::cout<<"Data_Reader::M_ReadFromString("<<parameter<<","<<inputstring<<")"<<std::endl;
#endif
    if (inputstring==noinputtag) {
      if (M_String()!=nullstring) {
	inputstring=M_String();
      }
      else {
	ATOOLS::msg.Error()<<"Data_Reader: No input string specified ! No default available !"<<std::endl
			   <<"            Abort reading."<<std::endl;
	return Default<Read_Type>();
      }
    }
    KillComments(inputstring);
    if (parameter==ATOOLS::nullstring) { 
      parameter=std::string("DUMMY_PARAMETER");
      inputstring=parameter+inputstring;
    }
    if(((pos=inputstring.find(parameter))!=-1)&&
       ((inputstring=inputstring.substr(pos+parameter.length())).length()>0)) {
      std::stringstream converter;
      converter<<inputstring;
      converter>>value;
#ifdef DEBUG__Data_Reader
      Type typeinfo;
      std::cout<<"   returning '"<<value<<"'"
	       <<" ( type = "<<typeinfo.GetType(value)<<" )"<<std::endl;
#endif
      return value;
    }
#ifdef DEBUG__Data_Reader
    Type typeinfo;
    std::cout<<"   returning default value '"<<value<<"'"
	     <<" ( type "<<typeinfo.GetType(value)<<" )"<<std::endl;
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
      if (M_FileName()!=nullstring) {
	filename=M_FileName();
      }
      else {
	ATOOLS::msg.Error()<<"Data_Reader: No input file specified ! No default available !"<<std::endl
			   <<"            Abort reading."<<std::endl;
	return value;
      }
    }
    if(!OpenFile(filename,std::ios_base::in)) {
      ATOOLS::msg.Out()<<"Data_Reader: Error opening "<<filename<<" !"<<std::endl;
      return value;
    }
    while (*M_File()) {
      getline(*M_File(),buffer);
      if((temp=M_ReadFromString<Read_Type>(parameter,buffer))!=Default<Read_Type>()) value=temp;
    }
    CloseFile();
    if (value==Default<Read_Type>()) {
      ATOOLS::msg.Out()<<"Data_Reader: Parameter "<<parameter<<" not specified in "<<filename<<" !"<<std::endl;
    }
    return value;
  }
  
  template <class Read_Type>
  std::vector<Read_Type> 
  Data_Reader::M_VectorFromString(std::string parameter, std::string inputstring,VectorType tempvtype)
  {
    if (tempvtype==VUnknown) tempvtype=M_VectorType();
    if (tempvtype==VUnknown) tempvtype=VVertical;
    std::string value;
    Type typeinfo;
    std::vector<Read_Type> values;
    int i = 0;
#ifdef DEBUG__Data_Reader
    std::cout<<"Data_Reader::M_VectorFromString("<<parameter<<","
	     <<inputstring<<","<<tempvtype<<")"<<std::endl;
#endif
    if (inputstring==noinputtag) {
      if (M_String()!=nullstring) {
	inputstring=M_String();
      }
      else {
	ATOOLS::msg.Error()<<"Data_Reader: No input string specified ! No default available !"<<std::endl
			   <<"            Abort reading."<<std::endl;
	return values;
      }
    }
    value = M_ReadFromString<std::string>(parameter,inputstring);
    for(;value != Default<std::string>();++i) {
#ifdef DEBUG__Data_Reader
      std::cout<<"   newline tags are ";
      for (unsigned int j=0;j<M_Seperator().size();std::cout<<"'"<<M_Seperator()[j++]<<"' ");
      std::cout<<std::endl;
#endif
      bool stop=false;
      for (unsigned int j=0;j<M_Seperator().size();++j) if (value==M_Seperator()[j]) stop=true;
      if (stop) break;
      values.push_back(M_ReadFromString<Read_Type>(nullstring,value));
      if (tempvtype==VVertical) break;
      inputstring=inputstring.substr(inputstring.find(value)+value.length());
      if ((typeinfo.GetType(value)!=Type::TChar) &&
	  (typeinfo.GetType(value)!=Type::TString)) {
	if ((inputstring.find(std::string("."))==0) ||
	    (inputstring.find(std::string("0"))==0)) {
	  inputstring = inputstring.substr(1);
	  for (;inputstring.find(std::string("0"))==0;
	       inputstring = inputstring.substr(1));
	}  
      }
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
  Data_Reader::M_VectorFromFile(std::string parameter, std::string filename,VectorType tempvtype)
  {
    if (tempvtype==VUnknown) tempvtype=M_VectorType();
    if (tempvtype==VUnknown) tempvtype=VVertical;
    std::vector<Read_Type> values, temp;
    std::string buffer;
    if (filename==noinputtag) {
      if (M_FileName()!=nullstring) {
	filename=M_FileName();
      }
      else {
	ATOOLS::msg.Error()<<"Data_Reader: No input file specified ! No default available !"<<std::endl
			   <<"            Abort reading."<<std::endl;
	return values;
      }
    }
    if (!OpenFile(filename,std::ios_base::in)) {
      ATOOLS::msg.Out()<<"Data_Reader: Error opening "<<filename<<" !"<<std::endl;
      return values;
    }
    while (*M_File()) {
      getline(*M_File(),buffer);
      temp = M_VectorFromString<Read_Type>(parameter,buffer,tempvtype);
      switch (tempvtype) {
      case VHorizontal:
	if (temp.size()!=0) values=temp;
	break;
      case VUnknown:
      case VVertical:
	for (unsigned int i=0; i<temp.size(); ++i) values.push_back(temp[i]);
	break;
      }
    }
    CloseFile();
    if (values.size() != 0) return values;
    ATOOLS::msg.Out()<<"Data_Reader: Parameter "<<parameter<<" not specified in "<<filename<<" !"<<std::endl;
    return values;
  }

  template <class Read_Type>
  std::vector< std::vector<Read_Type> > 
  Data_Reader::M_ArrayFromString(std::string parameter,std::string inputstring,MatrixType tempmtype)
  {
    if (tempmtype==MUnknown) tempmtype=M_MatrixType();
    if (tempmtype==MUnknown) tempmtype=MNormal;
    std::vector<Read_Type> value;
    std::vector< std::vector<Read_Type> > transposedvalues;
    std::string name;
    int sep, i = 0;
    if (inputstring==noinputtag) {
      if (M_String()!=nullstring) {
	inputstring=M_String();
      }
      else {
	ATOOLS::msg.Error()<<"Data_Reader: No input string specified ! No default available !"<<std::endl
			   <<"            Abort reading."<<std::endl;
	return transposedvalues;
      }
    }
    value=M_VectorFromString<Read_Type>(parameter,inputstring,VHorizontal);
    for(;value.size()!=0;++i) {
      transposedvalues.push_back(value);
      bool foundseperator=false;
      for (unsigned int i=0;i<M_Seperator().size();++i) {
	if ((sep=inputstring.find(M_Seperator()[i]))!=-1) foundseperator=true;
      }
      if(foundseperator) inputstring=inputstring.substr(sep+1);
      else inputstring=nullstring;
      value=M_VectorFromString<Read_Type>(nullstring,inputstring,VHorizontal);
    }
    if (transposedvalues.size()!=0) {
      if (tempmtype==MNormal) {
	std::vector< std::vector<Read_Type> > normalvalues;     
	for (unsigned int i=0;i<transposedvalues[0].size();++i) 
	  normalvalues.push_back(std::vector<Read_Type>(transposedvalues.size()));
	for (unsigned int j=0;j<transposedvalues.size();++j) {
	  for (unsigned int k=0;k<transposedvalues.size();++k) {
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
  Data_Reader::M_ArrayFromFile(std::string parameter,std::string filename,MatrixType tempmtype)
  {
    if (tempmtype==MUnknown) tempmtype=M_MatrixType();
    if (tempmtype==MUnknown) tempmtype=MNormal;
    std::vector< std::vector<Read_Type> > transposedvalues, temp;
    std::string buffer;
    if (filename==noinputtag) {
      if (M_FileName()!=nullstring) {
	filename=M_FileName();
      }
      else {
	ATOOLS::msg.Error()<<"Data_Reader: No input file specified ! No default available !"<<std::endl
			  <<"            Abort reading."<<std::endl;
	return transposedvalues;
      }
    }
    if (!OpenFile(filename,std::ios_base::in)) {
      ATOOLS::msg.Out()<<"Data_Reader: Error opening "<<filename<<" !"<<std::endl;
      return transposedvalues;
    }
    while (*M_File()) {
      getline(*M_File(),buffer);
      temp = M_ArrayFromString<Read_Type>(parameter,buffer,MTransposed);
      for (unsigned int i=0; i<temp.size(); ++i) {
	transposedvalues.push_back(temp[i]);
      }
    }
    CloseFile();
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
    ATOOLS::msg.Out()<<"Data_Reader: Parameter "<<parameter<<" not specified in "<<filename<<" !"<<std::endl;
    return transposedvalues;
  }
  
  bool Data_Reader::ReadFromFile(int &result,std::string parameter,std::string filename) 
  { if ((result=M_ReadFromFile<int>(parameter,filename))!=Default<int>()) return true; else return false; }
  bool Data_Reader::ReadFromFile(long int &result,std::string parameter,std::string filename) 
  { if ((result=M_ReadFromFile<long int>(parameter,filename))!=Default<long int>()) return true; else return false; }
  bool Data_Reader::ReadFromFile(float &result,std::string parameter,std::string filename) 
  { if ((result=M_ReadFromFile<float>(parameter,filename))!=Default<float>()) return true; else return false; }
  bool Data_Reader::ReadFromFile(double &result,std::string parameter,std::string filename) 
  { if ((result=M_ReadFromFile<double>(parameter,filename))!=Default<double>()) return true; else return false; }
  bool Data_Reader::ReadFromFile(std::string &result,std::string parameter,std::string filename)
  { if ((result=M_ReadFromFile<std::string>(parameter,filename))!=Default<std::string>()) return true; else return false; }
  
  bool Data_Reader::ReadFromString(int &result,std::string parameter,std::string filename) 
  { if ((result=M_ReadFromString<int>(parameter,filename))!=Default<int>()) return true; else return false; }
  bool Data_Reader::ReadFromString(long int &result,std::string parameter,std::string filename) 
  { if ((result=M_ReadFromString<long int>(parameter,filename))!=Default<long int>()) return true; else return false; }
  bool Data_Reader::ReadFromString(float &result,std::string parameter,std::string filename) 
  { if ((result=M_ReadFromString<float>(parameter,filename))!=Default<float>()) return true; else return false; }
  bool Data_Reader::ReadFromString(double &result,std::string parameter,std::string filename) 
  { if ((result=M_ReadFromString<double>(parameter,filename))!=Default<double>()) return true; else return false; }
  bool Data_Reader::ReadFromString(std::string &result,std::string parameter,std::string filename) 
  { if ((result=M_ReadFromString<std::string>(parameter,filename))!=Default<std::string>()) return true; else return false; }
  
  bool Data_Reader::VectorFromFile(std::vector<int> &result,std::string parameter,
				   std::string filename, VectorType tempvtype) 
  { if ((result=M_VectorFromFile<int>(parameter,filename,tempvtype)).size()!=0) return true; else return false; }
  bool Data_Reader::VectorFromFile(std::vector<long int> &result,std::string parameter,
				   std::string filename, VectorType tempvtype) 
  { if ((result=M_VectorFromFile<long int>(parameter,filename,tempvtype)).size()!=0) return true; else return false; }
  bool Data_Reader::VectorFromFile(std::vector<float> &result,std::string parameter,
				    std::string filename, VectorType tempvtype) 
  { if ((result=M_VectorFromFile<float>(parameter,filename,tempvtype)).size()!=0) return true; else return false; }
  bool Data_Reader::VectorFromFile(std::vector<double> &result,std::string parameter,
				    std::string filename, VectorType tempvtype) 
  { if ((result=M_VectorFromFile<double>(parameter,filename,tempvtype)).size()!=0) return true; else return false; }
  bool Data_Reader::VectorFromFile(std::vector<std::string> &result,std::string parameter,
				    std::string filename, VectorType tempvtype) 
  { if ((result=M_VectorFromFile<std::string>(parameter,filename,tempvtype)).size()!=0) return true; else return false; }
  
  bool Data_Reader::VectorFromString(std::vector<int> &result,std::string parameter,
				      std::string filename, VectorType tempvtype) 
  { if ((result=M_VectorFromString<int>(parameter,filename,tempvtype)).size()!=0) return true; else return false; }
  bool Data_Reader::VectorFromString(std::vector<long int> &result,std::string parameter,
				      std::string filename, VectorType tempvtype) 
  { if ((result=M_VectorFromString<long int>(parameter,filename,tempvtype)).size()!=0) return true; else return false; }
  bool Data_Reader::VectorFromString(std::vector<float> &result,std::string parameter,
				      std::string filename, VectorType tempvtype) 
  { if ((result=M_VectorFromString<float>(parameter,filename,tempvtype)).size()!=0) return true; else return false; }
  bool Data_Reader::VectorFromString(std::vector<double> &result,std::string parameter,
				      std::string filename, VectorType tempvtype) 
  { if ((result=M_VectorFromString<double>(parameter,filename,tempvtype)).size()!=0) return true; else return false; }
  bool Data_Reader::VectorFromString(std::vector<std::string> &result,std::string parameter,
				      std::string filename, VectorType tempvtype) 
  { if ((result=M_VectorFromString<std::string>(parameter,filename,tempvtype)).size()!=0) return true; else return false; }
  
  bool Data_Reader::ArrayFromFile(std::vector<std::vector<int> > &result,std::string parameter,
				  std::string filename, MatrixType tempmtype) 
  { if ((result=M_ArrayFromFile<int>(parameter,filename,tempmtype)).size()!=0) return true; else return false; }
  bool Data_Reader::ArrayFromFile(std::vector<std::vector<long int> > &result,std::string parameter,
				  std::string filename, MatrixType tempmtype) 
  { if ((result=M_ArrayFromFile<long int>(parameter,filename,tempmtype)).size()!=0) return true; else return false; }
  bool Data_Reader::ArrayFromFile(std::vector<std::vector<float> > &result,std::string parameter,
				  std::string filename, MatrixType tempmtype) 
  { if ((result=M_ArrayFromFile<float>(parameter,filename,tempmtype)).size()!=0) return true; else return false; }
  bool Data_Reader::ArrayFromFile(std::vector<std::vector<double> > &result,std::string parameter,
				  std::string filename, MatrixType tempmtype) 
  { if ((result=M_ArrayFromFile<double>(parameter,filename,tempmtype)).size()!=0) return true; else return false; }
  bool Data_Reader::ArrayFromFile(std::vector<std::vector<std::string> > &result,std::string parameter,
				  std::string filename, MatrixType tempmtype) 
  { if ((result=M_ArrayFromFile<std::string>(parameter,filename,tempmtype)).size()!=0) return true; else return false; }
  
  bool Data_Reader::ArrayFromString(std::vector<std::vector<int> > &result,std::string parameter,
				    std::string filename, MatrixType tempmtype) 
  { if ((result=M_ArrayFromString<int>(parameter,filename,tempmtype)).size()!=0) return true; else return false; }
  bool Data_Reader::ArrayFromString(std::vector<std::vector<long int> > &result,std::string parameter,
				    std::string filename, MatrixType tempmtype) 
  { if ((result=M_ArrayFromString<long int>(parameter,filename,tempmtype)).size()!=0) return true; else return false; }
  bool Data_Reader::ArrayFromString(std::vector<std::vector<float> > &result,std::string parameter,
				    std::string filename, MatrixType tempmtype) 
  { if ((result=M_ArrayFromString<float>(parameter,filename,tempmtype)).size()!=0) return true; else return false; }
  bool Data_Reader::ArrayFromString(std::vector<std::vector<double> > &result,std::string parameter,
				    std::string filename, MatrixType tempmtype) 
  { if ((result=M_ArrayFromString<double>(parameter,filename,tempmtype)).size()!=0) return true; else return false; }
  bool Data_Reader::ArrayFromString(std::vector<std::vector<std::string> > &result,std::string parameter,
				    std::string filename, MatrixType tempmtype) 
  { if ((result=M_ArrayFromString<std::string>(parameter,filename,tempmtype)).size()!=0) return true; else return false; }

} // end of namespace ATOOLS

#endif

