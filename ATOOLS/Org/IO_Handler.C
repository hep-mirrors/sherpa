#include <algorithm>
#include "ATOOLS/Org/IO_Handler.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Math/MyComplex.H"
#include "ATOOLS/Org/Message.H"

using namespace ATOOLS;
using namespace std;

IO_Handler::IO_Handler() 
{
  m_seps.push_back(';');
  m_coms.push_back('#');
  m_filename=std::string(""); 
}

IO_Handler::~IO_Handler() { 
  if (!(m_filename==std::string(""))) {
    m_file.close();
  }
}
    
// set output filename
int IO_Handler::SetFileName(std::string _name) {
  if (!(m_filename==std::string(""))) {
    m_file.close();
  }
  m_filename=_name;
  m_file.open(m_filename.c_str(),ios::out);

  if (!(m_file.good())) {
    msg_Out()<<" ERROR: opening "<<m_filename<<endl;
    return 0;
  }
  m_file.precision(15);
  return 1;
}
    
// set input filename
int IO_Handler::SetFileNameRO(string _name) {
  if (!(m_filename==std::string(""))) {
    m_file.close();
  }
  m_filename=_name;
  m_file.open(m_filename.c_str(),ios::in);

  if (!(m_file.good())) {
    msg_Error()<<"ERROR: opening "<<m_filename<<endl;
    return 0;
  }
  return 1;
}


void IO_Handler::SetSeparator(char c)
{
  m_seps.clear();
  m_seps.push_back(c);
}

void IO_Handler::AddSeparator(char c)
{
  m_seps.push_back(c);
}

void IO_Handler::SetComment(char c)
{
  m_coms.clear();
  m_coms.push_back(c);
}

void IO_Handler::AddComment(char c)
{
  m_coms.push_back(c);
}

    
// output file (compare rpa, etc.)
template <class Type> 
IO_Handler & IO_Handler::operator<<(const Type & value) {
  m_file<<" filename = "<<m_filename<<endl;
  m_file<<value;

  return *this;
}

template <class Type> 
void IO_Handler::MatrixOutput(const std::string name,Type ** const  values,const int nx, const int ny) {
  if (name!=std::string("")) 
    m_file<<" "<<name<<" = "<<endl;

  m_file<<"["<<nx<<";"<<ny<<"]";
  m_file<<"{";
  if (nx>0) ArrayOutput("", values[0],ny,0);
  for (int i=1;i<nx;++i) {
    m_file<<";"<<endl;
    ArrayOutput("", values[i],ny,0);
  }
  m_file<<"}"<<endl;
  m_nx=nx;
  m_ny=ny;  
}

template <class Type> 
void IO_Handler::ArrayOutput(const std::string name,const Type * values,const int nx, bool writesize) {
  if (name!=std::string("")) 
    m_file<<" "<<name<<" = "<<endl;

  if (writesize) m_file<<"["<<nx<<"]";
  m_file<<"{";
  if (nx>0) m_file<<values[0];
  for (int i=1;i<nx;++i) {
    if (i%10==0)
      m_file<<";"<<endl<<values[i];
    else
      m_file<<";"<<values[i];
  }
  m_file<<"}";
  if (writesize) {
    m_file<<endl;
    m_nx=nx;
  }
}
template <class Type> 
Type * IO_Handler::ArrayInput(const std::string name,int nx) {
  MyStrStream str; 
//   if (m_buffer.length()==0) {
//     getline(m_file,m_buffer); 
//     if (m_buffer.length()==0) {
//       getline(m_file,m_buffer); 
//     }
//   }
  do {
    if (m_buffer.length()==0) 
      getline(m_file,m_buffer); 
    for (unsigned int i=0; i<m_coms.size();++i) {
      int beg = m_buffer.find(m_coms[i]);
      if (beg >=0) {
	m_buffer=m_buffer.substr(0,beg);
      }
    }
  } while (m_buffer.length()==0);

  if (nx<0) {
    int beg = m_buffer.find("[");
    int end = m_buffer.find("]");
    if (beg==-1 || end==-1) {
      nx=0;
    }
    else {
      string ssize = m_buffer.substr(beg+1,end-1);
      str<<ssize;
      str>>nx;
    }
  }

  Type * values = new Type[nx];

  int x=0;
  string::iterator sit1=m_buffer.begin();
  string::iterator sit2=find(sit1,m_buffer.end(),'{');
  string::iterator send=find(sit1,m_buffer.end(),'}');
  sit1=++sit2;
  sit2=send;
  for (unsigned int i=0;i<m_seps.size();++i) {
    string::iterator si=find(sit1,send,m_seps[i]);
    if (si<sit2) sit2=si;
  }
  for (;x<nx ;++x) {
    MyStrStream helpstr;
    string value(sit1,sit2);
    helpstr<<value;
    helpstr>>values[x];

    sit1=++sit2;
    if ((sit1==send)) {
      getline(m_file,m_buffer); 
      sit1=m_buffer.begin();
      send=find(sit1,m_buffer.end(),'}');
    }
    sit2=send;
    for (unsigned int i=0;i<m_seps.size();++i) {
      string::iterator si=find(sit1,send,m_seps[i]);
      if (si<sit2) sit2=si;
    }
  }
  m_buffer=string(++sit2,m_buffer.end());

  m_nx=nx;  

  return values;
}

template <class Type> 
Type ** IO_Handler::MatrixInput(const std::string name,int nx, int ny) {
  MyStrStream str; 
  do {
    if (m_buffer.length()==0) 
      getline(m_file,m_buffer); 
    for (unsigned int i=0; i<m_coms.size();++i) {
      int beg = m_buffer.find(m_coms[i]);
      if (beg >=0) {
	m_buffer=m_buffer.substr(0,beg);
      }
    }
  } while (m_buffer.length()==0);
  if (nx<0) {
    int beg = m_buffer.find("[");
    int end = m_buffer.find("]");
    if (beg==-1 || end==-1) {
      cout<<" Error size not found "<<endl;
      nx=0; ny=0;
    }
    else {
      string ssize = m_buffer.substr(beg+1,end-1);
      int hit = ssize.find(";");
      str<<ssize.substr(0,hit);
      str>>nx;
      str<<ssize.substr(hit+1);
      str>>ny;
    }
  }

  int  hit=m_buffer.find('{');
  m_buffer=m_buffer.substr(hit+1);

  Type ** m = new Type*[nx];
  for(int i=0;i<nx;++i) {
    m[i]=ArrayInput<Type>("",ny);
    if (i<nx-1) getline(m_file,m_buffer);
  }

  m_buffer=string("");

  m_nx=nx;  
  m_ny=ny;  

  return m;
}

template <class Type> 
void IO_Handler::Output(const std::string name,const Type & value) {
 if (name!=std::string("")) 
   m_file<<" "<<name<<" = "<<value<<endl;
 else
   m_file<<value<<endl;
}

template <class Type> 
Type IO_Handler::Input(const std::string name) {
  if (name!=std::string("")) {}
 else {
  MyStrStream str; 
  do {
    if (m_buffer.length()==0) 
      getline(m_file,m_buffer); 
    for (unsigned int i=0; i<m_coms.size();++i) {
      int beg = m_buffer.find(m_coms[i]);
      if (beg >=0) {
	m_buffer=m_buffer.substr(0,beg);
      }
    }
  } while (m_buffer.length()==0);

  str<<m_buffer;
  m_buffer=std::string("");
  Type value;
  str>>value;
  return value;

  /*
   Type value;
   m_file>>value;
   return value;
  */
 }
 return (Type)0.0;
}

template <class Type> 
int IO_Handler::ValueInput(std::string name, Type & value) {
  if (m_vars.size()==0) {
    // create variable map
    for (int i=0;m_file;++i) {       
      getline(m_file,m_buffer);
      FillIn(m_buffer);
    }
  } 

  // looking for name
  Variable_Map::const_iterator cit=m_vars.find(name);
  if (cit==m_vars.end()) {
    return 0;
  } 
  else {
    std::string svalue = m_vars[name];
    if (svalue.length()==0) {
      return 0;
    }
    MyStrStream str;  
    
    str<<svalue;
    str>>value;
    return 1; 
  }

  value=0;
}
		      
void IO_Handler::FillIn(const std::string & m_buffer) {
  int hit = m_buffer.find(std::string("="));
  if (hit!=-1) {
    std::string name = m_buffer.substr(0,hit);
    Shorten(name);
    std::string value = m_buffer.substr(hit+1);
    Shorten(value);
    m_vars[name]=value;
  } 
}


void IO_Handler::Shorten(std::string& str) {
  for (;;) {    
    if (int(str[0])==32 || int(str[0])==9) str = str.substr(1);
    else break;
  }
  for (;;) {    
    if (int(str[str.length()-1])==32 ||
	//Tabulator
	int(str[str.length()-1])==9) str = str.substr(0,str.length()-1);
    else break;
  }
}


// read in class from file 
template <class Type> 
IO_Handler & IO_Handler::operator>>(Type & value) {
  m_file>>value;
  return *this;
}

template IO_Handler & IO_Handler::operator<< (const double &);
template IO_Handler & IO_Handler::operator>> (double &);
template int IO_Handler::ValueInput(std::string,double &);
// int
template void IO_Handler::Output(const std::string,const int &);
template void IO_Handler::ArrayOutput(const std::string,const int *,const int ,bool);
template void IO_Handler::MatrixOutput(const std::string ,int ** const ,const int ,const int );
template int IO_Handler::Input<int>(const std::string);
template int  * IO_Handler::ArrayInput<int>(const std::string, int);
template int ** IO_Handler::MatrixInput<int>(const std::string, int, int);
// double
template void IO_Handler::Output(const std::string,const double &);
template void IO_Handler::ArrayOutput(const std::string,const double *,const int ,bool);
template void IO_Handler::MatrixOutput(const std::string ,double ** const ,const int ,const int );
template double IO_Handler::Input<double>(const std::string);
template double  * IO_Handler::ArrayInput<double>(const std::string, int);
template double ** IO_Handler::MatrixInput<double>(const std::string, int, int);
// complex
template void IO_Handler::Output(const std::string,const Complex &);
template void IO_Handler::ArrayOutput(const std::string,const Complex *,const int ,bool);
template void IO_Handler::MatrixOutput(const std::string ,Complex ** const ,const int ,const int );
template Complex IO_Handler::Input<Complex>(const std::string);
template Complex  * IO_Handler::ArrayInput<Complex>(const std::string, int);
template Complex ** IO_Handler::MatrixInput<Complex>(const std::string, int, int);


//template <> void IO_Handler::operator<< (const Histogram &);
