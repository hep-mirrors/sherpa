#include <algorithm>
#include "IO_Handler.H"
#include "MyStrStream.H"
#include "Data_Return.H"
#include "MyComplex.H"
#include "Message.H"

using namespace AORGTOOLS;
using namespace std;

IO_Handler::IO_Handler() { filename=std::string(""); };

IO_Handler::~IO_Handler() { 
  if (!(filename==std::string(""))) {
    file.close();
  }
};
    
// set output filename
int IO_Handler::SetFileName(std::string _name) {
  if (!(filename==std::string(""))) {
    file.close();
  }
  filename=_name;
  msg.Tracking()<<" opened file "<<filename<<endl;
  file.open(filename.c_str(),ios::out);

  if (!(file.good())) {
    msg.Tracking()<<" ERROR: opening "<<filename<<endl;
    return 0;
  }
  file.precision(15);
  msg.Tracking()<<" done "<<endl;
  return 1;
}
    
// set input filename
int IO_Handler::SetFileNameRO(string _name) {
  if (!(filename==std::string(""))) {
    msg.Tracking()<<"before closing"<<endl;
    msg.Tracking()<<file;
    file.close();
  }
  filename=_name;
  msg.Tracking()<<" opened file "<<filename<<endl;
  file.open(filename.c_str(),ios::in);

  if (!(file.good())) {
    msg.Error()<<" ERROR: opening "<<filename<<endl;
    return 0;
  }
  msg.Tracking()<<" done "<<endl;
  return 1;
}
    
// output file (compare rpa, etc.)
template <class Type> 
IO_Handler & IO_Handler::operator<<(const Type & value) {
  file<<" filename = "<<filename<<endl;
  file<<value;

  return *this;
}

template <class Type> 
void IO_Handler::MatrixOutput(const std::string name,Type ** const  values,const int nx, const int ny) {
  msg.Tracking()<<" output "<<name<<endl;
  if (name!=std::string("")) 
    file<<" "<<name<<" = "<<endl;

  file<<"["<<nx<<";"<<ny<<"]";
  file<<"{";
  if (nx>0) ArrayOutput("", values[0],ny,0);
  for (int i=1;i<nx;++i) {
    file<<";"<<endl;
    ArrayOutput("", values[i],ny,0);
  }
  file<<"}"<<endl;
  
}

template <class Type> 
void IO_Handler::ArrayOutput(const std::string name,const Type * values,const int nx, bool writesize) {
  if (name!=std::string("")) 
    file<<" "<<name<<" = "<<endl;

  if (writesize) file<<"["<<nx<<"]";
  file<<"{";
  if (nx>0) file<<values[0];
  for (int i=1;i<nx;++i) {
    if (i%10==0)
      file<<";"<<endl<<values[i];
    else
      file<<";"<<values[i];
  }
  file<<"}";
  if (writesize) file<<endl;
  
}
template <class Type> 
Type * IO_Handler::ArrayInput(const std::string name,int nx) {
  MyStrStream str; 
  if (buffer.length()==0) {
    getline(file,buffer); 
    if (buffer.length()==0) {
      getline(file,buffer); 
    }
  }
  if (nx<0) {
    int beg = buffer.find("[");
    int end = buffer.find("]");
    if (beg==-1 || end==-1) {
      msg.Tracking()<<" Error size not fount "<<endl;
      nx=0;
    }
    else {
      string ssize = buffer.substr(beg+1,end-1);
      msg.Tracking()<<" ssize="<<ssize<<endl;
      str<<ssize;
      str>>nx;
    }
  }

  Type * values = new Type[nx];

  int x=0;
  string::iterator sit1=buffer.begin();
  string::iterator sit2=find(sit1,buffer.end(),'{');
  string::iterator send=find(sit1,buffer.end(),'}');
  sit1=++sit2;
  sit2=find(sit1,send,';');
  for (;x<nx ;++x) {
    MyStrStream helpstr;
    string value(sit1,sit2);
    helpstr<<value;
    helpstr>>values[x];

    sit1=++sit2;
    if ((sit1==send)) {
      getline(file,buffer); 
      sit1=buffer.begin();
      send=find(sit1,buffer.end(),'}');
    }
    sit2=find(sit1,send,';');
  }
  buffer=string(++sit2,buffer.end());
  return values;
}

template <class Type> 
Type ** IO_Handler::MatrixInput(const std::string name,int nx, int ny) {
  MyStrStream str;  
  getline(file,buffer); 
  if (nx<0) {
    int beg = buffer.find("[");
    int end = buffer.find("]");
    if (beg==-1 || end==-1) {
      msg.Tracking()<<" Error size not found "<<endl;
      nx=0; ny=0;
    }
    else {
      string ssize = buffer.substr(beg+1,end-1);
      msg.Tracking()<<" ssize="<<ssize<<endl;
      int hit = ssize.find(";");
      msg.Tracking()<<" hit="<<hit<<endl;
      str<<ssize.substr(0,hit);
      str>>nx;
      str<<ssize.substr(hit+1);
      str>>ny;
      msg.Tracking()<<" nx="<<nx<<" ny="<<ny<<endl;
    }
  }

  int  hit=buffer.find('{');
  buffer=buffer.substr(hit+1);

  Type ** m = new Type*[nx];
  for(int i=0;i<nx;++i) {
    m[i]=ArrayInput<Type>("",ny);
    if (i<nx-1) getline(file,buffer);
  }

  buffer=string("");
  return m;
}

template <class Type> 
void IO_Handler::Output(const std::string name,const Type & value) {
 if (name!=std::string("")) 
   file<<" "<<name<<" = "<<value<<endl;
 else
   file<<value<<endl;
}

template <class Type> 
Type IO_Handler::Input(const std::string name) {
 if (name!=std::string("")) 
   msg.Tracking()<<" "<<name<<" =  ?????????"<<endl;
 else {
   Type value;
   file>>value;
   return value;
 }
}

template <class Type> 
int IO_Handler::ValueInput(std::string name, Type & value) {
  if (vars.size()==0) {
    // create variable map
    msg.Tracking()<<file.gcount()<<endl;
    for (int i=0;file;++i) {       
      getline(file,buffer);
      msg.Tracking()<<i<<"## "<<buffer<<" ##"<<endl;
      FillIn(buffer);
    }
  } 

  // looking for name
  Variable_Map::const_iterator cit=vars.find(name);
  if (cit==vars.end()) {
    value =  NotDefined<Type>();
    return 0;
  } 
  else {
    std::string svalue = vars[name];
    if (svalue.length()==0) {
      value=NotDefined<Type>();
      return 0;
    }
    MyStrStream str;  
    
    str<<svalue;
    str>>value;
    return 1; 
  }

  value=0;
}
		      
void IO_Handler::FillIn(const std::string & buffer) {
  int hit = buffer.find(std::string("="));
  if (hit!=-1) {
    std::string name = buffer.substr(0,hit);
    Shorten(name);
    std::string value = buffer.substr(hit+1);
    Shorten(value);
    vars[name]=value;
  } 
};


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


// readin class from file 
template <class Type> 
 IO_Handler &  IO_Handler::operator>>(Type & value) {
  file>>value;
}

template IO_Handler & IO_Handler::operator<< (const double &);
template IO_Handler & IO_Handler::operator>> (double &);
template int IO_Handler::ValueInput(std::string,double &);
// int
template void IO_Handler::Output(const std::string,const int &);
template void IO_Handler::MatrixOutput(const std::string ,int ** const ,const int ,const int );
template int IO_Handler::Input<int>(const std::string);
template int ** IO_Handler::MatrixInput<int>(const std::string, int, int);
// double
template void IO_Handler::Output(const std::string,const double &);
template void IO_Handler::MatrixOutput(const std::string ,double ** const ,const int ,const int );
template double IO_Handler::Input<double>(const std::string);
template double ** IO_Handler::MatrixInput<double>(const std::string, int, int);
// complex
template void IO_Handler::Output(const std::string,const Complex &);
template void IO_Handler::MatrixOutput(const std::string ,Complex ** const ,const int ,const int );
template Complex IO_Handler::Input<Complex>(const std::string);
template Complex ** IO_Handler::MatrixInput<Complex>(const std::string, int, int);


//template <> void IO_Handler::operator<< (const AMATOOLS::Histogram &);
