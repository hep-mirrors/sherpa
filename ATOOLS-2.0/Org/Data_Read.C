#include "Data_Read.H"
#include "Data_Return.H"
#include "Message.H"
#include "MyStrStream.H"
#include "Exception.H"
#include "Type.H"
#include <iomanip>

using namespace ATOOLS;
using namespace std;
// static
Parameter_Map ATOOLS::Data_Read::s_commandlineparameters;

template <class Type> 
const Type Data_Read::ReturnData(const std::string &name,const Type type) 
{
  ATOOLS::msg.LogFile()<<name<<" \t= \t"<<type<<" \t! "
		       <<"Type<"<<ATOOLS::Type::GetType(type)<<"> "<<std::endl;
  return type;
}

void Data_Read::SetValue(std::string name, std::string value) {
  Shorten(name);
  Shorten(value);
  m_parameters[name]=value;
}

// definition
template <class Type>
Type  Data_Read::GetValue(std::string name, Type default_value) {
  Shorten(name);
  if (name.length()==0) {
    msg.LogFile()<<"Could not find any value for empty name. Return "<<default_value<<"."<<endl;
    return ReturnData(name,default_value);
  }
  Type dummy = GetValue<Type>(name);
  if (dummy!=NotDefined<Type>()) { 
    return ReturnData(name,dummy); 
  }
  msg.LogFile()<<"WARNING: Could not find any allowed value for "<<name
	       <<". Return "<<default_value<<"."<<endl;
  MyStrStream str;      
  std::string default_value_str;
  str<<default_value;
  str>>default_value_str;

  m_parameters[name]=default_value_str;
  return ReturnData(name,default_value);
}

template <class Type>
Type  Data_Read::GetValue(std::string name) {
  Shorten(name);
  Type invar;
  if (name.length()==0) return ReturnData(name,NotDefined<Type>());
  Parameter_Map::const_iterator cit=m_parameters.find(name);
  if (cit==m_parameters.end()) return ReturnData(name,NotDefined<Type>());
  std::string value = m_parameters[name];
  if (value.length()==0) return ReturnData(name,NotDefined<Type>());
  MyStrStream str;      
  str<<value;
  str>>invar;
  return ReturnData(name,invar);
}

Data_Read::Data_Read(std::string filename, bool ignoremissingfile) { 
  m_fileexists=true;
  ReadIn(filename,ignoremissingfile); 
}


void Data_Read::FillIn(std::string buffer) {
  if (buffer.length()>0 && buffer[0]!='!') {
    int hit = buffer.find(std::string("="));
    if (hit!=-1) {
      std::string name = buffer.substr(0,hit);
      Shorten(name);
      std::string value = buffer.substr(hit+1);
      int hit = value.find(std::string("!"));
      if (hit!=-1) value = value.substr(0,hit);
      Shorten(value);
      m_parameters[name]=value;
    }
  }
}

void Data_Read::ReadIn(std::string filename, bool ignoremissingfile) {
  std::ifstream file;
  file.open(filename.c_str());
  if (!file.good()) {
    if (ignoremissingfile) {
      msg_Tracking()<<" WARNING parameter file "<<filename<<" does not exist ! "<<std::endl;
      m_fileexists=false;
    }
    else {
    throw(Exception(ex::critical_error,std::string("Cannot open file '")+filename+std::string("'"),
		    "Data_Read","ReadIn"));
    }
  }
  std::string dummy;
      
  for (;file;) {
    getline(file,dummy);
    FillIn(dummy); 
  }
  file.close();

  AddCommandLine();
}

void Data_Read::SetCommandLine(std::string name, std::string value)
{
  Shorten(name);
  s_commandlineparameters[name]=value;
}

void Data_Read::AddCommandLine()
{
  for (Parameter_Iterator it = s_commandlineparameters.begin(); it!=s_commandlineparameters.end() ; ++it) {
    m_parameters[it->first]=it->second;
  }
}

// definition  (specialisation), explicit instanciation
template <> std::string Data_Read::GetValue<std::string>(std::string name) 
{
  Shorten(name);
  std::string invar;
  Parameter_Map::const_iterator cit=m_parameters.find(name);
  if (cit==m_parameters.end()) return ReturnData(name,NotDefined<std::string>());
  std::string value=m_parameters[name];
  if (value.length()==0) return ReturnData(name,NotDefined<std::string>());
  return ReturnData(name,value);
}

template <> Switch::code Data_Read::GetValue<Switch::code>(std::string name) {
  Shorten(name);
  Parameter_Map::const_iterator cit=m_parameters.find(name);
  if (cit==m_parameters.end()) return ReturnData(name,NotDefined<Switch::code>());
  std::string value = m_parameters[name];
  if (value.length()==0) return ReturnData(name,NotDefined<Switch::code>());
  if (value==std::string("On")) return ReturnData(name,Switch::On);
  if (value==std::string("Off")) return ReturnData(name,Switch::Off);
  msg.Error()<<"Error in Data_Read::GetValue<Switch::code>:"<<endl
	     <<"   Unknown Switch::code "<<name<<" = "<<value<<"."<<endl;
  return ReturnData(name,NotDefined<Switch::code>());
}

// Beams
template <> Beam_Type::code Data_Read::GetValue<Beam_Type::code>(std::string name) {
  Shorten(name);
  Parameter_Map::const_iterator cit=m_parameters.find(name);
  if (cit==m_parameters.end()) return ReturnData(name,NotDefined<Beam_Type::code>());
  std::string value = m_parameters[name];
  
  if (value==std::string("Monochromatic")) return ReturnData(name,Beam_Type::Monochromatic); 
  if (value==std::string("Gaussian")) return ReturnData(name,Beam_Type::Gaussian);    
  if (value==std::string("Laser_Backscattering")) return ReturnData(name,Beam_Type::Laser_Back);    
  if (value==std::string("Spectrum_Reader")) return ReturnData(name,Beam_Type::Spec_Read);    
  if (value==std::string("Simple_Compton")) return ReturnData(name,Beam_Type::Simple_Compton);
  msg.Error()<<"Error in Data_Read::GetValue<Beam_Type::code>:"<<endl
	     <<"   Unknown Beam type  "<<name<<" = "<<value<<"."<<endl;
  return ReturnData(name,NotDefined<Beam_Type::code>());
}

template <> Beam_Generator::code Data_Read::GetValue<Beam_Generator::code>(std::string name) {
  Shorten(name);
  Parameter_Map::const_iterator cit=m_parameters.find(name);
  if (cit==m_parameters.end()) return ReturnData(name,NotDefined<Beam_Generator::code>());
  std::string value = m_parameters[name];
  if (value==std::string("Internal")) return ReturnData(name,Beam_Generator::Internal);
  msg.Error()<<"Error in Data_Read::GetValue<Beam_Generator::code>:"<<endl
	     <<"Unknown Beam generator  "<<name<<" = "<<value<<"."<<endl;
  return ReturnData(name,NotDefined<Beam_Generator::code>());
}

template <> Beam_Shape::code Data_Read::GetValue<Beam_Shape::code>(std::string name) {
  Shorten(name);
  Parameter_Map::const_iterator cit=m_parameters.find(name);
  if (cit==m_parameters.end()) return ReturnData(name,NotDefined<Beam_Shape::code>());
  std::string value = m_parameters[name];
  if (value==std::string("Cylinder")) return ReturnData(name,Beam_Shape::Cylinder);
  if (value==std::string("Gaussian_Cylinder")) return ReturnData(name,Beam_Shape::Gaussian_Cylinder);
  msg.Error()<<"Error in Data_Read::GetValue<Beam_Shape::code>:"<<endl
	     <<"Unknown Beam shape  "<<name<<" = "<<value<<" !!!"<<endl;
  return ReturnData(name,NotDefined<Beam_Shape::code>());
}

template <> ISR_Type::code Data_Read::GetValue<ISR_Type::code>(std::string name) {
  Shorten(name);
  Parameter_Map::const_iterator cit=m_parameters.find(name);
  if (cit==m_parameters.end()) return ReturnData(name,NotDefined<ISR_Type::code>());

  std::string value = m_parameters[name];
  
  if (value==std::string("No"))              return ReturnData(name,ISR_Type::No);    
  if (value==std::string("simple Struct"))   return ReturnData(name,ISR_Type::Simple_Struc);    
  if (value==std::string("extended Struct")) return ReturnData(name,ISR_Type::Extended_Struc);    
  if (value==std::string("Projection"))      return ReturnData(name,ISR_Type::PA);    
  if (value==std::string("mod. Proj."))      return ReturnData(name,ISR_Type::MPA);    
  if (value==std::string("Extrapolation"))   return ReturnData(name,ISR_Type::EA);    
  if (value==std::string("KoralZ"))          return ReturnData(name,ISR_Type::KoralZ);    
  if (value==std::string("Pythia"))          return ReturnData(name,ISR_Type::Pythia);   
  if (value==std::string("KKMC"))            return ReturnData(name,ISR_Type::KKMC);    
  
  msg.Error()<<"Error in Data_Read::GetValue<ISR_Type::code>:"<<endl
	     <<"   Unknown ISR type  "<<name<<" = "<<value<<" !!!"<<endl;
  return ReturnData(name,NotDefined<ISR_Type::code>());
}

template <> String_Type::code Data_Read::GetValue<String_Type::code>(std::string name) {
  Shorten(name);
  Parameter_Map::const_iterator cit=m_parameters.find(name);
  if (cit==m_parameters.end()) return ReturnData(name,NotDefined<String_Type::code>());

  std::string value = m_parameters[name];

  if (value==std::string("NoString")) return ReturnData(name,String_Type::NoString);    
  if (value==std::string("String"))   return ReturnData(name,String_Type::String);    
  if (value==std::string("Library"))  return ReturnData(name,String_Type::Library);    

  msg.Error()<<"Error in Data_Read::GetValue<String_Type::code>:"<<endl
	     <<"   Unknown String_Type  "<<name<<" = "<<value<<" !!!"<<endl;
  return ReturnData(name,NotDefined<String_Type::code>());
}

template <> Model_Type::code Data_Read::GetValue<Model_Type::code>(std::string name) {
  Shorten(name);
  Parameter_Map::const_iterator cit=m_parameters.find(name);
  if (cit==m_parameters.end()) return ReturnData(name,NotDefined<Model_Type::code>());
  std::string value = m_parameters[name];
  if (value==std::string("pure_QCD")) return ReturnData(name,Model_Type::pure_QCD);
  if (value==std::string("QCD"))      return ReturnData(name,Model_Type::QCD);
  if (value==std::string("pure_EW"))  return ReturnData(name,Model_Type::pure_EW);
  if (value==std::string("SM"))       return ReturnData(name,Model_Type::SM);
  if (value==std::string("MSSM"))     return ReturnData(name,Model_Type::MSSM);
  if (value==std::string("THDM"))     return ReturnData(name,Model_Type::THDM);
  if (value==std::string("ADD"))      return ReturnData(name,Model_Type::ADD);
  if (value==std::string("SMHL"))     return ReturnData(name,Model_Type::SMHL);
  msg.Error()<<"Error in Data_Read::GetValue<Model_Type::code>:"<<endl
	     <<"   Unknown Model "<<name<<" = "<<value<<" !!!"<<endl;
  return ReturnData(name,NotDefined<Model_Type::code>());
}

template <>  Flavour Data_Read::GetValue<Flavour>(std::string name) {
  Shorten(name);
  
  Parameter_Map::const_iterator cit=m_parameters.find(name);
  if (cit==m_parameters.end()) return ReturnData(name,NotDefined<Flavour>());
  if (!kf_table.IsInitialised()) {
    msg.Error()<<"Warning in Flavour Data_Read::GetValue."<<endl
	       <<"   kf table not initialized yet. Return undefined flavour."<<endl;
    return ReturnData(name,NotDefined<Flavour>());
  }

  std::string value = m_parameters[name];
  bool anti= 0;  // 0 = particle;  1 = anti-particle
                 // looking for "anti-" statement
  int hit = value.find(std::string("anti-"));
  if (hit!=-1) {
    value = value.substr(hit+5);
    anti=1;
  }
  
  kf::code kfc;
  kfc = kf_table.FromString(value);
  if (kfc!=kf::none) return ReturnData(name,Flavour(kfc));
  else {
    hit = value.find(std::string("+"));
    if (hit!=-1) {
      value[hit] = '-';
      anti       = 1;
    } 
    else {
      hit = value.find(std::string("-"));
      if (hit!=-1) {
	value[hit] = '+';
	anti       = 1;
      }
    }
  }
  kfc = kf_table.FromString(value);
  if (kfc!=kf::none) {
    if (anti) return ReturnData(name,(Flavour(kfc).Bar()));
    return Flavour(kfc);
  }
  int kfci = GetValue<int>(name);
  Flavour fl = Flavour(kf::code(abs(kfci)));
  if (kfci<0) fl = fl.Bar();
  return ReturnData(name,fl);
}
 
int Data_Read::Crossfoot(string name) {
  int sum = 0;
  for (size_t i=0;i<name.length();++i)
    sum+=int(name[i]);
  return sum;
}

string Data_Read::GenerateKey() {
  // Hexadeximal Number
  //  nn_m_dddddd
  //  v v  version main and subnumber of code (not included yet)
  //  nn  number of parameters readin
  //  m   model number
  //  dddddd = Sum_i (cross sum of name * cross sum of value above)
  //  output hex,

  // possible extensions: 
  // * output 0-9,A-Z stat hex (130 times as many states))
  // * additional rotation by position
  // * make partitions of parameters that have been used and those that have not been used.

  int sum=0;
  for (Parameter_Iterator it = m_parameters.begin(); it!=m_parameters.end() ; ++it) {
    int id= ((Crossfoot(it->first)&0xff)*Crossfoot(it->second))&0xffffff;
    sum^=id;
  }
  sum=(sum&0xffffff); //|(0x1000000*(parameters.size() & 255));

  MyStrStream str;  
  string key;
  str.setf(ios::hex, ios::basefield);
  str<<(m_parameters.size() & 255);
  Model_Type::code m = GetValue<Model_Type::code>("MODEL");
  str<<"_";
  if (m!=Model_Type::Unknown) {
    str<<int(m);
  }
  else {
    str<<"_";
  }
  str<<"_";
  str<<sum;

  str>>key;
  return key;
}

void Data_Read::WriteOut(std::string filename,int flag) {
  std::fstream file;

#ifdef __GNUC__
#if __GNUC__ > 2 
  // GNU gcc 3.x.x C++ Compiler
  std::_Ios_Openmode flagc = _Ios_Openmode(flag);
  file.open(filename.c_str(),flagc);
#else
  // GNU gcc 2.95.x C++ Compiler
  file.open(filename.c_str(),flag);
#endif
#else
  // All others
  file.open(filename.c_str(),flag);
#endif
  
  // add a header
  file<<"!======================================== "<<endl;
  file<<"! File: "<<filename<<endl;
  file<<"!======================================== "<<endl;
  
  // write out map content
  for (Parameter_Iterator it = m_parameters.begin(); it!=m_parameters.end() ; ++it) {
    file<<" "<<it->first<<" = "<<it->second<<endl;
  }
  file.close();
}

void Data_Read::Shorten(std::string& str) {
  //kill initial spaces
  for (;str.length()>0;) {    
    if (int(str[0])==32 || int(str[0])==9) str = str.substr(1);
    else break;
  }
  //kill final spaces
  for (;str.length()>0;) {    
    if (int(str[str.length()-1])==32 ||
	//Tabulator
	int(str[str.length()-1])==9) str = str.substr(0,str.length()-1);
    else break;
  }
  size_t pos=std::string::npos;
  while ((pos=str.find(" "))!=std::string::npos) str[pos]='_'; 
}

// explicit instanciation for standard types:
template int              Data_Read::GetValue<int>(std::string);
template long             Data_Read::GetValue<long>(std::string);
template float            Data_Read::GetValue<float>(std::string);
template double           Data_Read::GetValue<double>(std::string);
//template std::string      Data_Read::GetValue<std::string>(std::string);

template int              Data_Read::GetValue<int>(std::string,int);
template long             Data_Read::GetValue<long>(std::string,long);
template float            Data_Read::GetValue<float>(std::string,float);
template double           Data_Read::GetValue<double>(std::string,double);
template std::string      Data_Read::GetValue<std::string>(std::string,std::string);
template Switch::code     Data_Read::GetValue<Switch::code>(std::string,Switch::code);


