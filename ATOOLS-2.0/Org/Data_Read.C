#include "Data_Read.H"
#include "Data_Return.H"
#include "Message.H"
#include "MyStrStream.H"

using namespace AORGTOOLS;  
using namespace APHYTOOLS;
using namespace std;

void Data_Read::SetValue(std::string name, std::string value) {
  Shorten(name);
  // define value
  Shorten(value);
  // insert name-value pair in list
  parameters[name]=value;

  cout<<" parameter "<<name<<" = "<<value<<endl;
}

// definition
template <class Type>
Type  Data_Read::GetValue(std::string name, Type default_value) {
  Shorten(name);
  if (name.length()==0) {
    msg.Events()<<"Could not find any value for empty name. Return "<<default_value<<"."<<endl;
    return default_value;
  }
  Type dummy = GetValue<Type>(name);
  if (dummy!=NotDefined<Type>()) { return dummy; }
  msg.Events()<<"Could not find any allowed value for "<<name<<". Return "<<default_value<<"."<<endl;
  return default_value;
}

template <class Type>
Type  Data_Read::GetValue(std::string name) {
  Shorten(name);
  Type invar;
  if (name.length()==0) return NotDefined<Type>();

  Parameter_Map::const_iterator cit=parameters.find(name);
  if (cit==parameters.end()) return  NotDefined<Type>();

  std::string value = parameters[name];
  if (value.length()==0) return NotDefined<Type>();
  
  MyStrStream str;      
  str<<value;
  str>>invar;
  return invar;
}

Data_Read::Data_Read(std::string filename) { ReadIn(filename); }


void Data_Read::FillIn(char * dummy) {
  if (dummy[0]!='!' && strlen(dummy)>0) {
    std::string buffer(dummy);
    int hit = buffer.find(std::string("="));
    if (hit!=-1) {
      // define name
      std::string name = buffer.substr(0,hit);
      Shorten(name);
      // define value
      std::string value = buffer.substr(hit+1);
      int hit = value.find(std::string("!"));
      if (hit!=-1) value = value.substr(0,hit);
      Shorten(value);
      // insert name-value pair in list
      parameters[name]=value;
    }
  }
}

void Data_Read::ReadIn(std::string filename) {
  std::ifstream file;
  file.open(filename.c_str());
  if (file.bad()) {
    msg.Error()<< " ERROR: opening " << filename <<endl;
    exit (-1);
  }
  char dummy[256];
      
  for (;file;) {
    file.getline(dummy,256);
    FillIn(dummy); 
  }
  file.close();
}





// definition  (specialisation), explicit instanciation
template <> std::string Data_Read::GetValue<std::string>(std::string name) {
  Shorten(name);
  Parameter_Map::const_iterator cit=parameters.find(name);
  if (cit==parameters.end()) return  NotDefined<std::string>();

  std::string value = parameters[name];
  std::string invar;

  if (value.length()==0)         return NotDefined<std::string>();
  /*
  MyStrStream str;      
  str<<value;
  str>>invar;
  if (value!=invar) cout<<" Warning:"<<value<<" -> "<<invar<<endl;
  return invar;
  */
  return value;
}


template <> Switch::code Data_Read::GetValue<Switch::code>(std::string name) {
  Shorten(name);
  Parameter_Map::const_iterator cit=parameters.find(name);
  if (cit==parameters.end()) return  NotDefined<Switch::code>();

  std::string value = parameters[name];

  if (value.length()==0)         return NotDefined<Switch::code>();
  if (value==std::string("On"))  return Switch::On;
  if (value==std::string("Off")) return Switch::Off;
  
  msg.Error()<<" unknown Switch::code "<<name<<" = "<<value<<" !!!"<<endl;
  return NotDefined<Switch::code>();
}

// Beams
template <> Beam_Type::code Data_Read::GetValue<Beam_Type::code>(std::string name) {
  Shorten(name);
  Parameter_Map::const_iterator cit=parameters.find(name);
  if (cit==parameters.end()) return  NotDefined<Beam_Type::code>();
  std::string value = parameters[name];
  
  if (value==std::string("Monochromatic"))        return Beam_Type::Monochromatic;    
  if (value==std::string("Gaussian"))             return Beam_Type::Gaussian;    
  if (value==std::string("Laser_Backscattering")) return Beam_Type::Laser_Back;    
    
  msg.Error()<<"Unknown Beam type  "<<name<<" = "<<value<<" !!!"<<endl;
  return NotDefined<Beam_Type::code>();
}

template <> Beam_Generator::code Data_Read::GetValue<Beam_Generator::code>(std::string name) {
  Shorten(name);
  Parameter_Map::const_iterator cit=parameters.find(name);
  if (cit==parameters.end()) return  NotDefined<Beam_Generator::code>();
  std::string value = parameters[name];
  
  if (value==std::string("Internal"))  return Beam_Generator::Internal;
    
  msg.Error()<<"Unknown Beam generator  "<<name<<" = "<<value<<" !!!"<<endl;
  return NotDefined<Beam_Generator::code>();
}




template <> Beam_Shape::code Data_Read::GetValue<Beam_Shape::code>(std::string name) {
  Shorten(name);
  Parameter_Map::const_iterator cit=parameters.find(name);
  if (cit==parameters.end()) return  NotDefined<Beam_Shape::code>();
  std::string value = parameters[name];
  
  if (value==std::string("Cylinder"))          return Beam_Shape::Cylinder;
  if (value==std::string("Gaussian_Cylinder")) return Beam_Shape::Gaussian_Cylinder;
  
  msg.Error()<<"Unknown Beam shape  "<<name<<" = "<<value<<" !!!"<<endl;
  return NotDefined<Beam_Shape::code>();
}



template <> ISR_Type::code Data_Read::GetValue<ISR_Type::code>(std::string name) {
  Shorten(name);
  Parameter_Map::const_iterator cit=parameters.find(name);
  if (cit==parameters.end()) return  NotDefined<ISR_Type::code>();

  std::string value = parameters[name];
  
  if (value==std::string("No"))                   return ISR_Type::No;    
  if (value==std::string("simple Struct"))        return ISR_Type::Simple_Struc;    
  if (value==std::string("extended Struct"))      return ISR_Type::Extended_Struc;    
  if (value==std::string("Projection"))           return ISR_Type::PA;    
  if (value==std::string("mod. Proj."))           return ISR_Type::MPA;    
  if (value==std::string("Extrapolation"))        return ISR_Type::EA;    
  if (value==std::string("KoralZ"))               return ISR_Type::KoralZ;    
  if (value==std::string("Pythia"))               return ISR_Type::Pythia;   
  if (value==std::string("KKMC"))                 return ISR_Type::KKMC;    
  
  msg.Error()<<"Unknown ISR type  "<<name<<" = "<<value<<" !!!"<<endl;
  return NotDefined<ISR_Type::code>();
}


template <> String_Type::code Data_Read::GetValue<String_Type::code>(std::string name) {
  Shorten(name);
  Parameter_Map::const_iterator cit=parameters.find(name);
  if (cit==parameters.end()) return  NotDefined<String_Type::code>();

  std::string value = parameters[name];

  if (value==std::string("NoString")) return String_Type::NoString;    
  if (value==std::string("String"))   return String_Type::String;    
  if (value==std::string("Library"))  return String_Type::Library;    

  msg.Error()<<"Unknown String_Type  "<<name<<" = "<<value<<" !!!"<<endl;
  return NotDefined<String_Type::code>();
}

// definition (specialisation), explicit instanciation
template <> Model_Type::code Data_Read::GetValue<Model_Type::code>(std::string name) {
  Shorten(name);
  Parameter_Map::const_iterator cit=parameters.find(name);
  if (cit==parameters.end()) return  NotDefined<Model_Type::code>();

  std::string value = parameters[name];
  if (value==std::string("pure_QCD")) return Model_Type::pure_QCD;
  if (value==std::string("QCD"))      return Model_Type::QCD;
  if (value==std::string("pure_EW"))  return Model_Type::pure_EW;
  if (value==std::string("SM"))       return Model_Type::SM;
  if (value==std::string("MSSM"))     return Model_Type::MSSM;
  if (value==std::string("THDM"))     return Model_Type::THDM;
  if (value==std::string("ADD"))      return Model_Type::ADD;

  msg.Error()<<"Unknown Model "<<name<<" = "<<value<<" !!!"<<endl;
  return NotDefined<Model_Type::code>();
}

template <>  Flavour Data_Read::GetValue<Flavour>(std::string name) {
  Shorten(name);
  
  Parameter_Map::const_iterator cit=parameters.find(name);
  if (cit==parameters.end()) return  NotDefined<Flavour>();
  if (!kf_table.IsInitialised()) {
    msg.Error()<<"Warning in Flavour Data_Read::GetValue."<<endl
	       <<"   kf table not initialized yet. Return undefined flavour."<<endl;
    return NotDefined<Flavour>();
  }

  std::string value = parameters[name];
  // compare strings
  bool anti= 0;  // 0 = particle;  1 = anti-particle
                 // looking for "anti-" statement
  int hit = value.find(std::string("anti-"));
  if (hit!=-1) {
    value = value.substr(hit+5);
    anti=1;
  }
  
  kf::code kfc;
  kfc = kf_table.FromString(value);
  if (kfc!=kf::none) return Flavour(kfc);
  else {
    // looking for "+" / "-"
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
    // here anti should always be true!!!!
    if (anti) return (Flavour(kfc).Bar());
    return Flavour(kfc);
  }
  
  
  // check if number
  int kfci = GetValue<int>(name);
  Flavour fl = Flavour(kf::code(abs(kfci)));
  if (kfci<0) fl = fl.Bar();
  return fl;
}
 
int Data_Read::Crossfoot(string name) {
  int sum = 0;
  for (int i=0;i<name.length();++i)
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
  for (Parameter_Iterator it = parameters.begin(); it!=parameters.end() ; ++it) {
    int id= ((Crossfoot(it->first)&0xff)*Crossfoot(it->second))&0xffffff;
    sum^=id;
  }
  sum=(sum&0xffffff); //|(0x1000000*(parameters.size() & 255));
//   cout.setf(ios::hex, ios::basefield);  
//   cout<<" globalID="<<sum<<endl;
//   cout.setf(ios::dec, ios::basefield); 

  MyStrStream str;  
  string key;
  str.setf(ios::hex, ios::basefield);
  str<<(parameters.size() & 255);
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
  msg.Events()<<"rpa_id="<<key<<endl;
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
  for (Parameter_Iterator it = parameters.begin(); it!=parameters.end() ; ++it) {
    file<<" "<<it->first<<" = "<<it->second<<endl;
  }
  file.close();
}

void Data_Read::Shorten(std::string& str) {
  //kill initial spaces
  for (;;) {    
    if (int(str[0])==32 || int(str[0])==9) str = str.substr(1);
    else break;
  }
  //kill final spaces
  for (;;) {    
    if (int(str[str.length()-1])==32 ||
	//Tabulator
	int(str[str.length()-1])==9) str = str.substr(0,str.length()-1);
    else break;
  }
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


