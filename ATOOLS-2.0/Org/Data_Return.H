#ifndef Data_Return_H
#define Data_Return_H

#include "Data_Read.H"
#include <iostream>


namespace ATOOLS {

#define UNDEF_FLOAT  __FLT_MAX__
#define UNDEF_DOUBLE __DBL_MAX__
#define UNDEF_INT    __INT_MAX__
#define UNDEF_LONG   __LONG_MAX__

  template <class Type>
  inline Type NotDefined() { return (Type)0; }

  template <class Type>
  inline Type IsDefined(Type d,const char * name) {
    if (d==NotDefined<Type>()) {
      std::cout<<"Variable "<<name<<" not yet defined!"<<std::endl;
      std::cout<<"check parameter file"<<std::endl;
      std::cerr<<"*Error: "<<name<<" not defined in parameter file!!"<<std::endl;
      abort();
    } 
    return d;
  }

  template <> inline double             NotDefined<double>() { return UNDEF_DOUBLE; } 
  template <> inline int                NotDefined<int>()    { return UNDEF_INT; } 
  template <> inline long               NotDefined<long>()   { return UNDEF_LONG; } 

  template <> inline std::string        NotDefined<std::string>() { 
    return std::string("Unknown parameter"); 
  } 
  
  template <> inline Switch::code       NotDefined<Switch::code>() { 
    return Switch::Unknown; 
  } 

  template <> inline Beam_Type::code     NotDefined<Beam_Type::code>() { 
    return Beam_Type::Unknown; 
  } 

  template <> inline ISR_Type::code     NotDefined<ISR_Type::code>() { 
    return ISR_Type::Unknown; 
  } 

  template <> inline String_Type::code  NotDefined<String_Type::code>() { 
    return String_Type::Unknown; 
  } 

  template <> inline Model_Type::code   NotDefined<Model_Type::code>() { 
    return Model_Type::Unknown; 
  } 

  template <> inline Flavour NotDefined<Flavour>() { 
    return Flavour(kf::none); 
  } 


} // namespace ATOOLS

#endif // Data_Return_H
