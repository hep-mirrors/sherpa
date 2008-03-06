#ifndef Data_Return_H
#define Data_Return_H

#include "Data_Read.H"
#include "My_Limits.H"
#include <iostream>

namespace ATOOLS {

  template <class Type> Type NotDefined() { return (Type)0; }

  template <class Type>
  Type IsDefined(Type d,const char * name) 
  {
    if (d!=NotDefined<Type>()) return d;
    std::cout<<"Variable "<<name<<" not yet defined !\n"
	     <<"check parameter file"<<std::endl;
    std::cerr<<"*Error: "<<name<<" not defined in parameter file!!"
	     <<std::endl;
    abort();
  }

  template <> double NotDefined<double>();   
  template <> int NotDefined<int>();      
  template <> long NotDefined<long int>(); 
  template <> std::string NotDefined<std::string>(); 
  template <> Switch::code NotDefined<Switch::code>(); 
  template <> Beam_Type::code NotDefined<Beam_Type::code>(); 
  template <> ISR_Type::code NotDefined<ISR_Type::code>(); 
  template <> String_Type::code NotDefined<String_Type::code>(); 
  template <> Flavour NotDefined<Flavour>(); 

}// end of namespace ATOOLS

#endif 
