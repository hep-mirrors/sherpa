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
}// end of namespace ATOOLS

  template <> double ATOOLS::NotDefined<double>();   
  template <> int ATOOLS::NotDefined<int>();      
  template <> long ATOOLS::NotDefined<long int>(); 
  template <> std::string ATOOLS::NotDefined<std::string>(); 
  template <> ATOOLS::Switch::code ATOOLS::NotDefined<ATOOLS::Switch::code>(); 
  template <> ATOOLS::Beam_Type::code ATOOLS::NotDefined<ATOOLS::Beam_Type::code>(); 
  template <> ATOOLS::ISR_Type::code ATOOLS::NotDefined<ATOOLS::ISR_Type::code>(); 
  template <> ATOOLS::String_Type::code ATOOLS::NotDefined<ATOOLS::String_Type::code>(); 
  template <> ATOOLS::Model_Type::code ATOOLS::NotDefined<ATOOLS::Model_Type::code>(); 
  template <> ATOOLS::Flavour ATOOLS::NotDefined<ATOOLS::Flavour>(); 


#endif 
