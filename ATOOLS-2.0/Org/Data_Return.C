#include "Data_Return.H"

using namespace ATOOLS;

template <> double ATOOLS::NotDefined<double>()   
{ 
  return std::numeric_limits<double>::max(); 
} 

template <> int ATOOLS::NotDefined<int>()      
{ 
  return std::numeric_limits<int>::max(); 
} 

template <> long ATOOLS::NotDefined<long int>() 
{ 
  return std::numeric_limits<long int>::max(); 
} 

template <> std::string ATOOLS::NotDefined<std::string>() 
{ 
  return std::string("Unknown parameter"); 
} 

template <> Switch::code ATOOLS::NotDefined<Switch::code>() 
{ 
  return Switch::Unknown; 
} 

template <> Beam_Type::code ATOOLS::NotDefined<Beam_Type::code>() 
{ 
  return Beam_Type::Unknown; 
} 

template <> ISR_Type::code ATOOLS::NotDefined<ISR_Type::code>() 
{ 
  return ISR_Type::Unknown; 
} 

template <> String_Type::code ATOOLS::NotDefined<String_Type::code>() 
{ 
  return String_Type::Unknown; 
} 

template <> Model_Type::code ATOOLS::NotDefined<Model_Type::code>() 
{ 
  return Model_Type::Unknown; 
} 

template <> Flavour ATOOLS::NotDefined<Flavour>() 
{ 
  return Flavour(kf::none); 
} 

