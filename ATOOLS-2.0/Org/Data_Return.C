#include "Data_Return.H"

using namespace ATOOLS;

template <> double NotDefined<double>()   
{ 
  return std::numeric_limits<double>::max(); 
} 

template <> int NotDefined<int>()      
{ 
  return std::numeric_limits<int>::max(); 
} 

template <> long NotDefined<long int>() 
{ 
  return std::numeric_limits<long int>::max(); 
} 

template <> std::string NotDefined<std::string>() 
{ 
  return std::string("Unknown parameter"); 
} 

template <> Switch::code NotDefined<Switch::code>() 
{ 
  return Switch::Unknown; 
} 

template <> Beam_Type::code NotDefined<Beam_Type::code>() 
{ 
  return Beam_Type::Unknown; 
} 

template <> ISR_Type::code NotDefined<ISR_Type::code>() 
{ 
  return ISR_Type::Unknown; 
} 

template <> String_Type::code NotDefined<String_Type::code>() 
{ 
  return String_Type::Unknown; 
} 

template <> Model_Type::code NotDefined<Model_Type::code>() 
{ 
  return Model_Type::Unknown; 
} 

template <> Flavour NotDefined<Flavour>() 
{ 
  return Flavour(kf::none); 
} 

