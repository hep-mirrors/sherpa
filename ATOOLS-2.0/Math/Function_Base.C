#include "Function_Base.H"

using namespace ATOOLS;

Function_Base::~Function_Base() {}

void Function_Base::SetDefault(double _m_defval)    
{ m_defval=_m_defval; }

void Function_Base::SetType(std::string _m_type) 
{ m_type=_m_type; }

void Function_Base::SetParameters(double *parameters)    
{ return; }

double Function_Base::GetValue(double x)        
{ return (*this)(x); }   

double Function_Base::operator()(double x)         
{ return m_defval; }

double Function_Base::operator()()               
{ return m_defval; }

std::string Function_Base::Type()                     
{ return m_type; }

