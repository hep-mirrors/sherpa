#include "Variable.H"

using namespace ATOOLS;
  
Variable::Variable()
{ SetType(Unknown); }

Variable::Variable(const TypeID _m_type)
{ SetType(_m_type); }

Variable::Variable(const unsigned int _m_type)
{ SetType((TypeID)_m_type); }

Variable::Variable(const std::string _m_type)
{ 
  SetType(StringToType(_m_type)); 
  if (_m_type!=std::string("")) SetName(_m_type);
}

std::string Variable::TypeToString(const TypeID type)
{ 
  switch (type) {
  case fourvec_r : return std::string("r_\\mu");
  case r_t       : return std::string("r_t");
  case r_x       : return std::string("r_x");
  case r_y       : return std::string("r_y");
  case r_z       : return std::string("r_z");
  case vec_r     : return std::string("\\vec{r}");
  case r_perp    : return std::string("r_\\perp");
  case r_abs     : return std::string("|\\vec{r}|");
  case r_abs_sqr : return std::string("r_\\mu^2");
  case fourvec_p : return std::string("p_\\mu");
  case p_t       : return std::string("p_t");
  case p_x       : return std::string("p_x");
  case p_y       : return std::string("p_y");
  case p_z       : return std::string("p_z");
  case vec_p     : return std::string("\\vec{p}");
  case p_perp    : return std::string("p_\\perp");
  case p_abs     : return std::string("|\\vec{p}|");
  case p_abs_sqr : return std::string("p_\\mu^2");
  case E         : return std::string("E");
  case E_z       : return std::string("E_z");
  case vec_E     : return std::string("\\vec{E}");
  case E_perp    : return std::string("E_\\perp");
  case m         : return std::string("m");
  case m_perp    : return std::string("m_\\perp");
  case s         : return std::string("s");
  case s_prime   : return std::string("s'");
  case s_i       : return std::string("s_i");
  case s_ij      : return std::string("s_{ij}");
  case y         : return std::string("y");
  case eta       : return std::string("\\eta");
  case angle     : return std::string("angle");
  case theta     : return std::string("\\theta");
  case theta_i   : return std::string("\\theta_i");
  case theta_ij  : return std::string("\\theta_{ij}");
  case phi       : return std::string("\\phi");
  case phi_ij    : return std::string("\\phi_{ij}");
  case Unknown   : return std::string("Unknown");
  }
  return std::string("");
}

Variable::TypeID Variable::StringToType(const std::string type)
{ 
  if (type=="r_\\mu")       return fourvec_r;
  if (type=="r_t")          return r_t;
  if (type=="r_x")          return r_x;
  if (type=="r_y")          return r_y;
  if (type=="r_z")          return r_z;
  if (type=="\\vec{r}")     return vec_r;
  if (type=="r_\\perp")     return r_perp;
  if (type=="|\\vec{r}|")   return r_abs;
  if (type=="r_\\mu^2")     return r_abs_sqr;
  if (type=="p_\\mu")       return fourvec_p;
  if (type=="p_t")          return p_t;
  if (type=="p_x")          return p_x;
  if (type=="p_y")          return p_y;
  if (type=="p_z")          return p_z;
  if (type=="\\vec{p}")     return vec_p;
  if (type=="p_\\perp")     return p_perp;
  if (type=="|\\vec{p}|")   return p_abs;
  if (type=="p_\\mu^2")     return p_abs_sqr;
  if (type=="E")            return E;
  if (type=="E_z")          return E_z;
  if (type=="\\vec{E}")     return vec_E;
  if (type=="E_\\perp")     return E_perp;
  if (type=="m")            return m;
  if (type=="m_\\perp")     return m_perp;
  if (type=="s")            return s;
  if (type=="s'")           return s_prime;
  if (type=="s_i")          return s_i;
  if (type=="s_{ij}")       return s_ij;
  if (type=="y")            return y;
  if (type=="\\eta")        return eta;
  if (type=="angle")        return angle;
  if (type=="\\theta")      return theta;
  if (type=="\\theta_i")    return theta_i;
  if (type=="\\theta_{ij}") return theta_ij;
  if (type=="\\phi")        return phi;
  if (type=="\\phi_{ij}")   return phi_ij;
  if (type=="Unknown")      return Unknown;
  return Unknown;
}

int Variable::TypeToSelectorID(const TypeID type)
{ 
  switch (type) {
  case p_perp    : return 12;
  case E         : return 11;
  case E_perp    : return 15;
  case m         : return 21;
  case y         : return 13;
  case eta       : return 16;
  case theta_i   : return 14;
  case theta_ij  : return 22; 
  default: return -1;
  }
  return -1;
}

Variable::TypeID Variable::SelectorIDToType(const int selectorid)
{ 
  switch (selectorid) {
  case 12: return p_perp;
  case 11: return E;
  case 15: return E_perp;
  case 21: return m;
  case 13: return y;
  case 16: return eta;
  case 14: return theta_i;
  case 22: return theta_ij;
  default: return Unknown;
  }
  return Unknown;
}

void Variable::SetName(const std::string _m_name)
{ m_name=_m_name; }

const std::string Variable::Name() const
{ return m_name; }

void Variable::SetType(const TypeID _m_type)
{ SetName(TypeToString(m_type=_m_type)); }

const Variable::TypeID Variable::Type() const
{ return m_type; }
  
unsigned int Variable::TypeToInt(const TypeID type)
{ return (int)type; }

Variable::TypeID Variable::IntToType(const unsigned int type)
{ return (TypeID)type; }




