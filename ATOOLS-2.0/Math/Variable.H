#ifndef Variable_H
#define Variable_H

#include <string>

namespace ATOOLS {

  class Variable {
  public:
    
    enum TypeID { fourvec_r = 100,
		  r_t       = 110,
		  r_x       = 111,
		  r_y       = 112,
		  r_z       = 113,
		  vec_r     = 120,
		  r_perp    = 130,
		  r_abs     = 140,
		  r_abs_sqr = 150,
		  fourvec_p = 200,
		  p_t       = 210,
		  p_x       = 211,
		  p_y       = 212,
		  p_z       = 213,
		  vec_p     = 220,
		  p_perp    = 230,
		  p_abs     = 240,
		  p_abs_sqr = 250,
		  E         = 300,
		  E_z       = 313,
		  vec_E     = 320,
		  E_perp    = 330,
		  m         = 400,
		  m_perp    = 430,
		  s         = 500,
		  s_prime   = 510,
		  s_i       = 520,
		  s_ij      = 530,
		  y         = 600,
		  eta       = 610,
		  angle     = 700,
		  theta     = 710,
		  theta_i   = 711,
		  theta_ij  = 712,
		  phi       = 720,
		  phi_ij    = 722,
		  Unknown   = 999 };

  private:
    
    TypeID m_type;

    std::string m_name;

  public:

    //constructors
    Variable();
    Variable(const TypeID _m_type);
    Variable(const unsigned int _m_type);
    Variable(const std::string _m_type);

    // setters
    void SetType(const TypeID _m_type);
    void SetName(const std::string _m_name);

    // getters
    const TypeID      Type() const;
    const std::string Name() const;

    // member functions
    unsigned int TypeToInt(const TypeID type);
    std::string  TypeToString(const TypeID type);

    TypeID IntToType(const unsigned int type);
    TypeID StringToType(const std::string type);

    static int    TypeToSelectorID(const TypeID type);
    static TypeID SelectorIDToType(const int selectorid);

  };

} // end of namespace ATOOLS

// #include "Variable.IC"

#endif
