#ifndef Variable_H
#define Variable_H

#include "Getter_Function.H"
#include <string>

namespace ATOOLS {

  class Vec3D;
  class Vec4D;

  template <class ValueType>
  class Variable_Base {
  public:

    typedef ValueType Value_Type;

  protected:
    
    std::string m_name;

    int m_selectorid;

  public:

    // constructors
    Variable_Base(const std::string &name);

    // destructor
    virtual ~Variable_Base();

    // member functions
    virtual Value_Type Value(const Vec3D *vectors);
    virtual Value_Type Value(const Vec4D *vectors);

    static void ShowVariables(const int mode=1);

    // inline functions
    inline const std::string &Name() const { return m_name; }

    inline int SelectorID() const { return m_selectorid; }

    inline Value_Type operator()(const Vec3D *vectors) 
    { return Value(vectors); }
    inline Value_Type operator()(const Vec4D *vectors) 
    { return Value(vectors); }

  };

  typedef Getter_Function<Variable_Base<double>,std::string> Variable_Getter;

} // end of namespace ATOOLS

#endif
