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
    const std::string &Name() const;

    int SelectorID() const;

    Value_Type operator()(const Vec3D *vectors) const;
    Value_Type operator()(const Vec4D *vectors) const;

    virtual Value_Type Value(const Vec3D *vectors) const;
    virtual Value_Type Value(const Vec4D *vectors) const;

    static void ShowVariables(const int mode=1);

  };

  typedef Getter_Function<Variable_Base<double>,std::string> Variable_Getter;

} // end of namespace ATOOLS

#endif
