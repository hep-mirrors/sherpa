#ifndef Variable_H
#define Variable_H

#include "Getter_Function.H"
#include "Vector.H"
#include <string>

namespace ATOOLS {


  template <class ValueType>
  class Variable_Base {
  public:

    typedef ValueType Value_Type;

  protected:
    
    std::string m_name, m_idname;

    int m_selectorid;

  public:

    // constructors
    Variable_Base(const std::string &name,
		  const std::string &idname="");

    // destructor
    virtual ~Variable_Base();

    // member functions
    const std::string &Name() const;
    const std::string &IDName() const;

    int SelectorID() const;

    Value_Type operator()(const Vec3D *vectors,const int &n=1) const;
    Value_Type operator()(const Vec4D *vectors,const int &n=1) const;

    virtual Value_Type Value(const Vec3D *vectors,const int &n=1) const;
    virtual Value_Type Value(const Vec4D *vectors,const int &n=1) const;

    static void ShowVariables(const int mode=1);

  };

  typedef Getter_Function<Variable_Base<double>,std::string> Variable_Getter;

} // end of namespace ATOOLS

#endif
