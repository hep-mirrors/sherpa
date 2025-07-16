#ifndef ATOOLS_Math_Scaling_H
#define ATOOLS_Math_Scaling_H

#include "ATOOLS/Org/Getter_Function.H"
#include <vector>

namespace ATOOLS {

  template <class ValueType>
  class Scaling_Base {
  public:

    typedef ValueType Value_Type;

  protected:

    std::string m_name;

  public:

    // destructor
    virtual ~Scaling_Base();

    // member functions
    virtual Value_Type operator()(const Value_Type &x);
    virtual Value_Type operator[](const Value_Type &y);

    static void ShowScalings(const int mode=1);

    const std::string Name() const;

  };// end of class Scaling_Base
  
  typedef Getter_Function<Scaling_Base<double>,std::string> Scaling_Getter;

}// end of namespace ATOOLS

#endif
