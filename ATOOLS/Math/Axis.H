#ifndef ATOOLS_Math_Axis_H
#define ATOOLS_Math_Axis_H

#include "ATOOLS/Math/Variable.H"
#include "ATOOLS/Math/Scaling.H"

namespace ATOOLS {

  template <class ValueType>
  class Axis {
  public:

    enum ScalingModeID { 
      Reference =  1,
      Identical =  2,
      Unknown   = 99 
    };

  public:

    typedef ValueType Value_Type;

  private:
    
    ScalingModeID m_scalingmode;

    Variable_Base<Value_Type> *p_variable;
    Scaling_Base<Value_Type>  *p_scaling;

  public:
    
    // constructor
    Axis();
    Axis(const Axis &ref);
    
    // destructor
    ~Axis();

    // member functions
    void SetScaling(const std::string &scalename);
    void SetVariable(const std::string &variablename);

    void SetScaling(Scaling_Base<Value_Type> *const scaling);
    void SetVariable(Variable_Base<Value_Type> *const variable);

    const Scaling_Base<Value_Type>  *Scaling() const;
    const Variable_Base<Value_Type> *Variable() const;

    void SetScalingMode(const ScalingModeID &scalingmode);

    ScalingModeID ScalingMode() const;

    Value_Type DisplayedValue(const Value_Type &realvalue,
			      ScalingModeID tempsmode) const;
    Value_Type RealValue(const Value_Type &displayedvalue,
			 ScalingModeID tempsmode) const;

    Value_Type operator()(const Value_Type &realvalue) const;
    Value_Type operator[](const Value_Type &displayedvalue) const;
    
  }; // end of class Axis

}// end of namespace ATOOLS

#endif
