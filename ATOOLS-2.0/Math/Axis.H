#ifndef Axis_H
#define Axis_H

#include "Smart_Pointer.H"
#include "Variable.H"
#include "Scaling.H"

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

    SP(Variable_Base<Value_Type>) p_variable;
    SP(Scaling_Base<Value_Type>)  p_scaling;

  public:
    
    // constructor
    Axis();
    Axis(const Axis &ref);
    
    // destructor
    ~Axis();

    // member functions
    void SetScaling(const std::string &scalename);
    void SetVariable(const std::string &variablename);

    Value_Type DisplayedValue(const Value_Type &realvalue,
			      ScalingModeID tempsmode);
    Value_Type RealValue(const Value_Type &displayedvalue,
			 ScalingModeID tempsmode);

    // inline functions
    inline void SetScalingMode(const ScalingModeID &scalingmode)
    { m_scalingmode=scalingmode; }

    inline void SetScaling(Scaling_Base<Value_Type> *const scaling)
    { p_scaling=scaling; }
    inline void SetVariable(Variable_Base<Value_Type> *const variable)
    { p_variable=variable; }

    inline ScalingModeID ScalingMode() const { return m_scalingmode; }

    inline SP(Scaling_Base<Value_Type>) &Scaling()   { return p_scaling;  }
    inline SP(Variable_Base<Value_Type>) &Variable() { return p_variable; }

    inline Value_Type operator()(const Value_Type &realvalue)
    { return DisplayedValue(realvalue,Unknown); }
    inline Value_Type operator[](const Value_Type &displayedvalue)
    { return RealValue(displayedvalue,Unknown); }
    
  }; // end of class Axis

}// end of namespace ATOOLS

#endif
