#ifndef Axis_H
#define Axis_H

#include "Variable.H"
#include "Scaling.H"

namespace ATOOLS {

  template <class Value_Type>
  class Axis {
  public:

    enum ScalingMode { Reference =  1,
		       Identical =  2,
		       Unknown   = 99 };

  public:

    typedef Value_Type ValueType;

  private:
    
    ATOOLS::Variable m_variable;
    ScalingMode      m_scalingmode;

    ATOOLS::Scaling_Base<ValueType> *p_scaling;

  public:
    
    // constructor
    Axis();
    
    // destructor
    ~Axis();

    // setters
    void SetVariable(ATOOLS::Variable _m_variable);
    void SetScalingMode(ScalingMode _m_scalingmode);

    void SetScaling(ATOOLS::Scaling_Base<ValueType> *_p_scaling);

    // getters
    ATOOLS::Variable GetVariable();
    ScalingMode      GetScalingMode();

    ATOOLS::Scaling_Base<ValueType> *GetScaling();

    // member functions
    ValueType DisplayedValue(ValueType realvalue,ScalingMode tempsmode);
    ValueType RealValue(ValueType displayedvalue,ScalingMode tempsmode);

    ValueType operator()(ValueType realvalue);
    ValueType operator[](ValueType displayedvalue);
    
  }; // end of class Axis

} // end of namespace ATOOLS

#include "Axis.IC"
#include "Axis.C"

#endif
