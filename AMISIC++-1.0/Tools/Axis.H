#ifndef Axis_H
#define Axis_H

#include "Variable.H"
#include "Scaling.H"

namespace ATOOLS {

  template <class Value_Type>
  class Axis {
  public:

    enum ScalingModeID { Reference =  1,
			 Identical =  2,
			 Unknown   = 99 };

  public:

    typedef Value_Type ValueType;

  private:
    
    ATOOLS::Variable m_variable;
    ScalingModeID    m_scalingmode;

    ATOOLS::Scaling_Base<ValueType> *p_scaling;

  public:
    
    // constructor
    Axis();
    
    // destructor
    ~Axis();

    // setters
    void SetVariable(const ATOOLS::Variable _m_variable);
    void SetScalingMode(const ScalingModeID _m_scalingmode);

    void SetScaling(ATOOLS::Scaling_Base<ValueType> *_p_scaling);
    void SetScaling(const std::string scalename);

    // getters
    const ATOOLS::Variable Variable() const;
    const ScalingModeID    ScalingMode() const;

    ATOOLS::Scaling_Base<ValueType> *Scaling();

    // member functions
    ValueType DisplayedValue(ValueType realvalue,ScalingModeID tempsmode);
    ValueType RealValue(ValueType displayedvalue,ScalingModeID tempsmode);

    ValueType operator()(ValueType realvalue);
    ValueType operator[](ValueType displayedvalue);
    
  }; // end of class Axis

} // end of namespace ATOOLS

#include "Axis.IC"
#include "Axis.C"

#endif
