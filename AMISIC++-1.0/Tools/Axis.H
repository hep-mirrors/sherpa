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

  //////////////////
  // Doxygen part //
  //////////////////

  /*
    \class Axis
    \brief Provides the properties of an axis in a graph
    
    The purpose of the class Axis is to provide other classes with the
    basic properties of an axis, as there are 

    \itemize
      \item the variable which is plotted and 
      \item the scaling which is used
    \enditemize

    Therefore Axis employs the classes Variable and Scaling.
  */

  /*
    \enum Axis::ScalingModeID
    \brief Defines a scaling mode

    There are two different types of scaling modes, namely

    \itemize
      \item Reference and 
      \item Identical
    \enditemize

    where Reference means that the current Scaling_Base is employed
    for the scaling of values on the axis and Identical means 
    that no scaling occurs.
  */

  /*
    \fn Axis::operator()
    \brief Performs the scaling of values
  */

  /*
    \fn Axis::operator[]
    \brief Performs the inverse scaling of values
  */

  /*
    \fn Axis::SetVariable(const ATOOLS::Variable _m_variable)
    \brief Sets the variable of the axis
  */

  /*
    \fn Axis::SetScalingMode(const ScalingModeID _m_scalingmode)
    \brief Sets the scaling mode
  */

  /*
    \fn Axis::SetScaling(const std::string scalename)
    \brief Sets the scaling function according to a tag
  */

  /*
    \fn Axis::Variable() const
    \brief Returns the variable of the axis
  */
  
  /*
    \fn Axis::ScalingMode() const
    \brief Returns the scaling mode
  */

  /*
    \fn Axis::Scaling()
    \brief Returns the current Scaling_Base
  */

} // end of namespace ATOOLS

#include "Axis.IC"
#include "Axis.C"

#endif
