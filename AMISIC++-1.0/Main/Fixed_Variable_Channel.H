#ifndef Fixed_Variable_Channel_H
#define Fixed_Variable_Channel_H

#include "Variable.H"
#include "Channel_Interface.H"

namespace PHASIC {

  class Fixed_Variable_Channel: public PHASIC::Channel_Interface {
  private:
    
    ATOOLS::Variable m_variable;
    double           m_value;

    bool m_trigger;
    
  public:
    
    // constructor
    Fixed_Variable_Channel(int nin,int nout,ATOOLS::Flavour *fl,
			   ATOOLS::Variable variable);
    
    // member functions
    void GeneratePoint(ATOOLS::Vec4D *p,double *ran);
    void GenerateWeight(ATOOLS::Vec4D *p);
    
    // inline functions
    ATOOLS::Variable Variable() const { return m_variable; }
    double           Value() const    { return m_value;    }
    bool             Trigger() const  { return m_trigger;  }
    
    void SetVariable(const ATOOLS::Variable variable) { m_variable=variable; }
    void SetValue(const double value)                 { m_value=value;       }
    void SetTrigger(const bool trigger)               { m_trigger=trigger;   }
    
  }; // end of class Fixed_Variable_Channel
  
} // end of namespace PHASIC

#endif
