#ifndef Fixed_Variable_Channel_H
#define Fixed_Variable_Channel_H

#include "Variable.H"
#include "Channel_Interface.H"

namespace PHASIC {

  class Fixed_Variable_Channel: public PHASIC::Channel_Interface {
  private:
    
    ATOOLS::Variable m_variable;
    double           m_value;
    std::string      m_chid;

    bool m_trigger;
    
  public:
    
    // constructor
    Fixed_Variable_Channel(int nin,int nout,ATOOLS::Flavour *fl,
			   ATOOLS::Variable variable);
    
    // member functions
    void GeneratePoint(ATOOLS::Vec4D *p,double *ran);
    void GenerateWeight(ATOOLS::Vec4D *p);

    std::string ChID();
    
    // inline functions
    inline ATOOLS::Variable Variable() const { return m_variable; }
    inline double           Value() const    { return m_value;    }
    inline bool             Trigger() const  { return m_trigger;  }
    
    inline void SetVariable(const ATOOLS::Variable variable) { m_variable=variable; }
    inline void SetValue(const double value)                 { m_value=value;       }
    inline void SetTrigger(const bool trigger)               { m_trigger=trigger;   }
    
    inline void SetChID(const std::string chid) { m_chid=chid; }
    
  }; // end of class Fixed_Variable_Channel
  
} // end of namespace PHASIC

#endif
