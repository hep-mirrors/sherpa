#ifndef AMISIC_Tools_Fixed_Variable_Channel_H
#define AMISIC_Tools_Fixed_Variable_Channel_H
/*!
  \file Fixed_Variable_Channel.H
  \brief Declares the class Fixed_Variable_Channel
*/


#include "ATOOLS/Math/Variable.H"
#include "PHASIC++/Channels/Channel_Interface.H"

namespace PHASIC {

  class Fixed_Variable_Channel: public PHASIC::Channel_Interface {
  private:
    
    ATOOLS::Variable_Base<double> *p_variable;

    double      m_value;
    std::string m_chid;

    bool m_trigger;
    
  public:
    
    // constructor
    Fixed_Variable_Channel(int nin,int nout,ATOOLS::Flavour *fl,
			   const std::string &variable);
    
    // member functions
    void GeneratePoint(ATOOLS::Vec4D *p,double *ran);
    void GenerateWeight(ATOOLS::Vec4D *p);

    std::string ChID();
    
    // inline functions
    inline ATOOLS::Variable_Base<double> *Variable() const { return p_variable; }

    inline double Value() const    { return m_value;    }
    inline bool   Trigger() const  { return m_trigger;  }
    
    inline void SetVariable(ATOOLS::Variable_Base<double> *const variable) 
    { p_variable=variable; }

    inline void SetValue(const double value)   { m_value=value;     }
    inline void SetTrigger(const bool trigger) { m_trigger=trigger; }
    
    inline void SetChID(const std::string chid) { m_chid=chid; }
    
  }; // end of class Fixed_Variable_Channel
  /*!
    \class Fixed_Variable_Channel
    \brief Special integration channel

    This class implements a special integration channel, which
    is designed for holding one phasespace variable fixed, while
    integrating out others. It is used during the event generation 
    step with Simple_Chain to generate events with fixed transverse 
    momenta for the outgoing partons.
  */

} // end of namespace PHASIC

#endif
