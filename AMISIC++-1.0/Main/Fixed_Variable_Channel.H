#ifndef Fixed_Variable_Channel_H
#define Fixed_Variable_Channel_H

#include "Channel_Interface.H"

namespace PHASIC {

    template <class Value_Type>
    class Fixed_Variable_Channel: public PHASIC::Channel_Interface {
    public:

      typedef Value_Type ValueType;

    private:
      
      ATOOLS::Variable m_variable;
      ValueType        m_value;

      bool m_trigger;

    public:

      // constructor
      Fixed_Variable_Channel(int _nin,int _nout,ATOOLS::Flavour *_fl,
			     ATOOLS::Variable _m_variable);

      // member functions
      void GeneratePoint(ATOOLS::Vec4D *_p,double *_ran);
      void GenerateWeight(ATOOLS::Vec4D *_p);

      // inline functions
      inline ATOOLS::Variable Variable() const { return m_variable; }
      inline ValueType        Value() const    { return m_value;    }
      inline bool             Trigger() const  { return m_trigger;  }

      inline void SetVariable(const ATOOLS::Variable variable)   
      { m_variable=variable; }
      inline void SetValue(const ValueType _m_value) 
      { m_value=_m_value; }
      inline void SetTrigger(const bool trigger) 
      { m_trigger=trigger; }

    }; // end of class Fixed_Variable_Channel
    
} // end of namespace PHASIC

#include "Fixed_Variable_Channel.C"

#endif
