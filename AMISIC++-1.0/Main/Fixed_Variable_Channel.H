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

    public:

      // constructor
      Fixed_Variable_Channel(int _nin,int _nout,ATOOLS::Flavour *_fl,
			     ATOOLS::Variable _m_variable);

      // setter
      void SetVariable(const ATOOLS::Variable _m_variable);
      void SetValue(const ValueType _m_value);

      // getter
      const ATOOLS::Variable Variable() const;
      const ValueType        Value() const;

      // member functions
      void GeneratePoint(ATOOLS::Vec4D *_p,double *_ran);
      void GenerateWeight(ATOOLS::Vec4D *_p);

    }; // end of class Fixed_Variable_Channel
    
} // end of namespace PHASIC

#include "Fixed_Variable_Channel.IC"
#include "Fixed_Variable_Channel.C"

#endif
