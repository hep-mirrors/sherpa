#ifndef Debugger_H
#define Debugger_H

#include <iostream>
#include <vector>
#include <map>

namespace ATOOLS {

  template <class ValueType>
  class Value_Carrier;

  class Value_Carrier_Base {
  protected:

    unsigned int m_type;

  public:
  
    // member functions
    template <class ValueType>
    Value_Carrier<ValueType> *UpCast();

    template <class ValueType> ValueType *Address();
    template <class ValueType> ValueType  Value();

    virtual void PrintAddress(std::ostream &ostr) const = 0;
    virtual void PrintValue(std::ostream &ostr) const = 0;

    unsigned int Type() const { return m_type; }

  };// end of class Value_Carrier_Base

  template <class ValueType>
  class Value_Carrier: public Value_Carrier_Base {
  public:

    typedef ValueType Value_Type;

  private:

    Value_Type *p_address;

  public:

    // constructor
    Value_Carrier(Value_Type &address);

    // member functions
    Value_Type *Address() const;
    Value_Type  Value() const;

    void PrintAddress(std::ostream &ostr) const;
    void PrintValue(std::ostream &ostr) const;

  };// end of class Value_Carrier

  class Debugger {
  private:

    static std::vector<Debugger*>  s_objects;
    static std::ostream           *s_output;

    long unsigned int m_address;
    std::string       m_name, m_type;

    std::map<std::string,Value_Carrier_Base*> m_values;

    friend class Value_Carrier_Base;

  public:
    
    // constructor
    Debugger(long unsigned int address,const std::string &name);

    // destructor
    ~Debugger();

    // members
    static void PrintStatus(const std::string &type="");

    template <class ValueType> void AddValue(const std::string &spec,ValueType &value);

    // inline functions
    inline std::map<std::string,Value_Carrier_Base*> &Values() { return m_values; }

    inline Value_Carrier_Base *Value(const std::string &spec) { return m_values[spec]; }

    template <class ValueType> inline ValueType *Address(const std::string &spec)
    { return m_values[spec]->Address<ValueType>(); }
    template <class ValueType> inline ValueType  Value(const std::string &spec)
    { return m_values[spec]->Value<ValueType>();   }
    
    inline const std::string &Type()              { return m_type; }
    inline void  SetType(const std::string &type) { m_type=type;   }

    inline static void SetOutput(std::ostream *const output) { s_output=output; }

    friend std::ostream &operator<<(std::ostream &ostr,const Debugger &debugger);

  };// end of class Debugger

  extern Debugger dbg;

}// end of namespace ATOOLS

#ifndef COMPILE__Debugger
#define Debugger() Debugger((long unsigned int)this,__PRETTY_FUNCTION__)
#endif

#include "Debugger.inl.H"

#endif
