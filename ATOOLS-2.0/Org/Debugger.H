#ifndef Debugger_H
#define Debugger_H

#include "Scaling.H"
#include <iostream>
#include <string>
#include <vector>
#include <map>

namespace ATOOLS {

  template <typename ValueType>
  class Value_Carrier;

  class Value_Carrier_Base {
  protected:

    std::string m_type;

  public:

    // constructor
    Value_Carrier_Base();
  
    // destructor
    virtual ~Value_Carrier_Base();
  
    // member functions
    template <typename ValueType>
    Value_Carrier<ValueType> *UpCast();

    template <typename ValueType> ValueType *Address();
    template <typename ValueType> ValueType  Value();

    virtual void PrintAddress(std::ostream &ostr) const = 0;
    virtual void PrintValue(std::ostream &ostr) const = 0;

    const std::string &Type() const { return m_type; }

  };// end of class Value_Carrier_Base

  template <typename ValueType>
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

    static unsigned int s_print;

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
    static void PrintStatus(const std::string &type="",const std::string &function="");
    static void PrintMethodInfo(const std::string &method);
    static void PrintLocalInfo(const std::string &method,const std::string &info);

    static void SetOutputMode(const unsigned int print);

    template <typename ValueType> void AddValue(const std::string &spec,ValueType &value);

    // inline functions
    inline std::map<std::string,Value_Carrier_Base*> &Values() { return m_values; }

    inline Value_Carrier_Base *Value(const std::string &spec) { return m_values[spec]; }

    template <typename ValueType> inline ValueType *Address(const std::string &spec)
    { return m_values[spec]->Address<ValueType>(); }
    template <typename ValueType> inline ValueType  Value(const std::string &spec)
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
#define dbg ATOOLS::dbg
#define SET_DEBUG_MODE dbg.SetOutputMode(std::numeric_limits<unsigned int>::max())
#endif

#define DEBUG_INFO dbg.PrintStatus("",__PRETTY_FUNCTION__)
#define SPECIAL_INFO(INFO) dbg.PrintStatus(INFO,__PRETTY_FUNCTION__)
#define METHOD_INFO dbg.PrintMethodInfo(__PRETTY_FUNCTION__)
#define LOCAL_INFO(INFO) dbg.PrintLocalInfo(__FUNCTION__,ToString(INFO))

#include "Debugger.inl.H"

#endif
