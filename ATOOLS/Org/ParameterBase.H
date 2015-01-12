#ifndef SHERPA_PARAMETER_BASE_H
#define SHERPA_PARAMETER_BASE_H
#include <iostream>
#include <iomanip>
#include <typeinfo>
#include <string>
#include <cmath>
#include "ATOOLS/Org/boost/lexical_cast.hpp"

//http://stackoverflow.com/questions/13980157/c-class-with-template-member-variable
class ParameterBaseBase
{
public:
    virtual ~ParameterBaseBase() {}
    template<class T> const T& get() const; //to be implimented after Parameter
    template<class T> const T& getDefault() const; //to be implimented after Parameter
    template<class T> const T& name() const; //to be implimented after Parameter
    template<class T, class U> void setValue(const U& rhs); //to be implimented after Parameter
    virtual bool isDefault() = 0; // Abstract method
};

template <typename T>
class ParameterBase : public ParameterBaseBase
{
public:
    ParameterBase(std::string name, const T& cur, const T& def, const bool allowChangeAfterLock=false) :
       m_value(cur),
       m_default(def),
       m_name(name), m_allowChangeAfterLock(allowChangeAfterLock) {
    }
    const T& get() const {return m_value;}
    const T& getDefault() const {return m_default;}
    void setValue(const T& cur) {
        if (!m_locked || m_allowChangeAfterLock) m_value=cur;
        else {std::cout << "Warning, not allowing change of locked parameter/switch '" << name() << "'" << std::endl;}
    }    
    void setValueAndLock(const T& cur) {setValue(cur);lock();}    
    void setDefault(const T& def) {m_default=def;}
    void lock() {m_locked=true;}
    void unlock() {m_locked=false;}
    bool isLocked() {return m_locked;}
    std::string name() { return m_name; }
    //virtual bool isDefault() {
        //bool dec = true;
        //std::string ttype = std::string(typeid(T).name());
        //if (ttype == "d") {
            //if (fabs(boost::lexical_cast<double>(m_value) - boost::lexical_cast<double>(m_default)) > 1e-15) {
             //dec = false;
            //}
        //}
        //if (ttype == "f") {
            //if (fabs(boost::lexical_cast<float>(m_value) - boost::lexical_cast<float>(m_default)) > 1e-7) {
             //dec = false;
            //}
        //}
        //else {
            //dec = boost::lexical_cast<std::string>(m_value) == boost::lexical_cast<std::string>(m_default);
        //}
        //return dec;
    //};
private:
    T m_value;
    T m_default;
    bool m_locked;
    bool m_allowChangeAfterLock;
    std::string m_name;
};

//Here's the trick: dynamic_cast rather than virtual
template<class T> const T& ParameterBaseBase::get() const
{ return dynamic_cast<const ParameterBase<T>&>(*this).get(); }
template<class T> const T& ParameterBaseBase::getDefault() const
{ return dynamic_cast<const ParameterBase<T>&>(*this).getDefault(); }
template<class T, class U> void ParameterBaseBase::setValue(const U& rhs)
{ return dynamic_cast<ParameterBase<T>&>(*this).setValue(rhs); }
#endif
