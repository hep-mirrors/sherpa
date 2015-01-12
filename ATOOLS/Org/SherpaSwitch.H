#ifndef SHERPASWITCH_H
#define SHERPASWITCH_H
#include <iostream>
#include <iomanip>
#include <typeinfo>
#include <string>
#include <vector>
#include <cmath>
#include "ATOOLS/Org/ParameterBase.H"

// The actual Switch class to be used
template <typename T>
class SherpaSwitch : public ParameterBase<T> {
    public:
        SherpaSwitch(std::string name, const T& cur, const T& def, std::vector<T>
                allowed, const bool allowChangeAfterLock=false) :
            ParameterBase<T>(name, cur, def, allowChangeAfterLock) {
                m_allowed=allowed;
    
    }
    bool isAllowed() {
        bool dec=false;
        for (unsigned int i=0;i<m_allowed.size();i++) {
            if (this->get() == m_allowed[i]) dec=true;
            break;
        }
        return dec;
    }
    virtual bool isDefault() {
        return (this->get() == this->getDefault());
    }
    void print() {
        std::cout << "Switch '" << this->name() << "' (type '" << std::string(typeid(T).name()) << "') current value: "  << this->get();
        if (this->isDefault()) { std::cout << " at default value --- "; }
        else { std::cout << " (default: " << this->getDefault() << ")"; }
        if (this->isAllowed()) { std::cout << " IS allowed setting ( "; }
        else { std::cout << " NOT an allowed setting! ( "; }
        for (unsigned int i=0;i<m_allowed.size();i++) { std::cout<<m_allowed[i] << " "; }
        std::cout << ")" << std::endl;
    }
    private:
    std::vector<T> m_allowed;
};
#endif
