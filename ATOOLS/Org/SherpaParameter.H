#ifndef SHERPAPARAMETER_H
#define SHERPAPARAMETER_H
#include <iostream>
#include <iomanip>
#include <typeinfo>
#include <string>
#include <cmath>
#include "ATOOLS/Org/ParameterBase.H"

// The actual Parameter class to be used
template <typename T>
class SherpaParameter : public ParameterBase<T> {
    public:
        SherpaParameter(std::string name, const T& cur, const T& def, double vmin, double vmax, const bool allowChangeAfterLock=false) : ParameterBase<T>(name, cur, def, allowChangeAfterLock) {
            assert(vmin < vmax);
            m_min=vmin;
            m_max=vmax;
    }
        bool inRange() {return (this->get() >= m_min && m_max >= this->get());}
        void print() {
            std::cout << "Parameter '" << this->name() << "' (type '" << std::string(typeid(T).name()) << "') current value: " << std::setprecision(16) << this->get();
            if (this->isDefault()) {
                std::cout << " at default value --- ";
            }
            else {
                std::cout << " (default: " << std::setprecision(16) << this->getDefault() << ")";
            }
            if (inRange()) {
                std::cout <<  "is IN valid range (" << this->m_min << " <= " << this->m_max << ")";
            }
            else {
                std::cout << "NOT in valid range (" << this->m_min << " <= " << this->m_max << ")";
            }
            std::cout << std::endl;
        }
    private:
        double m_min;
        double m_max;
};
#endif
