#ifndef SHERPAVECTOR_H
#define SHERPAVECTOR_H
#include <iostream>
#include <iomanip>
#include <typeinfo>
#include <string>
#include <vector>
#include <cmath>
#include "ATOOLS/Org/ParameterBase.H"
#include "ATOOLS/Org/boost/assign/std/vector.hpp" // for 'operator+=()'

using namespace boost::assign; // Bring 'operator+=()' into scope

// Vectors, such as selectors
template <typename T>
class SherpaVector : public ParameterBase<std::vector<T> > {
    public:
        SherpaVector(std::string name, const std::vector<T>& cur, const std::vector<T>& def, const bool allowChangeAfterLock=false)
           : ParameterBase<std::vector<T> >(name, cur, def,allowChangeAfterLock) {
               assert(cur.size() == def.size());
           }
                       
        virtual bool isDefault() {
            // Element-wise comparison
            bool dec=true;
            for (unsigned int i=0;i<this->getDefault().size();i++) {
                if (this->getDefault()[i] != this->get()[i]) dec=false;
            }
            return dec;
        };
                       
        void print() {
            std::cout << "Vector '" << this->name() ;
            if (this->isDefault()) {
                std::cout << " at default value [ ";
                for (unsigned int i=0;i<this->get().size();i++) {std::cout << this->get()[i] << " ";}
                std::cout << "]" <<std::endl;
            }
            else {
                std::cout << " DIFFERS from default [ ";
                for (unsigned int i=0;i<this->get().size();i++) {
                    if (this->get()[i] == this->getDefault()[i]) {std::cout << this->get()[i] << " ";}
                    else  {std::cout << this->get()[i] << "(" << this->getDefault()[i] << ") ";}
                }
                std::cout << "]" <<std::endl;
            }
        }
};
#endif
