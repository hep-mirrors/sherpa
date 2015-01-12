#ifndef SHERPAMATRIX
#define SHERPAMATRIX
#include <iostream>
#include <iomanip>
#include <typeinfo>
#include <string>
#include <vector>
#include <cmath>
#include "ATOOLS/Org/ParameterBase.H"
#include "ATOOLS/Org/boost/assign/std/vector.hpp" // for 'operator+=()'

using namespace boost::assign; // bring 'operator+=()' into scope

// Matrices
template <typename T>
class SherpaMatrix : public ParameterBase< std::vector< std::vector <T> > > {
    public:
        SherpaMatrix(std::string name, const std::vector< std::vector <T> >& cur, const std::vector< std::vector <T> >& def, const bool allowChangeAfterLock=false)
           : ParameterBase<std::vector< std::vector <T> > >(name, cur, def,allowChangeAfterLock) {
           }
                       
        virtual bool isDefault() {
            // Element-wise comparison
            bool dec=true;
            for (unsigned int i=0;i<this->getDefault().size();i++) {
                for (unsigned int j=0;j<this->getDefault()[i].size();j++) {
                    if (this->getDefault()[i][j] != this->get()[i][j]) dec=false;
                }
            }
            return dec;
        };
        void print() {
            std::cout << "Matrix '" << this->name() ;
            if (this->isDefault()) {
                std::cout << " at default value" << std::endl <<"[" << std::endl;
                for (unsigned int i=0;i<this->get().size();i++) {
                    for (unsigned int j=0;j<this->get()[i].size();j++) {
                    std::cout << this->get()[i][j] << " ";
                    }
                    std::cout << std::endl;
                }
                std::cout << "]" <<std::endl;
            }
            else {
                std::cout << " DIFFERS from default" << std::endl <<"[" <<std::endl;
                for (unsigned int i=0;i<this->get().size();i++) {
                    for (unsigned int j=0;j<this->get()[i].size();j++) {
                        if (this->get()[i][j] == this->getDefault()[i][j]) {std::cout << this->get()[i][j] << " ";}
                        else {std::cout << this->get()[i][j] << "(" << this->getDefault()[i][j] <<") ";}
                    }
                    std::cout << std::endl;
                }
                std::cout << "]" <<std::endl;
            }
        }
};
#endif
