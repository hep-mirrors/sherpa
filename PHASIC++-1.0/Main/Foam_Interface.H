#ifndef Foam_Interface_H
#define Foam_Interface_H

#include "TFOAM_INTEGRAND.h"
#include "TPSEMAR.h"
#include "TFOAM.h"
#include <vector>

namespace PHASIC {

  class Phase_Space_Handler;

  class Foam_Interface: public TFOAM_INTEGRAND {
  private:

    TPSEMAR *p_ran;
    TFOAM   *p_foam;

    Phase_Space_Handler *p_pshandler;

    std::string m_key;
    bool        m_initialize;

    std::vector<double> m_values;
    double              m_integral;

  public:

    // constructor 
    Foam_Interface(Phase_Space_Handler *const pshandler,
		   const std::string &key,const size_t dim);

    // destructor
    ~Foam_Interface();

    // member functions
    bool Initialize();

    double MCEvent();
    double Density(int ndim,double *val);    

    void SetValues(double *val);

    // inline functions
    inline TFOAM *const Foam() const { return p_foam; }

    inline const std::string &Key() const { return m_key; }

    inline const std::vector<double> &Values() const { return m_values; }

  };// end of class Foam_Interface

}// end of namespace PHASIC

#endif
