#ifndef DiHiggs_Virtual_H
#define DiHiggs_Virtual_H

#include <Python.h>
#include "PHASIC++/Process/Virtual_ME2_Base.H"

extern "C" {
  extern PyObject* grid_initialize(const char* grid_name);
  extern double grid_virt(PyObject* grid, double s, double t);
  extern void python_printinfo();
}

namespace MODEL {
  class Running_AlphaS;
}
 
namespace DiHiggs {

  class DiHiggs_Virtual : public PHASIC::Virtual_ME2_Base {

  public:

    DiHiggs_Virtual(const PHASIC::Process_Info& pi,
		    const ATOOLS::Flavour_Vector& flavs,
		    const std::string& grid_path);

    ~DiHiggs_Virtual();
  
    bool SetColours(const ATOOLS::Vec4D_Vector& momenta) { return true; }

    double Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& p);
  
    void Calc(const ATOOLS::Vec4D_Vector& momenta,
	      const double& born);

    void Calc(const ATOOLS::Vec4D_Vector& momenta);

  private:

    double m_symfac;

    PyObject* p_grid;

    MODEL::Running_AlphaS* p_as;
    
    double m_beta0;

  };

}

#endif
