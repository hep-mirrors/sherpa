#ifndef COMIX_Cluster_Color_Setter_H
#define COMIX_Cluster_Color_Setter_H

#include "PHASIC++/Main/Color_Integrator.H"

namespace COMIX {

  class Cluster_Algorithm;
  class Single_Process;

  class Color_Setter {
  private:

    Cluster_Algorithm *p_ca;

    Single_Process *p_xs;

    int m_cmode;

    bool SetRandomColors(const PHASIC::Int_Vector &ci,
			 const PHASIC::Int_Vector &cj);
    bool SetLargeNCColors(const PHASIC::Int_Vector &ci,
			  const PHASIC::Int_Vector &cj);

  public:

    // constructor
    Color_Setter(Cluster_Algorithm *const ca);

    // member functions
    bool SetColors(Single_Process *const xs);

  };//end of class Color_Setter

}// end of namespace COMIX

#endif
