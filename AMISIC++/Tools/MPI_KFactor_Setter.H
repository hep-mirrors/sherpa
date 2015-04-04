#ifndef AMISIC__Tools__MPI_KFactor_Setter_H
#define AMISIC__Tools__MPI_KFactor_Setter_H

#include "PHASIC++/Scales/KFactor_Setter_Base.H"

namespace AMISIC {

  class MPI_KFactor_Setter: public PHASIC::KFactor_Setter_Base {
  private:

    static double s_pt0;

  public:

    MPI_KFactor_Setter(const PHASIC::KFactor_Setter_Arguments &args);

    double KFactor();

    inline static void SetPT0(const double &pt0) { s_pt0=pt0; }

  };// end of class MPI_KFactor_Setter

}// end of namespace AMISIC

#endif
