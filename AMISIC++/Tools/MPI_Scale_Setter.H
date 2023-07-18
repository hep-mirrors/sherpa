#ifndef AMISIC_Tools_MPI_Scale_Setter_H
#define AMISIC_Tools_MPI_Scale_Setter_H

#include "PHASIC++/Scales/Scale_Setter_Base.H"

namespace AMISIC {
  class MPI_Scale_Setter: public PHASIC::Scale_Setter_Base {
  public:
    MPI_Scale_Setter(const PHASIC::Scale_Setter_Arguments &args);

    double Calculate(const std::vector<ATOOLS::Vec4D> &p,
		     const size_t &mode) final;
  };
}

#endif
