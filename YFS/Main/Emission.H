#ifndef YFS_Main_Emission_H
#define YFS_Main_Emission_H

#include "YFS/Main/Photon.H"

namespace YFS {

  struct EmissionResult {
    Photon_Vector photons;
    Photon_Vector me_photons;
    double weight = 1.;
    ATOOLS::Vec4D photon_sum;

    enum class Failure { none, initialize, makefsr, masswgt, formfactor };
    Failure fail = Failure::none;
    bool ok() const { return fail == Failure::none; }
  };

}// end of namespace YFS

#endif
