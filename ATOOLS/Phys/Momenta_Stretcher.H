#ifndef ATOOLS_Phys_Momenta_Stretcher_H
#define ATOOLS_Phys_Momenta_Stretcher_H

#include <vector>
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Math/Vector.H"

namespace ATOOLS {
  
  class Momenta_Stretcher {
  public:
    Momenta_Stretcher() {}
    bool MassThem(const int, const int, Vec4D *, const double *,const double=-1.e-14);
    bool MassThem(const int n0, std::vector<Vec4D>& moms, std::vector<double> masses,const double=-1.e-14);
    bool ZeroThem(const int, const int, Vec4D *,const double=-1.e-14);
    bool ZeroThem(const int n0, std::vector<Vec4D>& moms,const double=1.e-14);
    bool StretchBlob( Blob* blob, std::vector<double> masses = std::vector<double>());
    bool StretchMomenta( const Particle_Vector& outparts, std::vector<Vec4D>& moms);
    bool StretchMomenta( const Particle_Vector& outparts, std::vector<double>& masses);
  };
}
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member methods for Momenta_Stretcher.
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \fn bool Momenta_Stretcher::MassThem(const int,const int,ATOOLS::Vec4D *,const double *)
    \brief Puts the momenta on new mass-shells, respecting four-momentum
    conservation of the overall ensemble and the original directions as far as
    possible.

    If this succeeds, "true" is returned; if, in contrast, the sum of all new
    masses is larger than the overall energy available, such a transformation
    cannot work out and the method yields "false". Momenta have to be Zero'd
    first!

    <var>n0</var> specifies the first array/vector index to stretch, all
    momenta below that index will be ignored when stretching.

    In general the momenta set need to be in its overall c.m. frame for this method to work.
    Then the corresponding algorithm of Rambo is employed. It shifts the masses and uses
    scale transformations in an iterative procedure until convergence is reached. 
    
    For two momenta, the system does not need to be in its c.m. frame. There, the system is boosted 
    into its c.m. frame, the overall energy is calculated and the three-momenta are adjusted, 
    respecting the original directions. Finally, the systemv is boosted back.
  */
  /*!
    \fn bool Momenta_Stretcher::ZeroThem(const int,ATOOLS::Vec4D *)
    \brief Puts the momenta on zero mass, respecting four-momentum conservation of the overall
    ensemble and the original directions as far as possible. Momenta have to be in CMS!

    For two momenta this is straight-forwardly done after going to the two vectors c.m. system
    and calculating the total energy. For an arbitrary number of momenta, the inverse of the
    Rambo algorithm used in Momenta_Stretcher::MassThem(const int,ATOOLS::Vec4D *,const double *) is
    employed.
  */
  /*!
    \fn bool Momenta_Stretcher::StretchBlob(Blob*, std::vector<double>)
    \brief Puts the momenta of the outgoing particles on new mass shells either
           given by a vector of doubles or (if not given) to their FinalMass()es

    Uses ZeroThem and MassThem to accomplish the stretching, and also cares for
    boosting the momenta into CMS in the beginning and back afterwards.
  */
  /*!
    \fn bool Momenta_Stretcher::StretchMomenta(
             const Particle_Vector&,std::vector<Vec4D>& )
    \brief Puts the given momenta on new mass shells given by the FinalMasses
           of the Particle_Vector.

    Uses ZeroThem and MassThem to accomplish the stretching, and also cares for
    boosting the momenta into CMS in the beginning and back afterwards.
  */

#endif
