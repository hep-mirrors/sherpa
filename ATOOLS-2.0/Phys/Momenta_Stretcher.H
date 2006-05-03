#ifndef Momenta_Stretcher_H
#define Momenta_Stretcher_H

namespace ATOOLS {

  class Vec4D;
  
  class Momenta_Stretcher {
    public:
      Momenta_Stretcher() {}
      bool MassThem(const int, const int, Vec4D *, const double *);
      bool ZeroThem(const int, const int, Vec4D *);
  };
}
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member methods for Momenta_Stretcher.
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \fn bool Momenta_Stretcher::MassThem(const int,const int,ATOOLS::Vec4D *,const double *)
    \brief Puts the momenta on new mass-shells, respecting four-momentum conservation of the overall
    ensemble and the original directions as far as possible. If this succeeds, "true" is returned; if, 
    in contrast, the sum of all new masses is larger than the overall energy available, such a 
    transformation cannot work out and the method yields "false".

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
    ensemble and the original directions as far as possible.

    For two momenta this is straight-forwardly done after going to the two vectors c.m. system
    and calculating the total energy. For an arbitrary number of momenta, the inverse of the
    Rambo algorithm used in Momenta_Stretcher::MassThem(const int,ATOOLS::Vec4D *,const double *) is
    employed.
  */

#endif
