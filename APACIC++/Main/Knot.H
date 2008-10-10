#ifndef Knot_H
#define Knot_H

#include "MathTools.H"
#include "Particle.H"
#include "Simple_Polarisation_Info.H"

namespace APACIC {

  class Knot {
  public: 

    Knot *prev, *left, *right, *decay;

    ATOOLS::Particle *part;
    ATOOLS::Vec4D     lp[2], cms;

    int    stat, kn_no, shower, maxjets, dir, oc[2];
    bool   didkin;
    double t, tout, tmax, z, zs, E2, costh, phi, thcrit, maxpt2, x, pt2lcm;
    double smaxpt2, sthcrit, minpt2, qjv, qljv, tmo, lz[2], lE2[2];
    ATOOLS::Simple_Polarisation_Info polinfo;

  public :

    // constructors
    Knot();
    Knot(Knot * k);
    
    // destructor
    inline ~Knot() { if (part!=NULL) delete part; }

    // inline functions
    void CopyData(const Knot *const k);

    void Copy(const Knot *const k);

    void Store(const int i=0);
    void Restore(const int i=0);

    bool CheckMomentumConservation(const bool force=false) const;

  };// end of class Knot

  std::ostream &operator<<(std::ostream&,const Knot &knot);

  /*!
    \file
    \brief contains the class APACIC::Knot
  */

  /*!
    \class Knot 
    \brief Defines one parton splitting (basic element of a Tree)
  */

}// end of namespace APACIC

#endif





