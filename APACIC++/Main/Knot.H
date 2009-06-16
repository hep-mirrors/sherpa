#ifndef APACIC_Main_Knot_H
#define APACIC_Main_Knot_H

#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Phys/Simple_Polarisation_Info.H"

namespace APACIC {

  class Knot {
  public: 

    Knot *prev, *left, *right, *decay;

    ATOOLS::Particle *part;
    ATOOLS::Vec4D     lp[2], cms;

    int    stat, kn_no, kn_id, shower, dir, oc[2];
    bool   didkin;
    double t, tout, tmax, z, zs, E2, costh, phi, thcrit, maxpt2, x, pt2lcm;
    double smaxpt2, sthcrit, minpt2, qjv, tmo, lz[2], lE2[2], asme;
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

  inline size_t IdCount(const size_t &id)
  {
    size_t ic(id), cn(0);
    for (size_t i(0);ic>0;++i) {
      size_t c(1<<i);
      if (ic&c) { ++cn;	ic-=c; }
    }
    return cn;
  }

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





