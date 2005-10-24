#ifndef Knot_H
#define Knot_H

#include "MathTools.H"
#include "Particle.H"
#include "Simple_Polarisation_Info.H"

namespace APACIC {

  class Knot {
  public: 

    Knot *prev, *left, *right;

    ATOOLS::Particle *part;

    int    stat, kn_no;
    bool   didkin;
    double t, tout, tmax, z, zs, E2, costh, phi, thcrit, maxpt2, x, pt2lcm;
    ATOOLS::Simple_Polarisation_Info polinfo;

  public :

    // constructors
    inline Knot():
      prev(NULL), left(NULL),right(NULL), part(NULL), stat(0), kn_no(-1), 
      didkin(false), t(0.0), tout(0.0), tmax(0.0), z(0.0), zs(0.0),
      E2(0.0), costh(0.0), phi(0.0), thcrit(M_PI),
      maxpt2(1.0e10), x(0.), pt2lcm(0.0) {}
    inline Knot(Knot * k):
      prev(k->prev), left(k->left), right(k->right),
      part(new ATOOLS::Particle(*k->part)), stat(k->stat), kn_no(k->kn_no), 
      didkin(k->didkin), t(k->t), tout(k->tout), tmax(k->tmax), z(k->z), zs(0.0),
      E2(k->E2), costh(k->costh), phi(k->phi), thcrit(k->thcrit), 
      maxpt2(k->maxpt2), x(k->x), pt2lcm(0.0), polinfo(k->polinfo) {}
    
    // destructor
    inline ~Knot() { if (part!=NULL) delete part; }

    // inline functions
    inline void CopyData(const Knot *const k) 
    {
      if (part!=NULL && k->part) part->Copy(k->part);
      stat=k->stat;
      kn_no=k->kn_no;
      didkin=k->didkin;
      t=k->t;
      tout=k->tout;
      tmax=k->tmax;
      E2=k->E2;
      z=k->z;
      zs=k->zs;
      x=k->x;
      maxpt2=k->maxpt2;
      costh=k->costh;
      phi=k->phi;
      thcrit=k->thcrit;
      polinfo=k->polinfo;
    }

    inline void Copy(const Knot *const k) 
    {
      prev=k->prev;
      left=k->left;
      right=k->right;
      CopyData(k);
    }

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





