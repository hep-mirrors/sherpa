#ifndef Knot_H
#define Knot_H

#include "MathTools.H"
#include "Particle.H"
#include "Simple_Polarisation_Info.H"

namespace APACIC {
  // defines one parton splitting (basic element of a Tree)
  class Knot {
  public: 
    Knot * left, * right, * prev;
    ATOOLS::Particle * part;
    int stat;
    double t, tout, tmax, E2, costh, phi, z, thcrit, maxpt2, x, pt2lcm;
    int kn_no;
    ATOOLS::Simple_Polarisation_Info polinfo;
  public :
    Knot():left(0),right(0),prev(0),part(0),stat(0),
           t(0.),tout(0.),tmax(0.),E2(0.),costh(0.),phi(0.),z(0.),thcrit(M_PI),
	   maxpt2(1.e10),x(0.),pt2lcm(0.0),kn_no(-1),
	   polinfo(ATOOLS::Simple_Polarisation_Info()) {};
    Knot(Knot * in);
    ~Knot();

    //! copy only the interior (not the links)
    void CopyData(Knot * in);
    //! copy the knot properties
    void Copy(Knot * in);
  };

  //----------------------------------------------------------------------
  //    inline functions
  //----------------------------------------------------------------------
  inline Knot::Knot(Knot * ink) {
    prev   = ink->prev;
    left   = ink->left;
    right  = ink->right;

    ATOOLS::Particle * newp;
    newp   = new ATOOLS::Particle(*ink->part);
    part   = newp;
    stat   = ink->stat;

    t      = ink->t;
    tout   = ink->tout;
    tmax   = ink->tmax;
    E2     = ink->E2;
    z      = ink->z;
    x      = ink->x;
    maxpt2 = ink->maxpt2;
    costh  = ink->costh;
    phi    = ink->phi;
    thcrit = ink->thcrit;
    
    kn_no  = ink->kn_no;
    polinfo= ink->polinfo;
  };

  inline void Knot::CopyData(Knot * ink) {
    if (part && ink->part) part->Copy(ink->part);
    stat   = ink->stat;

    t      = ink->t;
    tout   = ink->tout;
    tmax   = ink->tmax;
    E2     = ink->E2;
    z      = ink->z;
    x      = ink->x;
    maxpt2 = ink->maxpt2;
    costh  = ink->costh;
    phi    = ink->phi;
    thcrit = ink->thcrit;
    
    kn_no  = ink->kn_no;
    polinfo= ink->polinfo;
  }

  inline void Knot::Copy(Knot * ink) {
    prev   = ink->prev;
    left   = ink->left;
    right  = ink->right;

    CopyData(ink);
  }


  inline Knot::~Knot() {
    if (part) delete part;
  }

  class Tree;

  std::ostream& operator<< (std::ostream&,Knot&);

  /*!
    \file
    \brief contains the class APACIC::Knot
  */

  /*!
    \class Knot 
    \brief Defines one parton splitting (basic element of a Tree)
  */
}

#endif





