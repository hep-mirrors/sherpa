#ifndef Knot_H
#define Knot_H

#include "MathTools.H"
#include "Parton.H"

namespace APACIC {
  // defines one parton splitting (basic element of a Tree)
  class Knot {
  public :
    //! constructor
    Knot():left(0),right(0),prev(0),part(0),stat(0),
           t(0.),tout(0.),E2(0.),costh(0.),phi(0.),z(0.),thcrit(M_PI),
	   maxpt2(1.e10),x(0.),kn_no(-1) {};
    //! copy constructor
    Knot(Knot * in);
    //! copy only the interior (not the links)
    void CopyData(Knot * in);
    //! copy the knot properties
    void Copy(Knot * in);

    ~Knot();
    //! parton yields :  flavour, momentum, spacepoint, lifetime
    APHYTOOLS::Parton * part;

    //! decay: virtuality, end virtuality, energy squared, angles, decay fraction 
    double t, tout, E2, costh, phi, z, thcrit;
    //! max p_T so far
    double maxpt2;
    //! x energy part of beam particle!
    double x;
    //! decay already fixed (=0) or not (=1)
    int stat;

    // tree structure:
    Knot* left;
    Knot* right;   
    Knot* prev;

    // reference number
    int kn_no;
  };

  //----------------------------------------------------------------------
  //    inline functions
  //----------------------------------------------------------------------
  inline Knot::Knot(Knot * ink) {
    prev   = ink->prev;
    left   = ink->left;
    right  = ink->right;

    APHYTOOLS::Parton * newp;
    newp   = new APHYTOOLS::Parton(ink->part);
    part   = newp;
    stat   = ink->stat;

    t      = ink->t;
    tout   = ink->tout;
    E2     = ink->E2;
    z      = ink->z;
    x      = ink->x;
    maxpt2 = ink->maxpt2;
    costh  = ink->costh;
    phi    = ink->phi;
    thcrit = ink->thcrit;
    
    kn_no  = ink->kn_no;
  };

  inline void Knot::CopyData(Knot * ink) {
    part->Copy(ink->part);
    stat   = ink->stat;

    t      = ink->t;
    tout   = ink->tout;
    E2     = ink->E2;
    z      = ink->z;
    x      = ink->x;
    maxpt2 = ink->maxpt2;
    costh  = ink->costh;
    phi    = ink->phi;
    thcrit = ink->thcrit;
    
    kn_no  = ink->kn_no;
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
}

#endif





