#ifndef Interface_Tools_H
#define Interface_Tools_H

#include "Particle_List.H"
#include "Blob.H"
#include "Tree.H"

namespace SHERPA {

  class Interface_Tools {

    APACIC::Tree **p_initrees; 
    APACIC::Tree  *p_fintree; 

    ATOOLS::Particle *p_dummy;

    bool IsColourConnected(ATOOLS::Particle *,ATOOLS::Particle *) ;

  public :

    // constructor
    Interface_Tools(APACIC::Tree **,APACIC::Tree *); 

    // destructor
    ~Interface_Tools(); 

    // member functions
    void InitializeIncoming(ATOOLS::Blob *,double,double,double,double,double,double);
    void InitializeOutGoing(ATOOLS::Blob *,double,double,double,double);

    double ColourAngle(ATOOLS::Particle *,ATOOLS::Blob *);

  };// end of class Interface_Tools

  /*!
    \file 
    \brief Contains the class SHERPA::Interface_Tools
  */

  /*!
    \class Interface_Tools
    \brief To interface simple $2\to 24 processes with the parton shower.

    This class contains the core methods to interface simple matrix elements with 
    the parton shower. It is employed by the SimpleXS_Apacic_Interface.
  */
  /*!
    \var APACIC::Tree ** Interface_Tools::p_initrees 
    A pointer to Apacic's two initial state trees that are to be filled.
  */
  /*!
    \var APACIC::Tree * Interface_Tools::p_fintree 
    A pointer to Apacic's final state tree that is to be filled.
  */
  /*!
    \fn bool Interface_Tools::IsColourConnected(ATOOLS::Particle *,ATOOLS::Particle *) 
    Check's the colour flow objects of both partons and returns .true. if some
    elements coincide, i.e. if they are colour connected.
  */
  /*!
    \fn Interface_Tools::Interface_Tools(APACIC::Tree **,APACIC::Tree *) 
    In the constructor of this class the pointers to the two trees are copied.
  */
  /*!
    \fn Interface_Tools::~Interface_Tools() 
    Nothing to be done in the destructor.
  */
  /*!
    \fn void Interface_Tools::InitializeIncoming(ATOOLS::Blob *,double,double,double)
    Given the hard scattering blob, the two initial state partons are extracted and
    put into newly initialized knots of the two respective initial state trees.
    The starting scale for the initial state parton shower is given globally into both
    knots, the respective angles determined by the angle to the colour connected partner
    is specific for each individual blob and also passed as an argument. Apart from starting
    scale and angle, the Bjorken-\f$x\f$ for both partons are constructed as well as
    the maximal \f$p_\perp\f$ given by the scale of the hard process. 
  */
  /*!
    \fn void Interface_Tools::InitializeOutGoing(ATOOLS::Blob *,double,double,double)
    Given the hard process blob, the two final state partons are extracted and put into
    newly initialized knots of the final state tree. Before that, a dummy knot has been
    constructed, which is basically colour less and has no flavour, but carries some 
    information concerning the starting scale. The new knots are given some angles
    and a mass, at the moment equal to their current mass, at which they are to be found after 
    the shower ends. Finally, mother-daughter relations are established for them with the dummy 
    knot serving as their mother.
  */
  /*!
    \fn double Interface_Tools::ColourAngle(ATOOLS::Particle *,ATOOLS::Blob *)
    Calculates the angle of a specific parton with respect to its colour connected partner.
    Basically, all partons in the blob are considered, if two partons are colour connected
    to the initial one, which might be the case for gluons, the larger angle is returned.
  */

}// end of namespace SHERPA

#endif




