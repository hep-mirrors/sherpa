#ifndef InterfaceTools_h
#define InterfaceTools_h

#include <set> 
#include <iterator> 
#include <iostream> 


#include "XS_Base.H"

#include "Tree.H"
#include "Knot.H"
#include "Particle_List.H"
#include "Particle.H"
#include "Flow.H"
#include "Blob.H"


namespace SHERPA {
  class Interface_Tools {
    APACIC::Tree ** p_initrees; 
    APACIC::Tree  * p_fintree; 

    bool IsColourConnected(ATOOLS::Particle *,ATOOLS::Particle *) ;
  public :
    Interface_Tools(APACIC::Tree **,APACIC::Tree *); 
    ~Interface_Tools(); 

    void InitializeIncoming(ATOOLS::Blob *,double,double,double,double,double,double);
    void InitializeOutGoing(ATOOLS::Blob *,double,double,double,double);

    double ColourAngle(ATOOLS::Particle *,ATOOLS::Blob *);
  };

  /*!
    \class Interface_Tools
    \brief
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
}
#endif




