#ifndef SimpleXS_Apacic_Interface_h
#define SimpleXS_Apacic_Interface_h

#include "Perturbative_Interface.H"
#include "Interface_Tools.H"

namespace SHERPA {
  class SimpleXS_Apacic_Interface : public Perturbative_Interface {
    Interface_Tools * p_tools;

    bool InitColours(ATOOLS::Blob *);
  public :
    SimpleXS_Apacic_Interface(Matrix_Element_Handler *,Shower_Handler *);
    ~SimpleXS_Apacic_Interface();
    bool   DefineInitialConditions(ATOOLS::Blob *);
    bool   FillBlobs(ATOOLS::Blob_List *); 
  };
  /*!
    \class SimpleXS_Apacic_Interface
    \brief The explicit realization of the Perturbative_Interface for the case of the internal
           \f$2\to 2\f$ processes and Apacic.

    This class provides the interface betwenn the internal explicit matrix elements for \f$2\to 2\f$
    processes with the parton shower of Apacic. It relies heavily on the Interface_Tools to calculate
    suitable starting conditions for the parton shower (like, e.g. the angle to the colour connected
    partner) and to fill in the newly instantiated partons into the binary tree structure used
    within the parton shower.
  */
  /*!
    \var Interface_Tools * SimpleXS_Apacic_Interface::p_tools
    A pointer to the Interface_Tools that are used throughout this class.
  */
  /*!
    \fn bool SimpleXS_Apacic_Interface::InitColours(ATOOLS::Blob *)
    This is the method that does the real work. It copies the colours from the EXTRA_XS::XS_Base
    into the flow-objects of the aprtons, finds the colour connected partners and calculates
    the initial opening angles for each parton and then initializes the binary trees for both the
    initial and final state parton shower through corresponding methods of the Interface_Tools.
  */
  /*!
    \fn SimpleXS_Apacic_Interface::SimpleXS_Apacic_Interface(Matrix_Element_Handler *,Shower_Handler *)
    Copies the pointers to the Shower_Handler and Matrix_Element_Handler and initializes the
    Interface_Tools.
  */
  /*!
    \fn SimpleXS_Apacic_Interface::~SimpleXS_Apacic_Interface()
    Deletes the Interface_Tools.
  */
  /*!
    \fn bool SimpleXS_Apacic_Interface::DefineInitialConditions(ATOOLS::Blob *)
    Checks whether the blob is of the form of a \f$2\to 2\f$ process, if no the program will
    terminate. Otherwise the blob (and all attached partons) is boosted in its c.m. frame and
    InitColours will initialize the parton showers by filling hte trees accordingly.
  */
  /*!
    \fn bool SimpleXS_Apacic_Interface::FillBlobs(ATOOLS::Blob_List *) 
    Returns .true. without any further manipulations. 
    \todo We have to check whether we need a FillBlobs method at all.
  */
}
#endif
