#ifndef Perturbative_Interface_h
#define Perturbative_Interface_h

#include "Matrix_Element_Handler.H"
#include "Shower_Handler.H"
#include "Blob.H"
#include "Blob_List.H"

namespace SHERPA {
  class Perturbative_Interface {
  protected :
    Matrix_Element_Handler * p_me;
    Shower_Handler         * p_shower;        
    bool                     m_ini,m_fin;
    int                      m_maxjetnumber;
    double                   m_weight;
    ATOOLS::Flavour        * p_fl;
    ATOOLS::Vec4D          * p_moms;

    int                      m_nin,m_nout;
  public :
    Perturbative_Interface(Matrix_Element_Handler * _me,Shower_Handler * _shower) :
      p_me(_me), p_shower(_shower), 
      m_ini(_shower->ISROn()), m_fin(_shower->FSROn()), m_maxjetnumber(_me->MaxJets()), m_weight(1.),
      p_fl(NULL), p_moms(NULL) {}

    virtual ~Perturbative_Interface() 
    {
      if (p_fl)   { delete p_fl;   p_fl   = NULL; }
      if (p_moms) { delete p_moms; p_moms = NULL; }
    }
    virtual bool   DefineInitialConditions(ATOOLS::Blob *) = 0;
    virtual bool   FillBlobs(ATOOLS::Blob_List *) = 0; 
    double GetWeight() { return m_weight; }
  };

  /*!
    \class Perturbative_Interface
    \brief The base class for all interfaces of matrix elements and parton showers.

  */
  /*!
    \var Matrix_Element_Handler * Perturbative_Interface::p_me
    A pointer to the Matrix_Element_Handler. It ensures that the information from there,
    like flavours and momenta of the partons involved as well as Feynman diagrams for the
    eventual reconstruction of parton shower histories is passed into the interface.
  */
  /*!
    \var Shower_Handler * Perturbative_Interface::p_shower        
    A pointer to the Shower_Handler to fill in the trees along which the parton shower
    evolves.
  */
  /*!
    \var bool Perturbative_Interface::m_ini
    A flag whether the initial state shower is enabled.
  */
  /*!
    \var bool Perturbative_Interface::m_fin
    A flag whether the final state shower is enabled.
  */
  /*!
    \var int Perturbative_Interface::m_maxjetnumber
    The maximal number of jets stemming from the matrix elements. This is important to
    start the parton shower properly and to steer the application of jet vetos inside the
    parton shower and the corresponding Sudakov weight inside the interface of matrix elements
    and the shower.
  */
  /*!
    \var double Perturbative_Interface::m_weight
    The weight for the combination procedure, mainly a product of ratios of \f$\alpha_s\f$ at
    different scales and some Sudakov factors.
  */
  /*!
    \var ATOOLS::Flavour * Perturbative_Interface::p_fl
    A pointer to the flavours of the partons coming from the matrix element.
  */
  /*!
    \var ATOOLS::Vec4D * Perturbative_Interface::p_moms
    A pointer to the momenta of the partons coming from the matrix element.
  */
  /*!
    \var int Perturbative_Interface::m_nin
    The number of incoming partons from the matrix element.
  */
  /*!
    \var int Perturbative_Interface::m_nout
    The number of outgoing partons from the matrix element.
  */
  /*!
    \fn Perturbative_Interface::Perturbative_Interface(Matrix_Element_Handler * _me,Shower_Handler * _shower)
    The constructor mainly copies the pointers to the essential objects.
  */
  /*!
    \fn virtual Perturbative_Interface::~Perturbative_Interface() 
    Deletes the poniters to flavours and momenta, p_fl and p_moms.
  */
  /*!
    \fn virtual bool Perturbative_Interface::Perturbative_Interface::DefineInitialConditions(ATOOLS::Blob *)
    This is the central method of all interfaces. After the (hard signal) blob is filled by the
    matrix element handler it is passed into the corresponding interface which in turn reconstructs
    a parton shower history and - by doing so - defines the start conditions for hte shower.
  */
  /*!
    \fn virtual bool Perturbative_Interface::FillBlobs(ATOOLS::Blob_List *)
  */
  /*!
    \fn double Perturbative_Interface::GetWeight()
    Returns the weight connected to the merging procedure.
  */
}
#endif







