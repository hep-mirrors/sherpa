#ifndef PDF_Main_ISR_Base_H
#define PDF_Main_ISR_Base_H

#include "PDF/Main/PDF_Base.H"

namespace PDF {

  class ISR_Base { 
  protected:

    PDF_Base *p_pdf;

    std::string     m_type;
    ATOOLS::Flavour m_bunch;

    double m_weight, m_exponent, m_xmax;
    bool   m_on, m_ewon;

  public:

    // constructor
    ISR_Base(PDF_Base *pdf);

    // destructor
    virtual ~ISR_Base();

    // member functions
    bool On() const   { return m_on;   }
    bool EWOn() const { return m_ewon; }

    void AssignKeys(ATOOLS::Integration_Info *const info) { p_pdf->AssignKeys(info); }

    // inline functions
    inline const std::string     Type() const    { return m_type;  }
    inline const ATOOLS::Flavour Flavour() const { return m_bunch; }

    inline bool Collinear(const double kp2) const { return p_pdf->Collinear(kp2); }

    inline double Cut(const std::string &type) const { return p_pdf->Cut(type); }

    inline double Exponent() const { return m_exponent; }
    inline double XMax() const     { return m_xmax;     } 

    inline PDF_Base *PDF() { return p_pdf; }

    inline void Extract(const ATOOLS::Flavour flavour,const double x) const 
    { p_pdf->Extract(flavour,x); }
    inline void Reset() const                 
    { p_pdf->Reset();    }

    // member functions
    virtual bool   CalculateWeight(double x,double z,double kp2,double Q2) = 0;
    virtual double Weight(ATOOLS::Flavour) = 0;

  };// end of class ISR_Base

  /*!
    \class ISR_Base
    \brief A mainly virtual class for different initial state radiation treatments
           and a wrapper around a corresponding PDF.
    
    This class is the base class for all possible treatments of initial state radiation
    (ISR). At the moment, two treatments are implemented, namely first of all the 
    case of no partonic substructure of the incoming bunch, realized through Intact,
    and, second, the case of partonic substructure as described by PDFs or
    structure functions, realized through Structure_Function. Therefore this class
    contains some information about the incoming bunch and in its realization
    as Structure_Function also concerning the PDF related to it.

    \todo It has to be sorted out what happens to the YFS framework.
  */
  /*!
    \var std::string ISR_Base::m_type
    The type of the ISR_Base, either "(None)" or "(SF)" for partonic bunches or
    a handling through PDFs/structure functions.

    \sa ATOOLS::ISR_Type
  */
  /*!
    \var ATOOLS::Flavour ISR_Base::m_bunch.
    The flavour of the incoming bunch particle.
  */
  /*!
    \var double ISR_Base::m_weight
    The internal weight for the spectrum. It is calculated through the method
    CalculateWeight(const double,const double).
  */
  /*!
    \fn bool ISR_Base::CalculateWeight(const double,const double)
    This is the method that will calculate the weights for the parton densities in
    dependence on the energy fraction of the outgoing particle w.r.t. the incoming one
    and in dependence on the scale of the process. 
  */
  /*!
    \fn double ISR_Base::Weight(ATOOLS::Flavour)
    Returns the weight after the flavours are checked.
  */
  /*!
    \fn bool ISR_Base::On()
    Returns .false. if no beam-handling is necessary, i.e. for monochromatic beams and
    .true. in all other cases.
  */
  /*!
    \fn double ISR_Base::Exponent() 
    This is a helper for the phase space integration. It returns a characteristic value -
    usually an exponent - for simple channel structures. 
  */
  /*!
    \fn double ISR_Base::Xmax()
    The maximally accessible energy fraction an outgoing particle can have.
  */
  /*!
    \fn PDF_Base * ISR_Base::PDF()
    Returns the internal pointer to the constituting PDF.
  */

}// end of namespace PDF

#endif
