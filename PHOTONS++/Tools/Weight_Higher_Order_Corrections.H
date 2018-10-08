#ifndef PHOTONS_Tools_Weight_Higher_Order_Corrections_H
#define PHOTONS_Tools_Weight_Higher_Order_Corrections_H

#include "PHOTONS++/Tools/Weight_Base.H"
#include "PHOTONS++/Main/Dipole_Type.H"


namespace PHOTONS {
  class PHOTONS_ME_Base;

  class Weight_Higher_Order_Corrections: public Weight_Base {
    private:
      ATOOLS::Particle_Vector m_olddipole;
      ATOOLS::Particle_Vector m_newdipole;
      ATOOLS::Particle_Vector m_softphotons;
      double          m_M;
      unsigned int    m_n;

      PHOTONS_ME_Base * p_pme;
      PHOTONS_ME_Base * p_app;

      virtual void CalculateWeight();
      virtual void CalculateMax();

      void CalculateWeightAndMaxWithME();

    public:
      Weight_Higher_Order_Corrections(const Particle_Vector_Vector&,
                                      const Particle_Vector_Vector&,
                                      Dipole_Type::code);
      virtual ~Weight_Higher_Order_Corrections();
  };


  

  /*!
    \file Weight_Higher_Order_Corrections.H
    \brief contains the class Weight_Higher_Order_Corrections
  */

  /*!
    \class Weight_Higher_Order_Corrections
    \brief calculates \f$ W_{\mathcal{C}} \f$ of the event given
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of the member variables for Weight_Higher_Order_Corrections
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \var Particle_Vector Weight_Higher_Order_Corrections::m_olddipole
    \brief contains the multipole before treatment
  */

  /*!
    \var Particle_Vector Weight_Higher_Order_Corrections::m_newdipole
    \brief contains the multipole after treatment
  */

  /*!
    \var Particle_Vector Weight_Higher_Order_Corrections::m_softphotons
    \brief contains all photons generated
  */

  /*!
    \var double Weight_Higher_Order_Corrections::m_M
    \brief contains the mass of the decaying particle
  */

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of the member methods for Weight_Higher_Order_Corrections
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \fn Weight_Higher_Order_Corrections::Weight_Higher_Order_Corrections(Particle_Vector, Particle_Vector, Particle_Vector, double, Dipole_Type::code)
    \brief initialises all variables, looks for MEs and calculates the weight

    The arguments to be passed are the multipole after treatment, the 
    multipole before treatment, the generated photons, the mass of the 
    decaying particle and the Dipole_Type.

    Calls <tt>Weight_Higher_Order_Corrections::ME_Selector</tt> to find
    out whether there is a generic or specific ME for the process. Only
    enabled if the USE_ME switch in Fragmentation.dat is set to 1. If 
    there are no MEs available, the internal dipole approximation will 
    be called.

    Calls <tt>Weight_Higher_Order_Corrections::CalculateWeightAndMax()</tt> 
    in case the ME-correction is envoked, <tt>Weight_Higher_Order_Corrections::CalculateWeight()</tt> 
    and <tt>Weight_Higher_Order_Corrections::CalculateMax()</tt> in 
    otherwise to calculate the weight and its maximum. The latter 
    value is the contribution to the maximum total weight but not 
    the maximum of \f$ W_{\mathcal{C}} \f$ by itself.
  */

  /*!
    \fn void Weight_Higher_Order_Corrections::CalculateWeight()
    \brief calculates the weight of the event in the 

    This function is not to be used in this class
  */

  /*!
    \fn void Weight_Higher_Order_Corrections::CalculateMax()
    \brief calculates the maximum weight of the configuration given

    The value returned is \f$ W_{\mathcal{C}} = 1 + \mathcal{O}(\alpha)_V \f$, 
    which is the value of \f$ K=0 \f$. It is not the maximum of 
    \f$ W_{\mathcal{C}} \f$, but its contribution to the maximum total weight.
  */

    
}

#endif
