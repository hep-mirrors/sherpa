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

      double RealCorrectionsOrder(int);
      double VirtualCorrectionsOrder(int);

      double Dmod(unsigned int, unsigned int, unsigned int);
      double Smod(unsigned int);
      double Kallen(double, double, double);

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

  /*!
    \var bool Weight_Higher_Order_Corrections::m_ME
    \brief specifies whether there is a generic or specific ME to be used
  */

  /*!
    \var ME_List::code Weight_Higher_Order_Corrections::m_which
    \brief specifies which generic or specific ME will be used
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
    \fn double Weight_Higher_Order_Corrections::RealCorrectionsOrder(int)
    \brief calculates the real corrections of \f$ \mathcal{O}(\alpha^n) \f$ in the collinear approximation

    The argument to be passed is the order of \f$ \alpha \f$ at which 
    the corrections should be calculated. Lower orders than the one 
    specified will not be included. The value returned for 
    \f$ \mathcal{O}(\alpha) \f$ is
    \f[
      \frac{1}{\tilde{\beta}_0^0}\sum_{i=1}^{n_\gamma}
        \frac{\tilde{\beta}_1^1(k_i)}{\tilde{S}(k_i)}
    \f]
    with
    \f[
      \tilde{\beta}_1^1(k) 
      = -\frac{\alpha}{4\pi^2}\tilde{\beta}_0^0
         \sum_{i<j}Z_iZ_j\theta_i\theta_j
                    (\bar{\mathcal{D}}_{ij}+\bar{\mathcal{D}}_{ji})
    \f]
  */

  /*!
    \fn double Weight_Higher_Order_Corrections::VirtualCorrectionsOrder(int)
    \brief calculates the virtual corrections of \f$ \mathcal{O}(\alpha^n) \f$ in the leading log approximation

    The argument to be passed is the order of \f$ \alpha \f$ at which 
    the corrections should be calculated. Lower orders than the one 
    specified will not be included. For \f$ \mathcal{O}(\alpha) \f$ these 
    terms will only be calculated for \f$ Z^0 \f$ and \f$ W^\pm \f$ decays 
    into fermion anti-fermion pairs since in all other cases the virtual 
    contribution is negligable. The value returned is finite if the charged 
    fermions are massive.
  */

  /*!
    \fn PHOTONS_ME_Base * Weight_Higher_Order_Corrections::MESelector(Particle_Vector_Vector,Particle_Vector_Vector)
    \brief returns an instatiation of the class of the ME specified by m_which
  */

  /*!
    \fn double Weight_Higher_Order_Corrections::Dmod(unsigned int, unsigned int, unsigned int)
    \brief calculates the infrared subtracted massive dipole splitting function \f$ \bar{\mathcal{D}}_{ij}(k) \f$

    The arguments to be passed are the multipole indices \f$ i \f$ 
    and \f$ j \f$ and the bremsstrahlung photon index \f$ k \f$ as 
    they label the particles in their respective Particle_Vector.

    The massive dipole splitting functions are depending on the 
    Dipole_Type of the subdipole \f$ (ij) \f$ of the multipole and 
    the spin of the emitting particle with label \f$ i \f$.

    At the moment only the infrared subtracted splitting functions 
    for massive fermions and spin-0 bosons are implemented. The ones 
    for massive spin-1 and spin-2 bosons and spin-3/2 fermions are 
    still missing.
  */

  /*!
    \fn double Weight_Higher_Order_Corrections::Smod(unsigned int)
    \brief calculates \f$ \tilde{S}(k) \f$

    The argument passed is the index of the photon in its Particle_Vector. 
    Returns 
    \f[ 
      \tilde{S}(k) 
      = \sum\limits_{i<j}\tilde{S}_{ij}(k) 
      = \sum\limits_{i<j}Z_iZ_j\theta_i\theta_j
          \left(\frac{p_1}{(p_1\cdot k_i)}-\frac{p_2}{(p_2\cdot k_i)}\right)^2 
    \f]
  */

  /*!
    \fn void Weight_Higher_Order_Corrections::CalculateWeight()
    \brief calculates the weight of the event in the 

    The value calculated is \f$ W_{\mathcal{C}} = 1 + 
    \mathcal{O}(\alpha)_V + \mathcal{O}(\alpha)_R \f$
  */

  /*!
    \fn void Weight_Higher_Order_Corrections::CalculateMax()
    \brief calculates the maximum weight of the configuration given

    The value returned is \f$ W_{\mathcal{C}} = 1 + \mathcal{O}(\alpha)_V \f$, 
    which is the value of \f$ K=0 \f$. It is not the maximum of 
    \f$ W_{\mathcal{C}} \f$, but its contribution to the maximum total weight.
  */

  /*!
    \fn void Weight_Higher_Order_Corrections::CalculateWeightAndMax(Particle_Vector_Vector, Particle_Vector_Vector)
    \brief calculates the weight and its contribution to the maximum weight for ME-corrections

    Takes the <tt>Particle_Vector_Vector</tt>s of the dressed and undressed 
    blobs to pass on to <tt>Weight_Higher_Order_Corrections::ME_Selector<tt> 
    and calculates the weight and its contribution to the maximum.
  */
    
}

#endif
