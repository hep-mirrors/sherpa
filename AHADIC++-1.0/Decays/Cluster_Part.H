#ifndef Cluster_Part_H
#define Cluster_Part_H


#include "Cluster.H"
#include "Histogram.H"
#include "Hadronisation_Parameters.H"
#include "Return_Value.H"
#include <map>


namespace AHADIC {
  class Cluster_Part {
  protected:
    int             m_4Qmode, m_ybar_mode;
    bool            m_ana;
    ATOOLS::Vec4D   m_momenta[4];

    Pair_Popper   * p_popper;
    double          m_fraction, m_ystar_expvalue, m_ystar_sigma;

    std::map<std::string,ATOOLS::Histogram *> m_histograms;

    bool   BuildKinematics(Cluster * const,ATOOLS::Flavour &);
    double GetYStar(const double,const double);
    double GetYBar(const double,const double);
  public:
    Cluster_Part();
    ~Cluster_Part();
    bool TestDecay(Cluster *);
    bool UpdateDecay(Cluster *,const int);
  };

  /*!
    \file Cluster_Part.H
    \brief
  */
  /*!
    \class Cluster_Part
    \brief
  */
  /*!
    \var int Cluster_Part::m_4Qmode
    \brief
  */
  /*!
    \var int Cluster_Part::m_ybar_mode
    \brief
  */
  /*!
    \var double Cluster_Part::m_fraction
    \brief
  */
  /*!
    \var double Cluster_Part::m_ystar_expvalue
    \brief
  */
  /*!
    \var double Cluster_Part::m_ystar_sigma
    \brief
  */
  /*!
    \var Pair_Popper * Cluster_Part::p_popper
    \brief
  */
  /*!
    \var ATOOLS::Vec4D Cluster_Part::m_momenta[4]
    \brief
  */
  /*!
    \var bool Cluster_Part::m_ana
    \brief
  */
  /*!
    \var std::map<std::string,ATOOLS::Histogram *> Cluster_Part::m_histograms
    \brief
  */
  /*!
    \fn Cluster_Part::Cluster_Part()
    \brief
  */
  /*!
    \fn Cluster_Part::~Cluster_Part()
    \brief
  */
  /*!
    \fn bool Cluster_Part::TestDecay(Cluster *)
    \brief
  */
  /*!
    \fn bool Cluster_Part::UpdateDecay(Cluster *,const int)
    \brief
  */
  /*!
    \fn bool Cluster_Part::BuildKinematics(Cluster * const,ATOOLS::Flavour &)
    \brief Builds the kinematics of a cluster->cluster+cluster decay in its own
    rest frame.

    The algorithm is as follows:
    - A flavour is selected through the AHADIC::Pair_Popper.  The upper limit
      thereby is given by \f$M_C-m_0-m_3\f$, where \f$M_C\f$ is the mass of the
      cluster and \f$m_{0,3}\f$ are the masses of its constituents.  The mass of 
      this popped new flavour will be used for \f$m_{1,2}\f$. 
    - Then the maximal momentum, \f$p_{\rm max}^2\f$ is determined through
      \f[
      p_{\rm max}^2 = \frac{[M_C^2+(m_0+m_1)^2-(m_2+m_3)^2]^2 - 
                            4(m_0+m_1)^2(m_2+m_3)^2}{4M_C^2}\,,
      \f]
      implying the minimal hadron masses to be \f$(m_0+m_1)\f$ and \f$(m_2+m_3)\f$.
      The maximal transverse momentum is given by
      \f[
      p_{\perp,{\rm max}}^2 = f p_{\rm max}^2
      \f]
      with \f$f\f$ given by Cluster_Part::m_fraction.
    - Equipped with these constraints, the basic kineamtic quantitites are calculated.
      - The actual \f$p_\perp\f$ is given through a call to 
        Pair_Popper::SelectPT(\f$p_{\perp,{\rm max}}^2\f$).
      - \f$y^*\f$ is defined through 
        Cluster_Part::GetYStar(\f$p_\perp^2+m_1^2,\,p_{\perp,{\rm max}}^2\f$)
      - The actual kinematical mass squared, \f$s_{qq}\f$, of the popped quark pair is 
        calculated and used to calculate \f$\bar y\f$ through a call to
	Cluster_Part::GetYBar(\f$p_{\rm max}^2,\,s_{qq}\f$).
      - These quantities are then used to calculate energy fractions and actual 
        energies.
      The procedure is repeated, until physical energies for the original cluster
      constituents (larger than their mass) emerge.  

  */
  /*!
    \fn bool Cluster_Part::GetYStar(const double,const double)
    \brief Returns a value for \f$y^*\f$ in dependence on \f$m_\perp^2\f$ and
    \f$m_{\rm max}^2\f$.
    
    The maximal value of \f$y^*\f$ is given by
    \f[
    y^*_{\rm max} = \frac12\log\frac{m_{\rm max}^2}{4m_\perp^2}\,,
    \f]
    and the shift \f$\Delta y^*\f$ is given by
    \f[
    \Delta y^* = \langle\delta y^*\rangle \cdot y^*_{\rm max}\,.
    \f]
    with \f$\langle\delta y^*\rangle \f$ given by Cluster_Part::m_ystar_expvalue.
    \f$y^*\f$ is then distributed according to a shifted Gaussian, i.e. according to
    \f[
    \exp\frac{-(y^*-\Delta y^*)^2}{\sigma_y}\,,
    \f]
    where \f$\sigma_y\f$ is Cluster_Part::m_ystar_sigma.
  */
  /*!
    \fn bool Cluster_Part::GetYBar(const double,const double)
    \brief Returns a value for \f$\bar y\f$ in dependence of \f$p_{\rm max}^2\f$ and
    \f$s_{qq}\f$.  

    The maximal value of \f$\bar y\f$ is given by
    \f[
    \bar y_{\rm max} = \frac12\log\frac{p_{\rm max}^2}{s_{qq}}
    \f]
    Then, in depedence of Cluster_Part::m_ybar_mode, \f$\bar y\f$ is chosen according to
    - Cluster_Part::m_ybar_mode=2 (Gaussian distribution)
      \f[
      d\bar y\sim\exp\frac{\bar y}{\bar y_{\rm max}}
      \f]
    - Cluster_Part::m_ybar_mode=1 (Hyperbolical)
      \f[
      d\cosh\bar y\sim 1
      \f]
    - Cluster_Part::m_ybar_mode=0 (Uniform)
      \f[
      d\bar y\sim 1
      \f]
  */
}
#endif
