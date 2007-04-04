#ifndef Gluon_Decayer_H
#define Gluon_Decayer_H

#include "Cluster.H"
#include "Blob.H"

namespace AHADIC {
  struct DecaySpecs { double weight,zmin,zmax; };
  typedef std::map<ATOOLS::Flavour,DecaySpecs *> FlavDecayMap;
  typedef FlavDecayMap::iterator                 FDIter;

  class Gluon_Decayer {
  private:
    FlavDecayMap        m_options;
    double              m_z,m_phi,m_mass;
    ATOOLS::Flavour     m_flav;
    ATOOLS::Vec4D       m_p1vec, m_p2vec;

    bool                Shift(Proto_Particle_List * pl);
    bool                DecayIt(Proto_Particle * part,
				Proto_Particle *& part1,Proto_Particle *& part2);
    void                SelectDecay(const ATOOLS::Vec4D &);
    double              ZWeight(const double, const double);
    void                SelectZ(const DecaySpecs*);
    void                SelectPhi();
    void                BuildKinematics(const ATOOLS::Vec4D &);
  public:
    Gluon_Decayer();
    ~Gluon_Decayer();
    bool DecayList(Proto_Particle_List * plin=NULL);
  };
  /*!
    \file Gluon_Decayer.H
    \brief Contains the class AHADIC::Gluon_Decayer and the stuct AHADIC::DecaySpecs. 
  */

  /*! 
    \struct DecaySpecs
    \brief This struct carries a flavour, a weight and a z-range for the decay of an individual gluon.
    It should be noted that the \f$z\f$ is an energy splitting variable, taken in the lab frame,
    and that quarks and diquarks are treated with the same splitting function.
  */
  /*! 
    \class Gluon_Decayer
    \brief This class is responsible for shifting the partons on their constituent mass-shells and
    for the non-pertubative gluon splittings, before clusters are formed and decayed.
  */

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member variables for Gluon_Decayer
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \var FlavDecayMap Gluon_Decayer::m_options
    \brief A list of flavours accessible in gluon decays, their relative flavour-dependent weight
    and the range of z-values (energy splitting in decay).
  */
  /*!
    \var double Gluon_Decayer::m_z
    \brief The energy splitting variable of the decay, from Gluon_Decayer::SelectZ(const DecaySpecs*).
  */
  /*!
    \var double Gluon_Decayer::m_phi
    \brief The azimuthal angle in the gluons decay, from Gluon_Decayer::SelectPhi().
  */
  /*!
    \var double Gluon_Decayer::m_mass
    \brief The constituent mass of the selected quark flavour in the individual, filled in by 
    Gluon_Decayer::SelectDecay(const Vec4D & mom).
  */
  /*!
    \var ATOOLS::Flavour Gluon_Decayer::m_flav
    \brief The selected quark flavour in the decay of the gluon, filled in by 
    Gluon_Decayer::SelectDecay(const Vec4D & mom).
  */
  /*!
    \var ATOOLS::Vec4D Gluon_Decayer::m_p1vec
    \brief Four-vector for one of the decay products of the actual gluon decay under consideration.
    Constructed by BuildKinematics(const ATOOLS::Vec4D &).
  */
  /*!
    \var ATOOLS::Vec4D Gluon_Decayer::m_p2vec
    \brief Four-vector for one of the decay products of the actual gluon decay under consideration.
    Constructed by BuildKinematics(const ATOOLS::Vec4D &).
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member methods for Gluon_Decayer
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \fn Gluon_Decayer::Gluon_Decayer()
    \brief The constructor fills the map of available decay options for the gluon. 

    The AHADIC::Constituents are iterated over, each consituent with twice its (constituent) mass
    smaller than the (constituent) gluon mass is put into the m_options.    
  */
  /*!
    \fn Gluon_Decayer::~Gluon_Decayer()
    \brief Empties and deletes the map of flavours and decay characteristics, m_options.
  */
  /*!
    \fn bool Gluon_Decayer::DecayList(Proto_Particle_List * plin=NULL)
    \brief This method steers the full procedure: First, all particles are put on their constituent
    mass shell, by calling Shift(Proto_Particle_List * pl), then all particles are iterated over and 
    eventually decayed, through DecayIt(ATOOLS::Particle *,ATOOLS::Particle *&,ATOOLS::Particle *&). 
    If the decaying particle is a gluon, it is removed from the list and the two offsprings are inserted 
    at the same place.
  */
  /*!
    \fn void Gluon_Decayer::Shift(Proto_Particle_List * pl)
    \brief Iterates over the particle list, extracts the momenta and finds the constituent masses from
    AHADIC::Constituents::Mass(ATOOLS::Flavour) and puts everything in two arrays. These are then passed
    to the AHADIC::Momenta_Stretcher.
  */
  /*!
    \fn bool Gluon_Decayer::DecayIt(Proto_Particle *,Proto_Particle *&,Proto_Particle *&)
    \brief Decays an individual gluon into two quarks and returns "true"

    First, the flavour of the quark pair is chosen by SelectDecay(const ATOOLS::Vec4D &), with
    the gluon momentum as argument. Then the decay kinematics is built, i.e. two four vectors
    m_p1vec and m_p2vec for the two quarks are constructed in BuildKinematics(const ATOOLS::Vec4D &).
    Then the two corresponding new particles are initialised and filled accordingly.

    \todo make the return value more elegant and straighten out the method 
    DecayList(Proto_Particle_List * plin=NULL).
  */
  /*!
    \fn void Gluon_Decayer::SelectDecay(const ATOOLS::Vec4D &)
    \brief Selects the actual decay mode (i.e. the quark flavour) for a gluon with given momentum.

    The decay mode is selected according to the phase space avaiable for it, i.e. by the weight
    of the corresponding z integration, given by ZWeight(const double zmin, const double zmax).
    The boundaries read
    \f[
    z_{\rm min, max} = \frac12\left(1\mp\frac{|\vec p_g|}{E_g}\,
                       \frac{\sqrt{(m_g^2-2m_q^2)^2-4m_q^4}}{m_g^2}\right)\,.
    \f]
    This weight is put into some DecaySpecs as the weight. 
  */
  /*!
    \fn double Gluon_Decayer::ZWeight(const double zmin, const double zmax)
    \brief The relevant part of the integrated gluon splitting function in the boundaries.

    In the interval \f$[z_{\rm min},\,z_{\rm max}]\f$, the integrated splitting function reads
    \f[
    W = (z_{\rm max}-z_{\rm min}) - (z_{\rm max}-z_{\rm min})^2 + 
        \frac23(z_{\rm max}-z_{\rm min})^3\,.
    \f]
    Ar the moment the same splitting function is being used for both gluons splittings into quarks and 
    gluons splitting into diquarks, the difference between the flavours is in the boundaries
    \f[
    z_{\rm min,max} = 
    \frac12\left(1\mp \frac{p_g}{E_g}\frac{\sqrt{(m_g^2-2.m_q^2)^2-4.m_q^4}}{m_g^2}\right)
    \f]
    only, where \f$p_g\f$ and \f$E_g\f$ denote the absolute value of the gluons three-momentum and its energy,
    both in the lab frame, and \f$m_g\f$ and \f$m_q\f$ are the gluon's and quarks mass, respectively.
  */
  /*!
    \fn void Gluon_Decayer::SelectZ(const DecaySpecs*)
    \brief Selects the energy splitting variable in the gluon decay according to the QCD
    splitting function.

    The distribution of the z values follows
    \f[
    P_{g\to q\bar q}(z)\sim z^2+(1-z)^2\,,
    \f]
    where \f$z\f$ is taken in the lab frame.  
  */
  /*!
    \fn void Gluon_Decayer::SelectPhi()
    \brief Selects the azimuthal angle in gluon decay, at the moment according to an isotropic
    distribution. 

    I.e. \f$\phi\f$ is chosen isotropically in \f$[0,2\pi]\f$.
  */
  /*!
    \fn void Gluon_Decayer::BuildKinematics(const ATOOLS::Vec4D &)
    \brief Constructs the decay kinematics for the gluon decay into two quarks, i.e. the vectors
    m_p1vec and m_p2vec.

    Given the gluon momentum as argument, and with help of z that has been given by 
    SelectZ(const DecaySpecs*), the two new energies and their three-momenta read
    \f[
    E_{q,\,\bar q} = \{z,\,(1-z)\} E_g\;,\;\;\;
    p_{q,\,\bar q} = \sqrt{E^2_{q,\,\bar q}-m^2_q}\,.
    \f]
    Consequently, the two opening angles w.r.t. the original gluon momentum are given by
    \f[
    \cos\theta_{q,\,\bar q} = \frac{p_g^2\pm(p_q^2-p_{\bar q}^2)}{2p_gp_{q,\,\bar q}}\,,
    \f]
    where in all cases the \f$p\f$ denote the absolute values of the corresponding three-momentum.
    With \f$\phi\f$ given by SelectPhi() the decay kinematics can now be fully
    reconstructed. Choosing the orthogonal orientations
    \f[
    \vec n_1 = \vec e_z\times \vec n_g
    \;\;\;\mbox{\rm and}\;\;\;
    \vec n_2 = \vec n_g\times\vec n_1
    \f]
    the four momenta read
    \f[
    p^\mu_{q,\,\bar q} = \left(E_{q,\,\bar q},\,
    p_{q,\,\bar q}\left[\cos\theta_{q,\,\bar q}\,\vec n_g
                        \mp\sin\theta_{q,\,\bar q}(\cos\phi\,\vec n_1+
                                                   \sin\phi\,\vec n_2)\right]\right)^T\,. 
    \f]
  */
}
#endif
