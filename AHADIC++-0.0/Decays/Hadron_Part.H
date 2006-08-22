#ifndef Hadron_Part_H
#define Hadron_Part_H

#include "Transitions.H"
#include "Cluster.H"


namespace AHADIC {
  struct chtrans {
    enum code {
      HH_only      = 0,
      CH_incl      = 1,
      CH_plus_Anni = 2
    };
  };

  struct hadsel {
    enum code {
      keep    = 0,
      newpair = 1
    };
  };

  class Hadron_Part {
  protected:
    chtrans::code m_cht;
    hadsel::code  m_hadsel;
    
    All_Single_Transitions * p_stransitions;
    All_Double_Transitions * p_dtransitions;
    virtual ATOOLS::Return_Value::code HHDecay(Cluster *,Part_List *,
					       ATOOLS::Flavour &,ATOOLS::Flavour &);
    virtual ATOOLS::Return_Value::code CHDecayOne(Cluster *,Part_List *,
						  Cluster *,ATOOLS::Flavour &);
    virtual ATOOLS::Return_Value::code CHDecayBoth(Cluster *,Part_List *,
						   ATOOLS::Flavour &,ATOOLS::Flavour &);
    virtual ATOOLS::Return_Value::code CheckDecayKinematics(Cluster *,
							    ATOOLS::Flavour &,ATOOLS::Flavour &);
    virtual ATOOLS::Return_Value::code CheckDecayKinematics(Cluster *,Cluster *,
							    ATOOLS::Flavour &,ATOOLS::Flavour &);
    virtual void TwoHadronDecay(Cluster *,Part_List *,
				ATOOLS::Flavour &,ATOOLS::Flavour &)=0;
  public:
    Hadron_Part();
    virtual ~Hadron_Part() {}
    virtual ATOOLS::Return_Value::code RedoDecay(Cluster *,Part_List *,int,
						 ATOOLS::Flavour &,ATOOLS::Flavour &);
    virtual ATOOLS::Return_Value::code ForcedDecay(Cluster *,Part_List *);
  };


  class Isotropic : public Hadron_Part {
  private:
    double m_smearparameter;
    void TwoHadronDecay(Cluster *,Part_List *,
			ATOOLS::Flavour &,ATOOLS::Flavour &);
  public:
    Isotropic();
    ~Isotropic() {}
  };

  class Retain : public Hadron_Part {
    double m_smearparameter;
    void TwoHadronDecay(Cluster *,Part_List *,
			ATOOLS::Flavour &,ATOOLS::Flavour &);
  public:
    Retain();
    ~Retain() {}
  };
  /*!
    \file Hadron_Part.H
    \brief Contains the virtual class AHADIC::Hadron_Part, its specific implementeations AHADIC::Isotropic and
    AHADIC::Retain and the structs chtrans and hadsel steering their mode of operation. 
  */

  /*! 
    \class Hadron_Part
    \brief The virtual mother class for binary cluster transitions involving hadrons. 

    It becomes effective, when in the cluster decays described by AHADIC::Cluster_Part one of the emerging
    objects is light. This is usually the case, when the emerging object is lighter than the heaviest hadron
    with a flavour wave function component equal to the flavour composition of the object in question, i.e.
    if one of the two emerging objects satisfies the 
  */
  /*! 
    \class Isotropic
    \brief A specific realisation of the virtual class AHADIC::Hadron_Part, in which the two emerging objects
    are distributed isotropically in space, independently of whether they are clusters or hadrons.
  */
  /*! 
    \class Retain
    \brief A specific realisation of the virtual class AHADIC::Hadron_Part, in which the two emerging objects
    are (roughly) distributed in space according to their distribution given by the AHADIC::Cluster_Part.
 
    \todo This class is not yet implemented.
  */
  /*! 
    \struct chtrans
    \brief Possible cluster transition modes.

    Cases:
    - 0 : \f$C\to HH'\f$ only;
    - 1 : \f$C\to HH\f$ and \f$C\to C'H\f$ (default choice);
    - 2 : like case 1, plus the option to have a full annihilation of a diquark-diquark cluster into two mesons.
          (not implemented yet).
  */
  /*! 
    \struct hadsel
    \brief Possible hadron selection modes in \f$C\to HH'\f$ transitions.

    Cases:
    - 0 : Keep flavour content from pair-popping;
    - 1 : Select a new hadron-pair to account for combined SU(6) weights (default choice).
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member variables for Hadron_Part
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \var chtrans::code Hadron_Part::m_cht
    \brief The mode for the allowed cluster decays involving hadrons.
  */
  /*!
    \var hadsel::code Hadron_Part::m_hadsel
    \brief The mode according to which hadrons are selected.
  */
  /*!
    \var All_Single_Transitions * Hadron_Part::p_stransitions
    \brief This class is used to determine which individual next lightest hadron with a known flavour wave 
    component is taken in case some decay kinematics does not work out. 

    The method for this is All_Single_Transitions::NextLightest(Cluster *, ATOOLS::Flavour &).
  */
  /*!
    \var All_Double_Transitions * Hadron_Part::p_dtransitions
    \brief  This class is used to determine which hadron pair emerges in a corresponding decay of a cluster.

    The method for this is All_Double_Transitions::IsoDecay(Cluster *, ATOOLS::Flavour &, ATOOLS::Flavour &).
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member methods for Hadron_Part
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \fn Hadron_Part::Hadron_Part()
    \brief Constructor, fixes the p_stransitions and p_dtransitions.
  */
  /*!
    \fn Hadron_Part::~Hadron_Part()
    \brief Destructor, no task so far.
  */
  /*!
    \fn virtual void Hadron_Part::RedoDecay(Cluster *,Part_List *,int,ATOOLS::Flavour &,ATOOLS::Flavour &)
    \brief Fixes the decay in case, one of the objects emerging in cluster decays falls into the hadron regime.

    There are different cases to consider:
    - If the cluster decay mode is m_cht=HHonly, i.e. decays into hadrons only, then the decay is performed through 
      sequential calls to All_Double_Transitions::IsoDecay(Cluster *, ATOOLS::Flavour &, ATOOLS::Flavour &), fixing the 
      two hadrons, and TwoHadronDecay(Cluster *,Part_List *,ATOOLS::Flavour &,ATOOLS::Flavour &) actually fixing the 
      kinematics for it and adding the merging hadrons into the particle list for the hadrons that will later on
      be attached to the fragmentation blob as outgoing particles.
    - If the cluster decay mode is not m_cht=HHonly, but the integer argument (mode) equals 3, i.e. both decay products 
      fall into the hadron region, then there are two options:
      -# If the hadron selection mode is m_hadsel=newpair, the sequence of actions basically is identical to
         the case considered above.
      -# If, in contrast, the hadron selection mode is not m_hadsel=newpair then the hadrons onbtained from the
         first trivial range checks are kept and it is only checked, whether the decay is kinematically allowed at
	 all. If not, the hadrons are replaced by lighter ones, see CheckDecayKinematics(Cluster *,Flavour &,Flavour &).
	 Having thus fixed the hadrons, the isotropic decay is performed through the same method as above,
	 TwoHadronDecay(Cluster *,Part_List *,ATOOLS::Flavour &,ATOOLS::Flavour &).
    - If the cluster decay mode is not m_cht=HHonly, but the integer argument (mode) equals 1 or 2, i.e. only one of the 
      decay products falls into the hadron region, then it is checked whether the decay is kinematically allowed,
      this time through CheckDecayKinematics(Cluster *,Cluster *,Flavour &,int). 

    \todo Some more work has to be done for \f$C\to C'H\f$ decays.
  */
  /*!
    \fn virtual bool Hadron_Part::ForcedDecay(Cluster *,Part_List *)
    \brief Forces the decay of a cluster into two hadrons in case, no suitable cluster decay kinematics can be 
    constructed in the Cluster_Part. For more information, cf.
    AHADIC::Hadron_Part::ForcedDecay(Cluster *,Part_List *,int,ATOOLS::Flavour &,ATOOLS::Flavour &)

    This method basically calls the All_Double_Transitions::IsoDecay(Cluster *, ATOOLS::Flavour &, ATOOLS::Flavour &) 
    method, providing two hadrons, and then performs the two-hadron decay with
    TwoHadronDecay(Cluster *,Part_List *,ATOOLS::Flavour &,ATOOLS::Flavour &).
    If no hadron pair can be found due to phase space constraints, the program aborts.

    \todo Here, something nicer should be found, like, e.g., a forced \f$\pi^0\f$ emission.
  */
  /*!
    \fn virtual void Hadron_Part::CheckDecayKinematics(Cluster *,ATOOLS::Flavour &,ATOOLS::Flavour &)
    \brief Checks whether a decay is kinematically allowed, i.e. whether \f$m_C\le m_1+m_2\f$.

    If this is not the case the two hadrons are sequentially replaced with the next-lighter ones with an
    identical relevant flavour wave function component, as given by 
    All_Single_Transitions::NextLightest(Cluster *, ATOOLS::Flavour &).
  */
  /*!
    \fn virtual int Hadron_Part::CheckDecayKinematics(Cluster *,Cluster *,ATOOLS::Flavour &,ATOOLS::Flavour &)
    \brief Here, more work is needed.
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member methods for Isotropic
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \fn Isotropic::Isotropic()
    \brief Constructor, fixes the two modes. At the moment, the default value is m_cht = CH_incl and m_hadsel = newpair.
  */
  /*!
    \fn Isotropic::~Isotropic()
    \brief Destructor, no task so far.
  */
  /*!
    \fn void Isotropic::TwoHadronDecay(Cluster *,Part_List *,ATOOLS::Flavour &,ATOOLS::Flavour &)
    \brief Isotropic decay of a cluster into two hadrons, which are then filled into the particle list to be
    attached to the fragmentaiton blob later on.

    The cluster is boosted into its c.m. system and two appropriate four vectors are constructed:
    - If no leading particle is involved, they are distributed isotropically in \f$4\pi\f$. 
    - If one of the two particles in the cluster is a leading particle, and if the smearing parameter is larger 
      than 0, then the decay is not isotropic anymore. Instead a Gaussian smearing of the form
      \f[
      {\cal P}(\cos\theta) \sim 
      \exp\left[-\langle\cos\theta^{(0)}\rangle^2(1-\cos\theta)^2\right]
      \f]
      is performed around the original direction of the two clusters constituent momenta. Here, \f$\langle\cos\theta^{(0)}\rangle\f$ 
      is the smearing parameter.

    The cluster and the two momenta are boosted back into the lab frame and then particles
    are initialised and filled accordingly (info flag "P"=primary hadron), which are then added to the Part_List.
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member methods for Retain
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \fn Retain::Retain()
    \brief Constructor, no task so far.
  */
  /*!
    \fn Retain::~Retain()
    \brief Destructor, no task so far.
  */
  /*!
    \fn void Retain::TwoHadronDecay(Cluster *,Part_List *,ATOOLS::Flavour &,ATOOLS::Flavour &)
    \brief Decay of a cluster into two hadrons, respecting the original directions in the cluster
    decay as given by the Cluster_Part. 

    - the cluster is boosted into its c.m. system and two appropriate four vectors are constructed,
      which are oriented along the three-momenta of the two offsprings obtained in the test decay in Clusterr_Part;
    - the cluster and the two momenta are boosted back into the lab frame and then particles
      are initialised and filled accordingly (info flag "P"=primary hadron), which are then added to the Part_List;
    - the leading particle treatment is identical to the one in Isotropic. A Gaussian smearing of the form
      \f[
      {\cal P}(\cos\theta) \sim 
      \exp\left[-\langle\cos\theta^{(0)}\rangle^2(1-\cos\theta)^2\right]
      \f]
      is performed around the original direction of the two clusters constituent momenta. Here, \f$\langle\cos\theta^{(0)}\rangle\f$ 
      is the smearing parameter.
  */
}
#endif
