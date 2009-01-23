#ifndef Soft_Cluster_Handler_H
#define Soft_Cluster_Handler_H

#include "Transitions.H"
#include "Hadronisation_Parameters.H"
#include "Strong_Coupling.H"
#include "Cluster.H"
#include "Histogram.H"
#include <map>

namespace AHADIC {
  struct PTdist {
    enum code {
      sinthet       = 1,
      alphadist     = 2
    };
  };

  class Soft_Cluster_Handler {
  private:
    Strong_Coupling    * p_as;
    Single_Transitions * p_singletransitions;
    Double_Transitions * p_doubletransitions;
    
    PTdist::code         m_ptmode;
    int                  m_isotropic;
    double               m_offset1,m_offset2,m_kappa,m_alpha,m_eta,m_photonenergy;
    double               m_pt2max;
    Cluster_List *       m_clin;

    bool                                      m_ana;
    std::map<std::string,ATOOLS::Histogram *> m_histograms;


    int    CheckCluster(SP(Cluster),bool,bool=false);
    int    EnforceDecay(SP(Cluster)); 
    double DecayWeight(SP(Cluster),ATOOLS::Flavour &,ATOOLS::Flavour &);
    double TransformWeight(SP(Cluster),ATOOLS::Flavour &,bool=false);
    double SelectPT(double PTmax,bool peaked=true);
    void   FixHHDecay(SP(Cluster),ATOOLS::Blob *);


    bool   CheckIfAllowed(Cluster_List * clin,double & E);
    bool   UpdateTransitions(Cluster_List * clin);
    bool   TryToEnforceTransition(Cluster_List * clin);
    bool   ShiftMomenta(Cluster_List *);
  public:
    Soft_Cluster_Handler(Strong_Coupling * as,
			 Single_Transitions * singletransitions,
			 Double_Transitions * doubletransitions,
			 const double offset1,const double offset2,
			 const double kappa,const double alpha,const double eta,
			 const double photonenergy,int isotropic,bool ana=false);
    ~Soft_Cluster_Handler();
    bool   TreatClusterList(Cluster_List *,ATOOLS::Blob *);

  };

  /*!
    \file Soft_Cluster_Handler.H
    \brief Contains the class AHADIC::Soft_Cluster_Handler and the structs
    AHADIC::stm and AHADIC::dtm, which steer the way clusters transform into hadrons. 
  */


  /*! 
    \struct stm 
    \brief To label the way cluster transformation weights are calculated.

    Available options are
    - simplemass      = 1,
    - massXwaves      = 2,
    - masswidth       = 11, and
    - masswidthXwaves = 12.
  */

  /*! 
    \struct dtm 
    \brief  To define the way cluster decay weights are calculated.

    Available options are
    - PS_only      = 0, and
    - waves_PS     = 1
  */
  
  /*! 
    \class Soft_Cluster_Handler
    \brief Transforms or decays clusters into hadrons.

    The basic idea here is the following: for each cluster in a list, a transformation
    probability and a decay probability is calculated, which, obviously, depends on the
    mass of the cluster and its flavor content.  If any of these probabilities yields
    a non-zero result, the cluster is treated accordingly.  If both are non-zero the
    decision is made probabilistically.  

    The respective probabilities are calculated through the methods
    Soft_Cluster_Handler::TransformWeight(Cluster *,ATOOLS::Flavour &,bool=false) and 
    Soft_Cluster_Handler::DecayWeight(Cluster *,ATOOLS::Flavour &,ATOOLS::Flavour &) and the
    corresponding algorithms are documented there.  

    The default values for the calculation of decay/transition probabilities are
    m_stmode=stm::masswidthXwaves and m_dtmode=dtm::waves_PS. 
  */

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member variables for Soft_Cluster_Handler
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \var stm::code Soft_Cluster_Handler::m_stmode
    \brief Defining the cluster to hadron transition probabilities
  */
  /*!
    \var dtm::code Soft_Cluster_Handler::m_dtmode
    \brief Defining the cluster to hadron decay probabilities
  */
  /*!
    \var Single_Transitions * Soft_Cluster_Handler::p_singletransitions
    \brief A pointer to a list of all transitions of clusters with given flavor content
    to lists of hadrons with fitting wave function component and their value.  
  */
  /*!
    \var Double_Transitions * Soft_Cluster_Handler::p_doubletransitions
    \brief A pointer to a list of all transitions of clusters with given flavor content
    to lists of pairs of hadrons with fitting wave function components and their value.  
  */
  /*!
    \var double Soft_Cluster_Handler::m_alpha
    \brief The (tunable) exponent \f$\alpha\f$ in the calculation of the decay weights.

    This parameter is set through 'C->HH_DECAY_EXPONENT'.
  */
  /*!
    \var double Soft_Cluster_Handler::m_offset1
    \brief The (tunable) mass shift \f$\delta_1\f$ in the determination of open transition
    channels.

    This parameter is set through 'TRANSITION_OFFSET'.
  */
  /*!
    \var double Soft_Cluster_Handler::m_offset2
    \brief The (tunable) mass shift \f$\delta_2\f$ in the determination of open decay
    channels.

    This parameter is set through 'DECAY_OFFSET'.
  */
  /*!
    \var double Soft_Cluster_Handler::m_photonenergy
    \brief The minimal photon energy in forced decays of clusters that must transform into 
    hadrons.

    This parameter is set through 'PHOTON_ENERGY'.
  */
  /*!
    \var std::map<Cluster *,std::vector<ATOOLS::Flavour> > Soft_Cluster_Handler::m_ctrans
    \brief A map connecting each cluster with a vector of potential offspring flavors.

    This vector should have length 0, 1, or 2 depending on whether the cluster transform or
    decays into hadrons (1 or 2) or not (0).
  */
  /*!
    \var std::map<Cluster *,std::vector<ATOOLS::Flavour> >::iterator Soft_Cluster_Handler::m_ctit
    \brief An iterator for the map connecting each cluster with a vector of potential 
    offspring flavors.
  */

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member methods for Soft_Cluster_Handler
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \fn Soft_Cluster_Handler::Soft_Cluster_Handler()
    \brief Constructor, no task so far, apart from setting parameter values and filling in
    the members. 
  */
  /*!
    \fn Soft_Cluster_Handler::~Soft_Cluster_Handler()
    \brief Destructor, no task so far.
  */
  /*!
    \fn bool Soft_Cluster_Handler::TreatClusterList(Cluster_List *,ATOOLS::Blob *)
    \brief Checks for clusters in lists, whether they have to transform or decay into
    hadrons.  

    - To start with, each single cluster is checked through the method 
      Soft_Cluster_Handler::CheckCluster(Cluster *,bool), where the bool indicates whether
      this cluster is the only one in the list.  If so, there's no other cluster which
      could, in principle, carry the recoil of the transformation for non-identical masses
      (the regular case) of cluster and emerging hadron.  By the call of CheckCluster,
      the map of emerging flavors related to that cluster is filled.
    - If the list of clusters contains only one cluster, and if there are two flavors attached to
      it, the cluster will decay into these tow particles.  The kinematics and other characteristics
      are filled by Soft_Cluster_Handler::FixHHDecay(Cluster *,ATOOLS::Flavour &,ATOOLS::Flavour &).
      For lists with only one cluster, the method ends here. 
    - If, in contrast, the list of clusters contains more than one cluster, the four momenta may 
      have to be reshuffled in order to accommodate for recoils in eventual cluster to hadron 
      transitions.  This is done by Soft_Cluster_Handler::ShiftMomenta(const int).  Then, again, 
      the clusters are iterated over.  In case there's two flavors attached to a cluster, it will 
      decay into these two particles.  The decay kinematics etc. are then determined through a call
      to Soft_Cluster_Handler::FixHHDecay(Cluster *,ATOOLS::Flavour &,ATOOLS::Flavour &).  If there's
      only one flavor, it must, by construction be a cluster, which is then added as outgoing
      particle to the initial blob.    
  */
  /*!
    \fn void Soft_Cluster_Handler::CheckCluster(Cluster *,bool)
    \brief Construct decay and transformation weights for clusters and corresponding channels. 

    The algorithm here proceeds as follows:
    - First of all, the cluster is checked whether it must disintegrate.  This is the case
      only for relatively light clusters, consisting of two diquarks.  The check is performed 
      and the emerging hadrons are determined through a call to the method
      Single_Transitions::MustDesintegrate(Cluster *,ATOOLS::Flavour &,ATOOLS::Flavour &).
    - Then decay weight and transformation weight plus eventually emerging hadrons are determined
      by calling Soft_Cluster_Handler::DecayWeight(Cluster *,ATOOLS::Flavour &,ATOOLS::Flavour &)
      and Soft_Cluster_Handler::TransformWeight(Cluster *,ATOOLS::Flavour &,bool=false), 
      respectively.  If both are zero, the method terminates here.
    - The respective channel is then chosen probabilistically and the flavor(s) are filled
      into the cluster-flavor map, Soft_Cluster_Handler::m_ctrans.  

    There's a caveat, however.  If a single hadron transition is selected for a cluster that 
    must decay (for instance, if it is the only one in its list), checks are performed to ensure that
    this transition can be performed:
    - The first is to check whether the hadron mass plus the minimal photon energy are smaller 
      than the cluster mass, in order for the cluster to decay into that hadron plus a photon. 
    - If this is not the case, a new transformation is chosen, such that the emerging
      hadron mass is constrained to be smaller than the cluster mass (this is not necessarily
      the case in the first trial).  If the weight for this transformation is again compared
      with the decay weight.  If still the transition wins the race, the cluster is decayed into
      a hadron and a photon.  
  */
  /*!
    \fn double Soft_Cluster_Handler::DecayWeight(Cluster *,ATOOLS::Flavour &,ATOOLS::Flavour &)
    \brief Calculates and returns the summed weight for all decays and selects one.  

    In order to calculate this weight, a cluster must fulfill some requirements:
    - If its mass is smaller than the mass of its two constituents plus twice the lightest
      possible constituent, any open decay channel is taken, irrespective of any other
      constraint.
    - If the combined mass of the heaviest possible hadron pair is smaller than the cluster
      mass plus the decay offset, the decay weight is 0 and no decay channel is selected.
    - If the combined mass of the lightest possible hadron pair is larger than the cluster
      mass, there's no open channel and the weight is 0.

    There's two options to define the total weight for each decay of a cluster with mass \f$M_c\f$
    made of the (di)quarks \f$q\f$ and \f$\bar q\f$, namely
    - m_dtm = PS_only:
      \f[
      wt_{PS} = \frac{\sqrt{[M_c^2-(m_1+m_2)^2][M_c^2-(m_1-m_2)^2]}}{16\pi M_c^3}\cdot
                \left[\frac{4m_1m_2}{M_c^2}\right]^\alpha
      \f]
    - m_dtm = waves_PS:
      \f[
      wt = wt_{PS}\cdot |\Psi_{q\bar q'}(h_1)|^2|\Psi_{q'\bar q}(h_2)|^2
                  \cdot wt_{q'}
      \f]
    Here, \f$m_{1,2}\f$ are the masses of the hadrons \f$h_{1, 2}\f$.  Their wave component
    with fitting flavors is given by the \f$\Psi\f$, and \f$wt_{q'}\f$ is the weight to
    pop the extra flavor \f$q'\f$.  This combined wave function - flavor popping
    probability part of the total weight is encoded in the corresponding entry of
    Double_Transitions.
  */
  /*!
    \fn double Soft_Cluster_Handler::TransformWeight(Cluster *,ATOOLS::Flavour &,bool=false)
    \brief Calculates and returns the weight for the most likely cluster-hadron transition.

    The algorithm here is as follows:
    - Clusters, for which the combined mass of the lightest decay channel is larger than
      their mass are always considered; they will always have a transformation weight
      and thus they will always be transformed.  
    - If the cluster mass is larger than the combined mass of the lightest decay, then
      a check is performed, whether the heaviest hadron in question plus the one-particle
      transition offset is lighter than the cluster mass.  If so, the transformation weight is
      zero.  
    In other words, only those clusters are considered for transitions, which are light enough,
    such that they cannot decay regularly into two hadrons or whose mass is smaller than
    the heaviest hadron (with a flavor component matching their own flavor constellation)
    plus an offset.  

    Then, there are four ways to define the transition weight:
    - m_stm = simplemass and m_stm = massXwaves: 
      These two are realized in 
      Soft_Cluster_Handler::SimpleMassCriterion(Single_Transition_List *,ATOOLS::Flavour &,
                                                const double,bool=false).
    - m_stm = masswidth and m_stm = masswidthXwaves
      These two are realised in 
      Soft_Cluster_Handler::MWCriterion(Single_Transition_List *,ATOOLS::Flavour &,
		                        const double,bool=false)

    In all cases one transition, if needed is selected and its weight returned.  
  */
  /*!
    \fn void Soft_Cluster_Handler::FixHHDecay(Cluster *,ATOOLS::Flavour &,ATOOLS::Flavour &)
    \brief Fixes the kinematics of a cluster decay into two hadrons and saves the
    information in the cluster's offsprings.

    With the \f$\theta\f$ selected according to \f$\sin^\eta\theta\f$ for the decay at hand  - 
    it is taken w.r.t. the original direction of the clusters constituents - the kinematics 
    is fully fixed up to the azimuthal angle which is picked flat at random.  
  */
  /*!
    \fn void Soft_Cluster_Handler::ShiftMomenta(const int)
    \brief This method is used to shift momenta of a given list of clusters in such
    a way that eventually new masses are taken into account.  
  */
  /*!
    \fn double Soft_Cluster_Handler::SimpleMassCriterion(Single_Transition_List *,ATOOLS::Flavour &,
			                                 const double,bool=false)
    \brief This criterion for cluster transitions relies on mass differences and, eventually
    wave function overlaps only.

    The weight for Soft_Cluster_Handler::m_stm=simplemass is given by
    \f[
    wt_M = [(M_c^2-m_h^2)^2+\epsilon]^{-\frac{1}{4}}\,,
    \f]
    where \f$epsilon = 1.e-8\f$ is used to guarantee finite weights.  

    In case the transition mode is Soft_Cluster_Handler::m_stm=massXwaves the weight reads
    \f[
    wt = wt_M \cdot |\Psi_{q\bar q}(h)|^2
    \f]
    for clusters with flavor content \f$q-\bar q\f$ and the hadron wave component \f$\Psi\f$
    corresponding to it.  
  */
  /*!
    \fn double Soft_Cluster_Handler::MWCriterion(Single_Transition_List *,ATOOLS::Flavour &,
		                                 const double,bool=false)
    \brief This criterion for cluster transitions relies on a Breit-Wigner form and, eventually
    wave function overlaps.  

    The weight for Soft_Cluster_Handler::m_stm=masswidth is given by
    \f[
    wt_{MW} = \frac{1}{2M_c}\cdot \frac{(m_h\Gamma_h)^2}{(M_c^2-m_h^2)^2+(m_h\Gamma_h)^2}\,,
    \f]
    where \f$\Gamma_h\f$ is the hadron's width.  If it's zero the smallest width of
    all hadrons in question is taken.

    The weight for Soft_Cluster_Handler::m_stm=masswidthXwaves is then given by
    \f[
    wt = wt_{MW} \cdot |\Psi_{q\bar q}(h)|^2
    \f]
    for clusters with flavor content \f$q-\bar q\f$ and the hadron wave component \f$\Psi\f$
    corresponding to it.
  */
}
#endif
