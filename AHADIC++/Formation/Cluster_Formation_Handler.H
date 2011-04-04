#ifndef AHADIC_Formation_Cluster_Formation_Handler_H
#define AHADIC_Formation_Cluster_Formation_Handler_H

#include "AHADIC++/Formation/Gluon_Decayer.H"
#include "AHADIC++/Formation/Cluster_Former.H"
#include "AHADIC++/Formation/Colour_Reconnections.H"
#include "AHADIC++/Tools/Soft_Cluster_Handler.H"
#include "AHADIC++/Tools/Cluster.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Math/Histogram.H"
#include <map>


namespace AHADIC {
  class Cluster_Formation_Handler {
  private:
    bool                                m_single_cr,m_double_cr;

    Gluon_Decayer                     * p_gludecayer;
    Cluster_Former                    * p_cformer;
    Colour_Reconnections              * p_recons;
    Soft_Cluster_Handler              * p_softclusters;

    Cluster_List                      * p_clulist;
    ListOfPPLs                          m_partlists;
    std::vector<Cluster_List*>          m_clulists;

    bool                                      m_analyse;
    std::map<std::string,ATOOLS::Histogram *> m_histograms;

    bool ShiftOnMassShells();
    bool ShiftList(Proto_Particle_List *);
    Proto_Particle_List * SelectFromList(ListOfPPLs * lppl,Proto_Particle_List * ppl=NULL);
    bool ExtractSinglets(ATOOLS::Blob *);
    void SetInfoTagForPrimaryParticle(Proto_Particle * proto,leading::code lead) const;
    bool FormOriginalClusters();
    bool ApplyColourReconnections();
    bool ClustersToHadrons(ATOOLS::Blob *);
    bool MergeClusterListsIntoOne();
  public:
    Cluster_Formation_Handler(Cluster_List*, bool=false);
    ~Cluster_Formation_Handler();

    int FormClusters(ATOOLS::Blob *);
    void Reset();
  };
  /*!
    \file Cluster_Formation_Handler.H
    \brief Contains the class AHADIC::Cluster_Formeration_Handler.
  */
  /*!
    \class Cluster_Formation_Handler
    \brief The Cluster_Formation_Handler organises the formation of clusters at the end of
    the parton shower phase.  

    The cluster formation employs the following steps:
    - First, singlets are extracted from the partons, which are the incoming particles in the
      fragmentation blob.  At this stage, each actual parton is translated into a Proto_Particle,
      a class with largely reduced functionality when compared to the proper particles.  Anyway,
      the extraction of singlet is done by calling the method  
      Cluster_Formation_Handler::ExtractSinglets(),
      which obviously fills a vector of Proto_Particle lists, Cluster_Formation_Handler::m_partlists. 
    - Then actual protoclusters are formed.  To this end, colour-singlet parton lists are merged, 
      if their invariant mass is less than the sum of the parton constituent masses (keep in mind
      that at this point the gluon is already considered massive).  The gluons are then decayed
      with help of the p_gludecayer.  These two steps are steered by
      Cluster_Formation_Handler::FormOriginalClusters().
    - Having formed the clusters, non-perturbative colour-reconnections may be modelled in different ways,
      encoded in the method Cluster_Formation_Handler::ApplyColourReconnections().
    - Those clusters that are light enough to be directly transformed or decayed into hadrons are taken 
      care of by the method Cluster_Formation_Handler::ClustersToHadrons(ATOOLS::Blob_List *), which adds 
      these primordial hadrons as outgoing particles to the blob.  If the light clusters decay into two
      hadrons, a corresponding blob is added to the blob list.
    - Finally, the remaining clusters are merged into one large list, Cluster_Formation_Handler::p_clulist, 
      through the method Cluster_Formation_Handler::MergeClusterListsIntoOne().

    Their sequence is organised in Cluster_Formation_Handler::FormClusters(ATOOLS::Blob *,ATOOLS::Blob_List *).
  */
  /*!
    \var bool Cluster_Formation_Handler::m_single_cr
    \brief Steering the way colour reconnections are modelled, colour reconnections inside a colour
    singlet list.

    If true (the default choice) colour reconnections inside an original colour singlet list are allowed, 
    i.e. quarks inside such a list (after the gluon decay) can be reordered, such that they end in different 
    clusters than the original colour assignment would suggest.
  */
  /*!
    \var bool Cluster_Formation_Handler::m_double_cr
    \brief Steering the way colour reconnections are modelled, colour reconnections between two colour
    singlet lists.

    The default value is false.
  */
  /*!
    \var Gluon_Decayer * Cluster_Formation_Handler::p_gludecayer
    \brief A pointer to the Gluon_Decayer, employed in the formation of clusters.
  */
  /*!
    \var Cluster_Former * Cluster_Formation_Handler::p_cformer
    \brief A pointer to the Cluster_Former, employed in the formation of clusters.
  */
  /*!
    \var Colour_Reconnections * Cluster_Formation_Handler::p_recons
    \brief A pointer to the Clour_Reconnections, employed in the formation of clusters.

    The way the Colour_Reconnections object is used, is steered through 
    Cluster_Formation_Handler::m_single_cr (default: true) and 
    Cluster_Formation_Handler::m_double_cr (default: false).  If the former is true, then colour 
    reconnections inside a colour singlet list are allowed, if the latter is true, also colour 
    reconnections between two colour singlet lists are taken into account.  
  */
  /*!
    \var Soft_Cluster_Handler * Cluster_Formation_Handler::p_softclusters
    \brief A pointer to the Soft_Cluster_Handler.  It is employed to check if transitions of the 
    primoridal clusters to hadrons are neccessary.
  */
  /*!
    \var ATOOLS::Blob * Cluster_Formation_Handler::p_blob
    \brief The fragmentation blob, containing all final state partons as incoming particles.  
    While AHADIC progresses, clusters and/or hadrons are added as outgoing particles.
  */
  /*!
    \var Cluster_List * Cluster_Formation_Handler::p_clulist
    \brief The final (merged) list of all clusters that need further decays.  

    This list emerges from the colour-sinlget cluster lists, Cluster_Formation_Handler::m_clulists,
    after the method Cluster_Formation_Handler::MergeClusterListsIntoOne() has been called.
  */
  /*!
    \var ListOfPPLs Cluster_Formation_Handler::m_partlists
    \brief This vector contains all colour-singlet lists of proto-particles.
  */
  /*!
    \var std::vector<Cluster_List *> Cluster_Formation_Handler::m_clulists
    \brief The vector of all cluster lists emerging from the colour-singlet lists of proto-particles.

    Ultimately, these lists are merged into one: Cluster_Formation_Handler::p_clulist.  In the meantime,
    all operations of cluster formation are applied on these individual lists: As far as possible, momenta 
    are shuffled inside these lists, if this becomes neccessary in forced transitions from light clusters 
    to hadrons, colur reconnections are by default applied inisde these lists, etc..
  */
  /*!
    \var bool Cluster_Formation_Handler::m_analyse
    \brief A switch to switch the monitoring analysis on.
  */
  /*!
    \var std::map<std::string,ATOOLS::Histogram *> Cluster_Formation_Handler::m_histograms
    \brief A map of all histograms, eventually being filled when cluster formation is monitored.
  */
  /*!
    \fn void Cluster_Formation_Handler::Reset()
    \brief Resets all lists, in particluar: Cluster_Formation_Handler::m_clulists, 
    Cluster_Formation_Handler::p_clulist, Cluster_Formation_Handler::m_partlists before a new event
    is hadronised.
  */
  /*!
    \fn bool Cluster_Formation_Handler::ExtractSinglets()
    \brief Reorganises the incoming parton list in lists of colour-singlet configurations.    

    In so doing, each full Sherpa particle is translated into a largely reduced Proto_Particle,
    which contains only flavour and momentum information.  
  */
  /*!
    \fn bool Cluster_Formation_Handler::FormOriginalClusters()
    \brief Forms the original clusters.

    This proceeds in the following steps:
    - First of all, all singlet lists are checked, whether their actual mass is larger than the sum
      of their respective constituent masses.  If this is not the case, the respective Proto_Particle 
      list is merged with another of list.
    - Then the gluons in each remaining Proto_Particle list are decayed through the method
      Gluon_Decayer::DecayList(Proto_Particle_List * plin=NULL).  In case, problems are encountered
      in the gluon decays, the method returns "false".
    - If this was successful, the only method of the Cluster_Former,
      Cluster_Former::ConstructClusters(Proto_Particle_List *,Cluster_List *) is invoked to build the 
      primoridal clusters by just adding bit-by-bit two colour-adjacent (di-)quarks into a cluster. 
    - Eventually, the mass and number histograms "Cluster_Mass_Formation" and "Cluster_Number_Formation"
      are now being filled.
  */
  /*!
    \fn bool Cluster_Formation_Handler::ApplyColourReconnections()
    \brief Here, the colourb reconnections are steered.

    In Ahadic, non-perturbative colour reconnections are modelled by reshuffling the constituents
    of clusters, i.e. by exchanging two (anti-)(di-)quarks of two different clusters.  By defautl, i.e.
    if Cluster_Formation_Handler::m_single_cr is switched on, this happens within clusrte lists emerging
    from parton colour singlet configurations only; it is realised by invoking the method
    Colour_Reconnections::Singlet_CR(Cluster_List * clin) for each of the clour singlet lists.
    In addition, such exchanges can also happen between colour singlets (this case has to be
    switched on explicitly), then for each possible pair of singlet lists, the method 
    Colour_Reconnections::Two_Singlet_CR(Cluster_List * cl1,Cluster_List * cl2) is used.

    Eventually, the mass and number histograms "Cluster_Mass_Reconnections" are now being filled.
  */
  /*!
    \fn bool Cluster_Formation_Handler::ClustersToHadrons(ATOOLS::Blob_List *)
    \brief Manages the forced transition of light clusters into hadrons.

    The procedure here is as follows:
    - First, all individual colour singlet cluster lists are iterated over.  For each list, forced 
      cluster-to-hadron transitions are taken care off by 
      Soft_Cluster_Handler::TreatClusterList(Cluster_List *,ATOOLS::Blob *).
    - Then all clusters of the list are iterated over.  If no hadron decay was induced for the respective
      cluster, it is added as outgoing particle to the fragmentation blob.  It should be understood that 
      this includes the option that a cluster-to-hadron transition was enforced.  If a mixed decay was 
      enforced, the so far void pointers to its offsprings are added to the cluster list, otherwise,
      the cluster is added as outgoing particle to the fragmentation blob and a cluster-to-hadron
      decay blob is added to the events Blob_List.
    - Eventually, the mass and number histograms "Cluster_Mass_Transformed" and "Cluster_Number_Transformed"
      are now being filled.

    \todo: Understand the first step better.  
  */
  /*!
    \fn bool Cluster_Formation_Handler::MergeClusterListsIntoOne()
    \brief Here, all cluster list are merged into the global list.
  */
  /*!
    \fn Cluster_Formation_Handler::Cluster_Formation_Handler(bool=false)
    \brief Initialises the following objects: Gluon_Decayer, Cluster_Former, Colour_Reconnections, 
    Soft_Cluster_Handler, and cluster list.  If, in addition, the argument is true, then the 
    constructor initialises a number of histograms monitoring the formation phase of the code.
  */
  /*!
    \fn Cluster_Formation_Handler::~Cluster_Formation_Handler()
    \brief Writes out the monitoring histograms, if existent, and deletes the Gluon_Decayer, the 
    Cluster_Former, the Colour_Reconnections, the Soft_Cluster_Handler, and the cluster list.
  */
  /*!
    \fn int Cluster_Formation_Handler::FormClusters(ATOOLS::Blob *,ATOOLS::Blob_List *)
    \brief The central method of the class, steers the formation of clusters through the
    steps outlined above (and, in short, below).  

    In short, the steps taken during cluster formation are:
    - the construction of colour singlet configurations from the parton list, realised through
      Cluster_Formation_Handler::ExtractSinglets();
    - the shift of the gluons on a new constituent mass shell and their subsequent decay, followed
      by the formaton of proto-clusters based on the singlets alone, realised by
      Cluster_Formation_Handler::FormOriginalClusters();
    - eventual colur-reconnections, resulting in a re-shuffling of the cluster constituents and
      momenta, done in Cluster_Formation_Handler::ApplyColourReconnections();
    - the formation of first hadrons from light clusters, and corresponding shifts of the
      four-momenta of clusters and hadrons; this is done in  
      Cluster_Formation_Handler::ClustersToHadrons(ATOOLS::Blob_List *);
    - finally the merging of the original now hadronised colour-singlet lists into one global
      list of clusters in Cluster_Formation_Handler::MergeClusterListsIntoOne().
  */
  /*!
    \fn Soft_Cluster_Handler * Cluster_Formation_Handler::GetSoftClusterHandler()
    \brief Returns the pointer to the Soft_Cluster_Handler, used in the decay phase.
  */
  /*!
    \fn Cluster_List * Cluster_Formation_Handler::GetClusters()
    \brief Returns the pointer to the Cluster_List, used in the decay phase.
  */
  /*!
    \fn ATOOLS::Blob * Cluster_Formation_Handler::GetFragmentationBlob()
    \brief Not used.
  */
}
#endif
