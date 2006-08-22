#ifndef Cluster_Formation_Handler_H
#define Cluster_Formation_Handler_H

#include "Gluon_Decayer.H"
#include "Cluster_Former.H"
#include "Colour_Reconnections.H"
#include "Cluster_Transformer.H"
#include "Cluster.H"
#include "Blob_List.H"
#include "Blob.H"
#include "Histogram.H"
#include "Return_Value.H"
#include <map>


namespace AHADIC {
  class Cluster_Formation_Handler {
  private:
    bool                         m_single_cr,m_double_cr;

    Gluon_Decayer              * p_gludecayer;
    Cluster_Former             * p_cformer;
    Colour_Reconnections       * p_recons;
    Cluster_Transformer        * p_ctransformer;

    ATOOLS::Blob               * p_blob;
    Cluster_List               * p_clulist;
    std::vector<Part_List *>     m_partlists;
    std::vector<Cluster_List *>  m_clulists;

    bool                                      m_analyse;
    std::map<std::string,ATOOLS::Histogram *> m_histograms;

    void Reset();

    ATOOLS::Return_Value::code ExtractSinglets(ATOOLS::Blob_List * bl);
    ATOOLS::Return_Value::code FormOriginalClusters();
    ATOOLS::Return_Value::code ApplyColourReconnections();
    ATOOLS::Return_Value::code ClustersToHadrons();
    ATOOLS::Return_Value::code MergeClusterListsIntoOne();
  public:
    Cluster_Formation_Handler(bool=false);
    ~Cluster_Formation_Handler();

    ATOOLS::Return_Value::code FormClusters(ATOOLS::Blob *,
					    ATOOLS::Blob_List * bl=NULL);

    //! general access routines
    //@{
    Cluster_Transformer * GetClusterTransformer() { return p_ctransformer; }
    Cluster_List        * GetClusters()           { return p_clulist;  }
    ATOOLS::Blob        * GetFragmentationBlob()  { return p_blob; }
    //@}
  };

  /*!
    \file Cluster_Formation_Handler.H
    \brief Contains the class AHADIC::Cluster_Formation_Handler. 
  */

  /*! 
    \class Cluster_Formation_Handler
    \brief This class steers the formation phase in cluster fragmentation, where the gluons
    are decayed into quark-antiquark pairs and the primordial clusters are formed.

    For these tasks, the Cluster_Formation_Handler uses a number of other objects, namely
    - AHADIC::Gluon_Decayer to move the partons on their non-perturbative (constituent) mass
      shells and to decay the gluons;
    - AHADIC::Cluster_Former to form the first, trial proto-clusters, which are treated by
    - AHADIC::Colour_Reconnections, which models non-perturbative effects on the original partonic
      clolur flow and, hence, on the proto-clusters, which may be revised in this process;
    - AHADIC::Cluster_Transformer, which checks the emerging proto-clusters and eventually
      translates the lightest ones directly into hadrons.

    The sequence of these actions is visible to the outside only through the public method
    ATOOLS::Blob * FormClusters(ATOOLS::Blob_List * bl=NULL), which eventually returns 
    a newly initialised and filled ATOOLS::Blob of the type ATOOLS::btp::Cluster_Formation.  
  */

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member variables for Cluster_Formation_Handler
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \var bool Cluster_Formation_Handler::m_single_cr
    \brief Flag which enables (true) or disables (false) the non-perturbative colour reconnections
    on the level of one list of partons that have been in a colour singlet before the gluons have
    been decayed, i.e. of the type \f$qggggg\bar q\f$. Default value is true.
  */
  /*!
    \var bool Cluster_Formation_Handler::m_double_cr
    \brief  Flag which enables (true) or disables (false) the non-perturbative colour reconnections
    on the level of two list of partons that have been in two colour singlets before the gluons have
    been decayed. Default value is false.
  */
  /*!
    \var Gluon_Decayer * Cluster_Formation_Handler::p_gludecayer
    \brief The pointer to the Gluon_Decayer object, which is responsible for the forced decay of
    the gluons into quark-antiquark pairs at the end of the parton shower. It is initialised 
    (and deleted) by the Cluster_Formation_Handler.
  */
  /*!
    \var Cluster_Former * Cluster_Formation_Handler::p_cformer
    \brief The pointer to the Cluster_Former object, which is responsible for the formation of 
    the protoclusters from the quarks and antiquarks emerging from the parton shower or the forced
    gluon decays. It is initialised (and deleted) by the Cluster_Formation_Handler.
  */
  /*!
    \var Colour_Reconnections * Cluster_Formation_Handler::p_recons
    \brief The pointer to the Colour_Reconnections object, which is responsible for the modelling
    of non-perturbative shuffles in the colour configuration of the partons, leading eventually to
    a modification of the protoclusters' quark content. Its actions are steered by the two flags
    m_single_cr and m_double_cr. This object is initialised (and deleted) by the 
    Cluster_Formation_Handler.
  */
  /*!
    \var Cluster_Transformer * Cluster_Formation_Handler::p_ctransformer
    \brief The pointer to the Cluster_Transformer object, which is responsible for the eventual
    transformation of light primordial clusters emerging in the formation process into hadrons.
    These hadrons are then added to the "Fragmentation" ATOOLS::Blob, which is being filled by
    the full cluster fragmentation program. Anyways, this object is initialised (and deleted) by 
    the Cluster_Formation_Handler.
  */
  /*!
    \var ATOOLS::Blob * Cluster_Formation_Handler::p_blob
    \brief The pointer to the ATOOLS::Blob, which carries the incoming on-shell (i.e. massless) 
    partons and the outgoing primordial hadrons. It is initialised on an event-by-event base in the
    Cluster_Formation_Handler and being filled by various methods inside this class or others inside 
    the AHADIC namespace. 
  */
  /*!
    \var Cluster_List * Cluster_Formation_Handler::p_clulist
    \brief The list of clusters that are formed by the actions of the Cluster_Formation_Handler.
    This list is handed over at the end of the formation stage to the decay classes which in the end
    decay clusters into hadrons and thus empty the list of "active" clusters.  
  */
  /*!
    \var std::vector<Part_List *> Cluster_Formation_Handler::m_partlists
    \brief A list of particle lists, each of which corresponds to a colour singlet. 

    These lists are being constructed on an event-by-event basis by the method
    Cluster_Formation_Handler::ExtractSinglets(ATOOLS::Blob_List * bl).
  */
  /*!
    \var std::vector<Cluster_List *> Cluster_Formation_Handler::m_clulists
    \brief A list of cluster lists, each of which emerges from one of the colour singlet particle lists
    in Cluster_Formation_Handler::m_partlists.
  */
  /*!
    \var bool Cluster_Formation_Handler::m_analyse
    \brief A flag that steers whether (true) or not (false) an analyses of the primordial
    clusters in being performed. It is set from the outside, through the constructor.
  */
  /*!
    \var std::map<std::string,ATOOLS::Histogram *> Cluster_Formation_Handler::m_histograms
    \brief If an analysis is being performed, it will act on a list of histograms, available through
    this member. 

    Right now, the following histograms are being filled:
    - masses of protoclusters, immediately after formation ("Cluster_Mass_Formation")
    - masses of protoclusters, after colour reconnections have taken place 
      ("Cluster_Mass_Reconnections")
    - masses of clusters after the transition single light clusters -> single hadrons have taken
      place ("Cluster_Mass_Transformed")
    - the number of protoclusters ("Cluster_Number_Formation"), and
    - the number of clusters actual entering the cluster decay phase ("Cluster_Number_Transformed").
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member methods for Cluster_Formation_Handler
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \fn Cluster_Formation_Handler::Cluster_Formation_Handler(bool=false)
    \brief Constructor, initialises nearly all objects in the Cluster_Formation_Handler, namely
    the Gluon_Decayer, the Cluster_Former, the Colour_Reconnections, the Cluster_Transformer, and
    if the argument is true, also the histograms in Cluster_Formation_Handler::m_histograms.
  */
  /*!
    \fn Cluster_Formation_Handler::~Cluster_Formation_Handler()
    \brief Deletes all objects initialised by the constructor, and if needed, writes out the
    histograms before deleting them.
  */
  /*!
    \fn ATOOLS::Blob * Cluster_Formation_Handler::FormClusters(ATOOLS::Blob_List * bl=NULL)
    \brief The central method that steers the sequence of actions in the cluster formation process.

    After resetting the relevant lists in the formation process and initialising a new Blob to
    attach incoming partons and outgoing primoridal hadrons, it calls, in order,
    - Cluster_Formation_Handler::ExtractSinglets(ATOOLS::Blob_List * bl)
    - Cluster_Formation_Handler::FormOriginalClusters()
    - Cluster_Formation_Handler::ApplyColourReconnections()
    - Cluster_Formation_Handler::ClustersToHadrons()
    - Cluster_Formation_Handler::MergeClusterListsIntoOne()
  */
  /*!
    \fn void Cluster_Formation_Handler::Reset();
    \brief Resets, i.e. empties all particle and clusters lists in 
    Cluster_Formation_Handler::m_partlists and Cluster_Formation_Handler::m_clulists.
  */
  /*!
    \fn void Cluster_Formation_Handler::ExtractSinglets(ATOOLS::Blob_List * bl)
    \brief This method fills the incoming partons into the fragmentation blob and constructs the
    (singlet) particle lists in m_partlists.

    This method has two purposes: 
    - First, it copies all relevant, i.e. coloured, partons in the entire ATOOLS::Blob_List 
      (usually the full event), which do not yet have a "decay"-blob. The original partons are then
      added as incoming ones to the fragmentation blob, and the copies are put into the internal 
      list of particles to work with.
    - Second, these internal particles are sorted into singlets. To do so, all particles are iterated 
      over until a particle with a colour configuration (a,0) is found, i.e. with colour in the
      triplet representation. This marks the begin of a singlet list, the corresponding particle
      is put into the new list and deleted from the original one. The new list is put into the list 
      of singlet lists, Cluster_Formation_Handler::m_partlists. Then the particles in the original
      list are iterated over, until a particle with colour configuration (b,a) is found; this particle
      is put into the new list and deleted from the old one. This is repeated until a particle
      with configuration (0,c) is found, marking the end of the new list. These new lists thus 
      represent linear lists in colour space, clealry not the only possible configuration for colour
      singlets. If particles are left in the original list, after all triplet starts have been 
      dealt with, the starting point for lsits become particles with configurations (a,b) and then 
      the singlets from "rings" in colour space.
  */
  /*!
    \fn void Cluster_Formation_Handler::FormOriginalClusters()
    \brief In this method the Gluon_Decayer is used to first shift all particles on thier constituent
    mass shell, before the gluons are decayed into quark-antiquark pairs. Then the emerging quarks
    are put into colour-neutral clusters.

    This method iterates over all singlet lists in Cluster_Formation_Handler::m_partlists and 
    puts them into the Gluon_Decayer::DecayList(Part_List *) method. There, the partons are put
    on their constituent mass shell and the gluons are decayed into quark-antiquark pairs, which
    in turn carry the corresponding colour indices. Through this decay the particle list is being
    modified: all gluons are erased and replaced by the quarks. Then, the particle list is 
    handed over to Cluster_Former::ConstructClusters(Part_List *,Cluster_List *), where the
    quarks and antiquarks are put into colour-neutral protoclusters. In so doing, new clusters are 
    initialised and put into newly created cluster lists, one such cluster list for each
    particle list. 
  */
  /*!
    \fn void Cluster_Formation_Handler::ApplyColourReconnections()
    \brief This method calls the Colour_Reconnections object to model colour rearrangements
    between two clusters inside one or two singlet lists.

    - If colour reconnections inside one singlet list are enabled, for each cluster list in 
      Cluster_Formation_Handler::m_clulists the method Colour_Reconnections::Singlet_CR(Cluster_List *) 
      is being called , which reshuffles quarks inside the clusters of one cluster list. 
    - If colour reconnections between two singlet lists are enabled, for each pair of cluster lists
      in Cluster_Formation_Handler::m_clulists the (not yet operational) method
      Colour_Reconnections::Two_Singlet_CR(Cluster_List *,Cluster_List *) is being called, 
      which reshuffles quarks in clusters of two cluster lists. 
  */
  /*!
    \fn void Cluster_Formation_Handler::ClustersToHadrons()
    \brief This method calls the Cluster_Transformer object in different ways to transform
    light clusters into primordial hadrons, where neccessary. 

    - First of all, the cluster lists are iterated over and those lists, that consist of one single
      cluster only, are dealt with. Whether this is necessary or not is determined through the method
      Cluster_Transformer::TreatSingleCluster(Cluster_List *,ATOOLS::Blob *), which returns true,
      if the cluster must decay into a hadron. If this is the case, this cluster list is added to the 
      shortest other cluster list; depending on the transformation mode in the Cluster_Transformer,
      the cluster may have been decayed. 
      -# In this case, the offsprings have already been added to the fragmentation blob, and the 
         cluster list is empty, thus not altering the behaviour of the other ones.
      -# Otherwise, the resulting cluster list has more than two entries and is dealt with in the next 
         round.
    - In a second step all non-empty cluster lists are iterated over and are treated with the
      method Cluster_Transformer::TreatClusterList(Cluster_List *,ATOOLS::Blob *). There, each 
      individual cluster is tested whether it must decay into a hadron, which is immediately selected.
      Then momenta are shuffled inside this cluster list to account for potentially altered masses 
      of clusters (and hadrons), without violating four-momentum conservation.
  */
  /*!
    \fn void Cluster_Formation_Handler::MergeClusterListsIntoOne()
    \brief Here, finally, the different cluster lists are merged into a final one, the 
    Cluster_Formation_Handler::p_clulist, which contains all clusters that need to be dealt with
    through appropriate cluster decays. The Cluster_Formation_Handler::m_clulists are thus emptied.
  */
}
#endif
