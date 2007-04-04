#ifndef Cluster_Formation_Handler_H
#define Cluster_Formation_Handler_H

#include "Gluon_Decayer.H"
#include "Cluster_Former.H"
#include "Colour_Reconnections.H"
#include "Soft_Cluster_Handler.H"
#include "Cluster.H"
#include "Blob_List.H"
#include "Blob.H"
#include "Histogram.H"
#include "Return_Value.H"
#include <map>


namespace AHADIC {
  class Cluster_Formation_Handler {
  private:
    bool                                m_single_cr,m_double_cr;

    Gluon_Decayer                     * p_gludecayer;
    Cluster_Former                    * p_cformer;
    Colour_Reconnections              * p_recons;
    Soft_Cluster_Handler              * p_softclusters;

    ATOOLS::Blob                      * p_blob;
    Cluster_List                      * p_clulist;
    VectorOfPPLs                        m_partlists;
    std::vector<Cluster_List *>         m_clulists;

    bool                                      m_analyse;
    std::map<std::string,ATOOLS::Histogram *> m_histograms;

    void Reset();
    ATOOLS::Return_Value::code ExtractSinglets();
    ATOOLS::Return_Value::code FormOriginalClusters();
    ATOOLS::Return_Value::code ApplyColourReconnections();
    ATOOLS::Return_Value::code ClustersToHadrons(ATOOLS::Blob_List *);
    ATOOLS::Return_Value::code MergeClusterListsIntoOne();
  public:
    Cluster_Formation_Handler(bool=false);
    ~Cluster_Formation_Handler();

    ATOOLS::Return_Value::code FormClusters(ATOOLS::Blob *,
					    ATOOLS::Blob_List *);

    //! general access routines
    //@{
    Soft_Cluster_Handler * GetSoftClusterHandler() { return p_softclusters; }
    Cluster_List         * GetClusters()           { return p_clulist;  }
    ATOOLS::Blob         * GetFragmentationBlob()  { return p_blob; }
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

    The way the Colour_Reconnections object is used, is steered through m_single_cr (default: true)
    and m_double_cr (default: false).  If the former is true, then colour reconnections inside a colour 
    singlet list are allowed, if the latter is true, also colour reconnections between two colour singlet 
    lists are taken into account.  
  */
  /*!
    \var Soft_Cluster_Handler * Cluster_Formation_Handler::p_softclusters
    \brief
  */
  /*!
    \var ATOOLS::Blob * Cluster_Formation_Handler::p_blob
    \brief
  */
  /*!
    \var Cluster_List * Cluster_Formation_Handler::p_clulist
    \brief
  */
  /*!
    \var VectorOfPPLs Cluster_Formation_Handler::m_partlists
    \brief
  */
  /*!
    \var std::vector<Cluster_List *> Cluster_Formation_Handler::m_clulists
    \brief
  */
  /*!
    \var bool Cluster_Formation_Handler::m_analyse
    \brief
  */
  /*!
    \var std::map<std::string,ATOOLS::Histogram *> Cluster_Formation_Handler::m_histograms
    \brief
  */
  /*!
    \fn void Cluster_Formation_Handler::Reset()
    \brief
  */
  /*!
    \fn ATOOLS::Return_Value::code Cluster_Formation_Handler::ExtractSinglets()
    \brief Reorganises the incoming parton list in lists of colour-singlet configurations.    

    In so doing, each full Sherpa particle is translated into a largely reduced Proto_Particle,
    which contains only flavour and momentum information.  
  */
  /*!
    \fn ATOOLS::Return_Value::code Cluster_Formation_Handler::FormOriginalClusters()
    \brief
  */
  /*!
    \fn ATOOLS::Return_Value::code Cluster_Formation_Handler::ApplyColourReconnections()
    \brief
  */
  /*!
    \fn ATOOLS::Return_Value::code Cluster_Formation_Handler::ClustersToHadrons(ATOOLS::Blob_List *)
    \brief
  */
  /*!
    \fn ATOOLS::Return_Value::code Cluster_Formation_Handler::MergeClusterListsIntoOne()
    \brief
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
    \fn ATOOLS::Return_Value::code Cluster_Formation_Handler::FormClusters(ATOOLS::Blob *,ATOOLS::Blob_List *)
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
    \brief
  */
  /*!
    \fn Cluster_List * Cluster_Formation_Handler::GetClusters()
    \brief
  */
  /*!
    \fn ATOOLS::Blob * Cluster_Formation_Handler::GetFragmentationBlob()
    \brief
  */
}
#endif
