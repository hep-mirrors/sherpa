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
      fragmentation blob.  At this stage, the actual partons are translated into Proto_Particles,
      a class with largely reduced functionality when compared to the proper particles.  Anyway,
      the extraction of singlet is done by calling the method  
      ATOOLS::Return_Value::code Cluster_Formation_Handler::ExtractSinglets(),
      which obviously fills a vector of Proto_Particle lists, Cluster_Formation_Handler::m_partlists. 
    - Then actual protoclusters are formed.  To this end, colour-singlet parton lists are merged, 
      if their invariant mass is less than the sum of the parton constituent masses (keep in mind
      that at this point the gluon is already considered massive).  The gluons are then decayed
      with help of the p_gludecayer.  These two steps are steered by
      ATOOLS::Return_Value::code Cluster_Formation_Handler::FormOriginalClusters().
    - Having formed the clusters, non-perturbative colour-reconnections may be modelled in different ways,
      encoded in the method ATOOLS::Return_Value::code Cluster_Formation_Handler::ApplyColourReconnections().
    - Those clusters that are light enough to be directly transformed or decayed into hadrons are taken 
      care of by the method 
      ATOOLS::Return_Value::code Cluster_Formation_Handler::ClustersToHadrons(ATOOLS::Blob_List *), which adds 
      these primordial hadrons as outgoing particles to the blob.  If the light clusters decay into two
      hadrons, a corresponding blob is added to the blob list.
    - Finally, the remaining clusters are merged into one large list, Cluster_Formation_Handler::p_clulist, 
      through the method ATOOLS::Return_Value::code Cluster_Formation_Handler::MergeClusterListsIntoOne().
  */
  /*!
    \var bool m_single_cr
    \bried
  */
  /*!
    \var bool m_double_cr
    \bried
  */
  /*!
    \var Gluon_Decayer * p_gludecayer
    \bried
  */
  /*!
    \var Cluster_Former * p_cformer
    \bried
  */
  /*!
    \var Colour_Reconnections * p_recons
    \bried
  */
  /*!
    \var Soft_Cluster_Handler * p_softclusters
    \bried
  */
  /*!
    \var ATOOLS::Blob * p_blob
    \bried
  */
  /*!
    \var Cluster_List * p_clulist
    \bried
  */
  /*!
    \var VectorOfPPLs m_partlists
    \bried
  */
  /*!
    \var std::vector<Cluster_List *> m_clulists
    \bried
  */
  /*!
    \var bool m_analyse
    \brief
  */
  /*!
    \var std::map<std::string,ATOOLS::Histogram *> m_histograms
    \brief
  */
  /*!
    \fn void Cluster_Formation_Handler::Reset()
    \brief
  */
  /*!
    \fn ATOOLS::Return_Value::code Cluster_Formation_Handler::ExtractSinglets()
    \brief
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
    \brief
  */
  /*!
    \fn Cluster_Formation_Handler::~Cluster_Formation_Handler()
    \brief
  */
  /*!
    \fn ATOOLS::Return_Value::code Cluster_Formation_Handler::FormClusters(ATOOLS::Blob *,ATOOLS::Blob_List *)

    \brief
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
  };
}
#endif
