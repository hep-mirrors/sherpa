#ifndef Cluster_Decay_Handler_H
#define Cluster_Decay_Handler_H

#include "Cluster_Part.H"
#include "Soft_Cluster_Handler.H"
#include "Blob_List.H"
#include "Cluster.H"
#include "Cluster_Decay_Analysis.H"
#include "Return_Value.H"


namespace AHADIC {
  class Cluster_Decay_Handler {
  private:
    bool                     m_cib;

    Soft_Cluster_Handler   * p_softclusters;
    Cluster_Part           * p_clus;

    Cluster_Decay_Analysis * p_analysis;
 
    ATOOLS::Blob *  DecayIt(Cluster *);
    void            AnalyseThis(ATOOLS::Blob *);
  public:
    Cluster_Decay_Handler(Soft_Cluster_Handler *,bool=false,bool=false);
    ~Cluster_Decay_Handler();

    ATOOLS::Return_Value::code DecayClusters(Cluster_List *,ATOOLS::Blob_List * = NULL);
  };
  /*!
    \file Cluster_Decay_Handler.H
    \brief Contains the class AHADIC::Cluster_Formeration_Handler.
  */
  /*!
    \class Cluster_Decay_Handler
    \brief The Cluster_Decay_Handler organises the decay of clusters at the end of
    the parton shower phase.  
  */
  /*!
    \var bool m_cib
    \brief
  */
  /*!
    \var Soft_Cluster_Handler * p_softclusters
    \brief
  */
  /*!
    \var Cluster_Part * p_clus
    \brief
  */
  /*!
    \var Cluster_Decay_Analysis * p_analysis
    \brief
  */
  /*!
    \fn ATOOLS::Blob * Cluster_Decay_Handler::DecayIt(Cluster *)
    \brief
  */
  /*!
    \fn void Cluster_Decay_Handler::AnalyseThis(ATOOLS::Blob *)
    \brief
  */
  /*!
    \fn Cluster_Decay_Handler::Cluster_Decay_Handler(Soft_Cluster_Handler *,bool=false,bool=false)
    \brief
  */
  /*!
    \fn Cluster_Decay_Handler::~Cluster_Decay_Handler()
    \brief
  */
  /*!
    \fn ATOOLS::Return_Value::code Cluster_Decay_Handler::DecayClusters(Cluster_List *,ATOOLS::Blob_List * = NULL)
    \brief
  */
}
#endif
