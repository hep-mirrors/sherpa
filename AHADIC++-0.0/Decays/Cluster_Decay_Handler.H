#ifndef Cluster_Decay_Handler_H
#define Cluster_Decay_Handler_H

#include "Cluster_Decayer_Base.H"
#include "Cluster_Decay_Analysis.H"
#include "Cluster_Transformer.H"
#include "Blob.H"
#include "Cluster.H"


namespace AHADIC {
  struct cdm {
    enum code {
      QoverM_Isotropic        = 11,
      QoverM_Retain           = 12,
      FourFermions_Isotropic  = 21,
      FourFermions_Retain     = 22
    };
  };

  class Cluster_Decay_Handler {
  private:
    cdm::code                m_cdm;
    Cluster_Decayer_Base   * p_decayer;
    Cluster_Decay_Analysis * p_analysis;
    Cluster_Transformer    * p_transformer;

    Part_List              * p_partlist;
    ATOOLS::Blob           * p_blob;

    bool                                      m_analyse;
    std::map<std::string,ATOOLS::Histogram *> m_histograms;

    bool DecayIt(Cluster *);
    void InitializeAnalysis();
    void AnalyseThis(ATOOLS::Blob *);
  public:
    Cluster_Decay_Handler(Cluster_Transformer *,bool=false);
    ~Cluster_Decay_Handler();
    void DecayClusters(Cluster_List *,ATOOLS::Blob *);
    Part_List    * GetParticles() { return p_partlist; }
  };
  /*!
    \file Cluster_Decay_Handler.H
    \brief Contains the class AHADIC::Cluster_Decay_Handler. 
  */

  /*! 
    \class Cluster_Decay_Handler
    \brief
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member variables for Cluster_Decay_Handler
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \var Cluster_Decay_Handler::
    \brief 
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member methods for Cluster_Decay_Handler
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \fn Cluster_Decay_Handler::Cluster_Decay_Handler()
    \brief Constructor, no task so far.
  */
  /*!
    \fn Cluster_Decay_Handler::~Cluster_Decay_Handler()
    \brief Destructor, no task so far.
  */
  /*!
    \fn void Cluster_Decay_Handler::
    \brief 
  */
}
#endif
