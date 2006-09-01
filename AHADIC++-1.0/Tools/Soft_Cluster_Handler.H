#ifndef Soft_Cluster_Handler_H
#define Soft_Cluster_Handler_H

#include "Transitions.H"
#include "Hadronisation_Parameters.H"
#include "Cluster.H"
#include "Blob.H"

namespace AHADIC {
  struct stm {
    enum code {
      simplemass      = 1,
      massXwaves      = 2,
      masswidth       = 11,
      masswidthXwaves = 12
    };
  };

  struct dtm {
    enum code {
      PS_only      = 0,
      waves_PS     = 1
    };
  };
  
  

  class Soft_Cluster_Handler {
  private:
    stm::code            m_stmode;
    dtm::code            m_dtmode;
    Single_Transitions * p_singletransitions;
    Double_Transitions * p_doubletransitions;
    
    Pair_Popper        * p_popper;
    double               m_alpha,m_fraction;
    double               m_offset1,m_offset2,m_photonenergy;
    
    std::map<Cluster *,std::vector<ATOOLS::Flavour> >           m_ctrans;
    std::map<Cluster *,std::vector<ATOOLS::Flavour> >::iterator m_ctit;

    void   CheckCluster(Cluster *,bool);
    double DecayWeight(Cluster *,ATOOLS::Flavour &,ATOOLS::Flavour &);
    double TransformWeight(Cluster *,ATOOLS::Flavour &,bool=false);
    void   FixHHDecay(Cluster *,ATOOLS::Flavour &,ATOOLS::Flavour &);


    void   ShiftMomenta(const int);
    double SimpleMassCriterion(Single_Transition_List *,ATOOLS::Flavour &,
			       const double,bool=false);
    double MWCriterion(Single_Transition_List *,ATOOLS::Flavour &,
		       const double,bool=false);
  public:
    Soft_Cluster_Handler();
    ~Soft_Cluster_Handler();
    bool           TreatClusterList(Cluster_List *,ATOOLS::Blob *);

  };

  /*!
    \file Soft_Cluster_Handler.H
    \brief Contains the class AHADIC::Soft_Cluster_Handler. 
  */

  /*! 
    \class Soft_Cluster_Handler
    \brief Constructs a cluster list from a particle list, through the (so far) only relevant method
    void ConstructClusters(Part_List * plin,Cluster_List * clout).
  */

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member variables for Soft_Cluster_Handler
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member methods for Soft_Cluster_Handler
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \fn Soft_Cluster_Handler::Soft_Cluster_Handler()
    \brief Constructor, no task so far.
  */
  /*!
    \fn Soft_Cluster_Handler::~Soft_Cluster_Handler()
    \brief Destructor, no task so far.
  */
  /*!
    \fn void Soft_Cluster_Handler::TransformClusters(Part_List * plin,Cluster_List * clout)
    \brief 
  */

}
#endif
