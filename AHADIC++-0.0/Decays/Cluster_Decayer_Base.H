#ifndef Cluster_Decayer_Base_H
#define Cluster_Decayer_Base_H

#include "Cluster_Part.H"
#include "Hadron_Part.H"
#include "Transitions.H"
#include "Cluster.H"


namespace AHADIC {
  class Cluster_Decayer_Base {
  private:
    All_Single_Transitions * p_stransitions;
    Cluster_Part           * p_cdecs;
    Hadron_Part            * p_chads;
 
    int             m_test;
    double          m_offset;
    ATOOLS::Flavour m_had1, m_had2;
    
    void TestOffSprings(Cluster *);
    void TreatHadDecay(Cluster *,Part_List *);
  public:
    Cluster_Decayer_Base(Cluster_Part *,Hadron_Part *);
    ~Cluster_Decayer_Base();
    bool Treat(Cluster *,Part_List *);
  };
  /*!
    \file Cluster_Decayer_Base.H
    \brief Contains the class AHADIC::Cluster_Decayer_Base. 
  */

  /*! 
    \class Cluster_Decayer_Base
    \brief This class is responsible for the organisation of individual cluster decays, for which it owns pointers
    to a AHADIC::Cluster_Part and a AHADIC::Hadron_Part object.

    The cluster decays are organised through a single method, namely Treat(Cluster *,Part_List *).
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member variables for Cluster_Decayer_Base
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \var Cluster_Decayer_Base::
    \brief
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member methods for Cluster_Decayer_Base
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \fn Cluster_Decayer_Base::Cluster_Decayer_Base()
    \brief Constructor, no task so far.
  */
  /*!
    \fn Cluster_Decayer_Base::~Cluster_Decayer_Base()
    \brief Destructor, no task so far.
  */
  /*!
    \fn void Cluster_Decayer_Base::
    \brief 
  */
}
#endif
