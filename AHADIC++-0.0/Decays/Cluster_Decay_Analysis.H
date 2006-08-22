#ifndef Cluster_Decay_Analysis_H
#define Cluster_Decay_Analysis_H

#include "Blob.H"
#include "Cluster.H"
#include "Histogram.H"
#include <map>


namespace AHADIC {
  class Cluster_Decay_Analysis {
  private:
    bool                                      m_analyse;
    std::map<std::string,ATOOLS::Histogram *> m_histograms;
  public:
    Cluster_Decay_Analysis();
    ~Cluster_Decay_Analysis();

    void AnalyseThis(ATOOLS::Blob *);
  };
  /*!
    \file Cluster_Decay_Analysis.H
    \brief Contains the class AHADIC::Cluster_Decay_Analysis. 
  */

  /*! 
    \class Cluster_Decay_Analysis
    \brief
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member variables for Cluster_Decay_Analysis
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \var Cluster_Decay_Analysis::
    \brief 
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member methods for Cluster_Decay_Analysis
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \fn Cluster_Decay_Analysis::Cluster_Decay_Analysis()
    \brief Constructor, no task so far.
  */
  /*!
    \fn Cluster_Decay_Analysis::~Cluster_Decay_Analysis()
    \brief Destructor, no task so far.
  */
  /*!
    \fn void Cluster_Decay_Analysis::
    \brief 
  */
}
#endif
