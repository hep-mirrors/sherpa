#ifndef AHADIC_Decays_Cluster_Decay_Analysis_H
#define AHADIC_Decays_Cluster_Decay_Analysis_H

#include "ATOOLS/Phys/Blob.H"
#include "AHADIC++/Tools/Cluster.H"
#include "ATOOLS/Math/Histogram.H"
#include <map>


namespace AHADIC {
  class Cluster_Decay_Analysis {
  private:
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
    \brief The class Cluster_Decay_Analysis contains and fills histograms concerning the
    composition of the sets of primordial hadrons emerging from fragmentation and their momenta.
  */
  /*!
    \fn Cluster_Decay_Analysis::Cluster_Decay_Analysis()
    \brief Constructor, initialises the histograms.
  */
  /*!
    \fn Cluster_Decay_Analysis::~Cluster_Decay_Analysis()
    \brief Destructor, writes out the histograms.
  */
  /*!
    \fn void Cluster_Decay_Analysis::
    \brief Fills the histograms.
  */
}
#endif
