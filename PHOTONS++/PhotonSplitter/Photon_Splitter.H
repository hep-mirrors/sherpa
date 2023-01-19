#ifndef PHOTONS_PhotonSplitter_Photon_Splitter_H
#define PHOTONS_PhotonSplitter_Photon_Splitter_H

#include "PHOTONS++/PhotonSplitter/Sudakov.H"
#include <cstddef>

namespace ATOOLS {
  class Blob;
}

namespace PHOTONS {

  class Photon_Splitter {

  private:
    int m_on;
    Sudakov m_sudakov;

  public:
    Photon_Splitter();
    Photon_Splitter(int mode);
    ~Photon_Splitter();
    bool SplitPhotons(ATOOLS::Blob *);
  };

  /*!
    \file Photon_Splitter.H
    \brief Contains the class PHOTONS::Photon_Splitter
  */

  /*! 
    \class Photon_Splitter
    \brief The steering class for the soft-photon splitting algorithm in Sherpa
  */
  /////////////////////////////////////////////////////////////////////////////////////////
  // Description of member variables for Photon_Splitter
  /////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \var static int m_on
    \brief whether the photon splitting algorithm is active 
  */
  /*!
    \var Sudakov m_sudakov
    \brief the Sudakov object which implements the veto algorithm 
  */

  /////////////////////////////////////////////////////////////////////////////////////////
  // Description of member methods for Photon_Splitter
  /////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \fn bool PhotonSplitter::SplitPhotons(Blob *)
    \brief processes blob and hands the relevant information to the Sudakov to be treated
  */
}
#endif
