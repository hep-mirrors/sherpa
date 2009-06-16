#ifndef PHOTONS_MEs_ME_List_H
#define PHOTONS_MEs_ME_List_H

#include "PHOTONS++/Main/Dipole_Type.H"

#include "PHOTONS++/MEs/PHOTONS_ME_Base.H"

#include "PHOTONS++/MEs/Vector_To_Fermion_Fermion.H"
#include "PHOTONS++/MEs/Scalar_To_Fermion_Fermion.H"
#include "PHOTONS++/MEs/Scalar_To_Scalar_Scalar.H"
#include "PHOTONS++/MEs/Z_To_Lepton_Lepton.H"
#include "PHOTONS++/MEs/W_To_Lepton_Neutrino.H"
#include "PHOTONS++/MEs/Tau_To_Lepton_Neutrinos.H"

namespace PHOTONS {
  class ME_List {
    private:
    public:
      enum code {
        // generic types 3-body
        Vector_To_Scalar_Scalar       =  1,
        Vector_To_Fermion_Fermion     =  2,
        Vector_To_Vector_Vector       =  3,

        Scalar_To_Scalar_Scalar       =  11,
        Scalar_To_Fermion_Fermion     =  12,
        Scalar_To_Vector_Vector       =  13,

        Fermion_To_Fermion_Vector     =  21,
        Fermion_To_Fermion_Scalar     =  22,


        // special types
        Z_To_Lepton_Lepton            =  101,
        W_To_Lepton_Neutrino          =  102,
        Tau_To_Lepton_Neutrinos       =  103,

        // unknown
        unknown                       =  0
      };
  };


  class ME_Browser {
    private:
      bool            m_at_all;
      ME_List::code   m_ME;

    public:
      ME_Browser(Particle_Vector_Vector);
      ~ME_Browser() {}

      inline bool GetMEInfo()         { return m_at_all; }
      inline ME_List::code WhichME()  { return m_ME; }
  };


  inline ME_Browser::ME_Browser(Particle_Vector_Vector pvv) {
    // Z_To_Lepton_Lepton
    if ((pvv.size() == 4) &&
        (pvv.at(0).size() == 0) &&
        (pvv.at(1).size() == 1) && (pvv.at(1).at(0)->Flav().Kfcode() == kf_Z) &&
        (pvv.at(2).size() == 2) && (pvv.at(2).at(0)->Flav().IsLepton() == true) &&
        (pvv.at(3).size() == 0))  {
                                    m_at_all = true;
//                                     m_ME = ME_List::Z_To_Lepton_Lepton;
                                    m_ME = ME_List::Vector_To_Fermion_Fermion;
                                  }
    // W_To_Lepton_Neutrino
    else if ((pvv.size() == 4) &&
             (pvv.at(0).size() == 1) && (pvv.at(0).at(0)->Flav().Kfcode() == kf_Wplus) &&
             (pvv.at(1).size() == 0) &&
             (pvv.at(2).size() == 1) && (pvv.at(2).at(0)->Flav().IsLepton() == true) &&
             (pvv.at(3).size() == 1) && (pvv.at(3).at(0)->Flav().IsLepton() == true)) 
                                  {
                                    m_at_all = true;
                                    m_ME = ME_List::W_To_Lepton_Neutrino;
                                  }
    // Tau_To_Lepton_Neutrinos
    else if ((pvv.size() == 4) &&
             (pvv.at(0).size() == 1) && (pvv.at(0).at(0)->Flav().Kfcode() == kf_tau) &&
             (pvv.at(1).size() == 0) &&
             (pvv.at(2).size() == 1) && (pvv.at(2).at(0)->Flav().IsLepton() == true) &&
             (pvv.at(3).size() == 2) && (pvv.at(3).at(0)->Flav().IsLepton() == true)) 
                                  {
                                    m_at_all = true;
                                    m_ME = ME_List::Tau_To_Lepton_Neutrinos;
                                  }
    // Vector_To_Fermion_Fermion (neutral vector)
    else if ((pvv.size() == 4) &&
             (pvv.at(0).size() == 0) &&
             (pvv.at(1).size() == 1) && (pvv.at(1).at(0)->Flav().IsVector() == true) &&
             (pvv.at(2).size() == 2) && (pvv.at(2).at(0)->Flav().IsFermion() == true) &&
             (pvv.at(3).size() == 0))
                                  {
                                    m_at_all = true;
                                    m_ME = ME_List::Vector_To_Fermion_Fermion;
                                  }
    // Scalar_To_Fermion_Fermion (neutral scalar)
    else if ((pvv.size() == 4) &&
             (pvv.at(0).size() == 0) &&
             (pvv.at(1).size() == 1) && (pvv.at(1).at(0)->Flav().IsScalar() == true) &&
             (pvv.at(2).size() == 2) && (pvv.at(2).at(0)->Flav().IsFermion() == true) &&
             (pvv.at(3).size() == 0))
                                  {
                                    m_at_all = true;
                                    m_ME = ME_List::Scalar_To_Fermion_Fermion;
                                  }
    // Scalar_To_Scalar_Scalar (neutral scalar)
    else if ((pvv.size() == 4) &&
             (pvv.at(0).size() == 0) &&
             (pvv.at(1).size() == 1) && (pvv.at(1).at(0)->Flav().IsScalar() == true) &&
             (pvv.at(2).size() == 2) && (pvv.at(2).at(0)->Flav().IsScalar() == true) &&
             (pvv.at(2).size() == 2) && (pvv.at(2).at(1)->Flav().IsScalar() == true) &&
             (pvv.at(3).size() == 0))
                                  {
                                    m_at_all = true;
                                    m_ME = ME_List::Scalar_To_Scalar_Scalar;
                                  }
    else                          {
                                    m_at_all = false;
                                    m_ME = ME_List::unknown;
                                  }
  }

  /*!
    \file ME_List.H
    \brief contains the classes ME_List and ME_Browser
  */

  /*!
    \class ME_List
    \brief contains an enum for all generic and specific infrared subtracted and squared MEs implemented
  */

  /*!
    \class ME_Browser
    \brief decides which, if any, ME is applicable
  */
  //////////////////////////////////////////////////////////////////////////////////////
  // members of ME_Browser
  //////////////////////////////////////////////////////////////////////////////////////
  /*!
    \var  bool ME_Browser::m_at_all
    \brief true, if there exist a suitable ME
  */

  /*!
    \var ME_List::code ME_Browser::m_ME
    \brief contains the information about which ME to be used
  */

  /*!
    \fn ME_Browser::ME_Browser(Particle_Vector_Vector)
    \brief checks, whether there is a suitable ME implemented and sets m_at_all and m_ME accordingly
  */

  /*!
    \fn ME_Browser::~ME_Browser()
    \brief destructor
  */

  /*!
    \fn bool ME_Browser::GetMEInfo()
    \brief returns true, if there is a suitable ME
  */

  /*!
    \fn ME_List::code ME_Browser::WhichME()
    \brief returns the information about which ME to be used
  */


}


#endif
